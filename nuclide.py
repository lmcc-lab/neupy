import pandas as pd
from neupy.databases import load_databases as ld
from pprint import pprint
import logging
import numpy as np
from typing import List, Union, Callable, Optional
import os
from neupy.config import *
from neupy.module_exceptions.nuclide_exceptions import DecayChainDepthWarning
from scipy.integrate import odeint
from neupy import path


class _Bateman:

    def __init__(self, n: int, lambdai: Union[np.ndarray, list], AZI: str, N10: float=1.0):
        self.n = n
        self.lambdai = lambdai
        self.N10 = N10     
        self.AZI = AZI
        self._debug_record = dict()   
    
    @staticmethod
    def bateman_equation_latex(n):
        for nn in range(n):
            prod_str = ''.join(['\lambda_{{{0}}}'.format(i+1) for i in range(nn)])
            frac_str = ''
            for i in range(nn+1):
                denom_str = ''.join(['(\lambda_{{{0}}}-\lambda_{{{1}}})'.format(j+1, i+1) for j in range(nn+1) if j != i])
                frac_str += '{1}{0} e^{{-\lambda_{{{2}}}t}}{3}'.format(r"\frac{" if nn>0 else "","" if i==0 else "+", i+1, "}}{{{0}}}".format(denom_str) if nn > 0 else "")
            print(r'N_{{{0}}}(t) &= N_1(0)\times{1}\times\left({2}\right)'.format(nn+1, prod_str, frac_str), r'\\')

    def bateman_equation(self, t):
        """
        lambdai is an array of length n, [lambda0, lambda1, ... lambdan], where we are calculating Nn

        Bateman equations are

        dN1(t)/dt = -λ1N1(t)
        dNi(t)/dt = -λiNi(t) + λ(i-1)N(i-1)(t)
        dNk(t)/dt = λ(k-1)N(k-1)(t)
        
        With general solutions
        Nn(t) = N1(0)(Π_{i=1}^{n-1}λi)Σ_{i=1}^{n} ( e^{-λit} / (Π_{j=1, j≠i}^{n} (λj - λi))
        """
        prod = np.prod(self.lambdai[:self.n])
        frac = 0
        for i, L in enumerate(self.lambdai[:self.n+1]):
            denom = np.prod([l-L for j, l in enumerate(self.lambdai[:self.n+1]) if i != j])
            exp_term = np.exp(-L*t)
            if denom == 0:
                if self._debug_record.get('bateman_denom_error', None) is None:
                    self._debug_record = {'bateman_denom_error': []}
                self._debug_record['bateman_denom_error'].append({'nuclide': self.AZI, "denom_vals": [(l,L) for j, l in enumerate(self.lambdai[:self.n+1]) if i != j]})
                continue
            frac += exp_term/denom
        return self.N10 * prod * frac

    def integral_bateman_equation(self, t):
        """
        Taking the integral of the bateman equation just divides the exponential terms by -λi.
        
        ∫Nn(t)dt = N1(0)(Π_{i=1}^{n-1})Σ_{i=1}^{n} ( ∫e^{-λit}dt / (Π_{j=1, j≠i}^{n} (λj - λi))
                 = -N1(0)/λi * (Π_{i=1}^{n-1})Σ_{i=1}^{n} ( e^{-λit} / (-λiΠ_{j=1, j≠i}^{n} (λj - λi))
        """
        prod = np.prod(self.lambdai[:self.n])
        frac = 0
        for i, L in enumerate(self.lambdai[:self.n+1]):
            denom = np.prod([l-L for j, l in enumerate(self.lambdai[:self.n+1]) if i != j])
            if L != 0:
                exp_term = np.exp(-L*t)/-L
            else:
                exp_term = t
            if denom == 0:
                self._debug_record['bateman_denom_error'].append({'nuclide': self.AZI, "denom_vals": [(l,L) for j, l in enumerate(self.lambdai[:self.n+1]) if i != j]})
                continue
            ci = prod/denom
            frac += ci * exp_term + ci/L if L != 0 else ci * exp_term
        return self.N10 * frac
    
    def derivative_bateman_equation(self, t):
        """
        Taking the derivative of the bateman equation just multiplies the exponential terms by -λi.
        
        dNn(t)/dt = N1(0)(Π_{i=1}^{n-1})Σ_{i=1}^{n} ( -λie^{-λit} / (Π_{j=1, j≠i}^{n} (λj - λi))
        """
        prod = np.prod(self.lambdai[:self.n])
        frac = 0
        for i, L in enumerate(self.lambdai[:self.n+1]):
            denom = np.prod([l-L for j, l in enumerate(self.lambdai[:self.n+1]) if i != j])
            exp_term = -L * np.exp(-L*t)
            if denom == 0:
                if self._debug_record.get('bateman_denom_error', None) is None:
                    self._debug_record = {'bateman_denom_error': []}
                self._debug_record['bateman_denom_error'].append({'nuclide': self.AZI, "denom_vals": [(l,L) for j, l in enumerate(self.lambdai[:self.n+1]) if i != j]})
                continue
            frac += exp_term/denom
        return self.N10 * prod * frac


class _BatemanProfile:

    def __init__(self, bateman_object: _Bateman, bateman_method: str) -> None:
        self.function = bateman_object.__getattribute__(bateman_method)
    
    def calculate(self, t):
        return self.function(t)

class _WeightedBatemanProfile(_BatemanProfile):

    def __init__(self, bateman_object: _Bateman, bateman_method: str, intensity: float) -> None:
        super().__init__(bateman_object, bateman_method)
        self.intensity = intensity
    
    def calculate(self, t):
        return self.intensity * self.function(t)



class Nuclide:

    def __init__(self, A: int=None, Z: int=None, level: int = 0, AZI:int=None, nubase: pd.DataFrame=None, fy: dict=None, config_file: dict=None, nubase_config: dict=None, search_in_fy=False,
                 force_level_0 = True, database_path = f"{path}\\databases\\") -> None:
        """
        Search for a nuclide in the Nubase and Fission databases. If no database is provided in arguements
        nubase and fy, then it will load whatever databases it can from the databases folder. 
        
        A and Z must be provided, and optionally level or NIA. If neither of these are provided then ValueError 
        will be raised. 

        # Params

        A: int = None, mass number
        Z: int = None, proton number
        level: int = 0, level=0 (gs); level=1,2 (isomers); level=3,4 (levels); level=5 (resonance); level=8,9 (IAS). level=3,4,5,6 can also indicate isomers (when more than two isomers are presented in a nuclide)
        AZI: int = None, the combined AAAZZZi value for searching
        nubase: pd.DataFrame = None, nubase database
        search_in_fy: bool = False, Search through fission yield databases to match input nuclide to database. Generally we run this nuclide class
            after getting a nuclide from the databases, so we should already have this information without having to research, therefore this is False
            by default. 


        """
        if (A is None or Z is None) and AZI is None:
            raise ValueError("A and Z must be provided or AZI")
        
        if AZI is not None:
            AZI = AZI[:6] + '0' if force_level_0 else AZI
            self.AZI = AZI
            self.A = int(AZI[:3])
            self.Z = int(AZI[3:6])
            self.level = int(AZI[6])
            self.AAA = AZI[:3]
            self.ZZZi = AZI[3:]
        elif Z is not None and A is not None:
            # convert A, Z, level into AZI format.
            self.AZI = self._azi_gen(A, Z, level)
            self.Z = Z
            self.A = A
            self.level = level
        else:
            raise ValueError(f"Z={Z} and A={A}, both of these need to be not None if AZI is None")
        

        # load nubase
        if nubase is None:
            try:
                nubase, config_file = ld.load_nubase(path=database_path, log_level=logging.ERROR)
                nubase_config = config_file['nubase2020']
            except FileNotFoundError:
                print("WARNING: nubase2020.txt not found, trying nubase2016.txt\n")
                try:
                    nubase, config_file = ld.load_nubase(path=database_path, filename='nubase2016.txt', log_level=logging.ERROR)
                    nubase_config = config_file['nubase2016']
                except FileNotFoundError:
                    raise FileNotFoundError("Nuclide requires nubase2020.txt or nubase2016.txt to be present in .databases/ directory. Please check that one of these files exists.")

        self.nubase = nubase
        self.config = config_file
        self.nubase_config = nubase_config

        if fy is None:
            fy = ld.load_all_fy_databases(path=database_path)
        
        self.fy = fy

        # Find AZI in fy
        self.nuclide_fission_info = {}
        if search_in_fy:
            for fy_type, db in self.fy.items():
                self.nuclide_fission_info[fy_type] = db.loc[db['AZI'] == self.AZI]

        # Find AZI in nubase
        self.nuclide_nubase_info = self.nubase.loc[self.nubase['AZI'] == self.AZI]
        self.found = True
        if self.nuclide_nubase_info.shape[0] == 0:
            self.found = False
        
        # Initilise some params
        self.decay_chain = pd.DataFrame([], columns=['nuclide', 'dm'])
        self.allow_energetically_possible_decay_paths = True
        self.br_intensity_extra_params = False
        self.decay_chain_made = False
        self.total_neutrino_profile = 0
        self.integral_total_neutrino_profile = 0
        self.derivative_total_neutrino_profile = 0
        self._focused_linear_chain = None
        self.fission_yeild = None
        self.total_neutrino_profile_with_xe_poisoning = 0


    @staticmethod
    def _azi_gen(A, Z, level=0):
        return ''.join(["0"]*(3-len(f"{A}")))+f"{A}" + ''.join(["0"]*(3-len(f"{Z}"))) + f"{Z}" + f"{level}"

    def all_nubase_data_on_nuclide(self):
        """
        Run this method if you don't want just the exact AZI entry in nubase, but instead want all the information
        given just the A and Z numbers. 

        This updates self.nuclide_nubase_info
        """
        self.nuclide_nubase_info = self.nubase.loc[(self.nubase['A'] == self.A) & (self.nubase['Z'] == self.Z)]

    def convert_intensity_to_decimal(self, intensity):
        if intensity in ['', '?']:
            return 0.0
        if isinstance(intensity, str) and '[' in intensity:
            intensity = intensity.split('[', 1)[0]
            self.br_intensity_extra_params = True
        return float(intensity)/100

    def make_decay_chain(self, include_isomer_decay=False, allow_energtically_possible=True, path_num='', generation=0):
        """
        Decay chains are generated using the decay mode information for the nuclide for each generation of the decay chain. 
        
        This method recursively generates the decay chain, only exiting the recursion if it meets a stable nuclide. 

        ## Params
        include_isomer_decay: bool=False, If the fission product is in a higher energy level, then you can choose whether to include 
                                            it's internal decay as part of the decay chain by setting include_isomer_decay to True, 
                                            default is False.
        allow_energetically_possible: bool=False, Some decay modes are energetically possible but not observed. You can choose to
                                                  allow or disallow these decay paths using this toggle. By default it's set to True.
        path_num: leave as '', it is passed on during the recursion.
        """
        if self.decay_chain_made:
            return

        generation += 1
        if generation == max_chain_depth:
            raise DecayChainDepthWarning(f"Max generations reached for nuclide {self.AZI}")
        
        self.allow_energetically_possible_decay_paths = allow_energtically_possible
        if not include_isomer_decay:                                    # If this flag is False, then set nuclide level to 0 and assume no isomer
            self.level = 0                                              # decay. Will make True case in future
            self.AZI = self._azi_gen(self.A, self.Z, self.level)
        br = self.nuclide_nubase_info.loc[self.AZI, 'BR']               # Get the decay modes from self.nuclide_nubase_info
        half_life = self.nuclide_nubase_info.loc[self.AZI, 'T #']
        paths = br.split(';')                                           # Break them into the different decay modes (seperated by ';')
        decay_mode_key = self.nubase_config['decay_mode_key']           # Get the key from the config file
        es = self.nubase_config['equality_symbols']                     # Get the equality symbols from the config file
        decay_paths = pd.DataFrame([], columns=['nuclide', 'dm', 'intensity', 'dintensity', 'AZI'])       # Initialise the decay_paths dataframe. Dataframe structure is | decay_path (Index) | nuclide | dm |  
        path_num_offset = 0
        for new_path_num, path in enumerate(paths):                     # Go through the decay modes of the nuclide, labeling the path num using new_path_num

            if ('?' in path and not '=?' in path) and not allow_energtically_possible: # If allow_energyetically_possible is set to True, then nothing will happen
                continue                                                               # Otherwise, if False, then it will check if only '?' is present in BR path
                                                                                       # If it is then we skip this one.
            for sym in es:                                              # This section is just to split the decay mode from the branch ratio, where we split based on the
                if sym not in path:                                     # different possible equality symbols used. For example, B-=100 or A~90 etc.
                    continue                                            # 
                dm, ratio = path.split(sym, 1)                          #
                dratio = 0
                if ' ' in ratio:
                    ratio, dratio = ratio.split(' ', 1)
                if half_life == 'stbl':                                 # We also check if the decay mode is IS, meaning the nuclide is stable and it's giving the abundance
                    return                                              # instead. We return out of this function if this is true.
                break                                                   #
            
            transform_key = decay_mode_key.get(dm, {"transform": 'None'})['transform']
            if transform_key == 'None':
                path_num_offset += 1
                continue
            dA, dZ = eval(transform_key)                                # Convert the nuclide transformation from string to tuple, representing (delta A, delta Z) for the given
                                                                        # decay mode
            next_gen = self._azi_gen(self.A + dA, self.Z + dZ)          # Convert the changed nuclide to AZI notation
            next_nuclide = Nuclide(AZI=next_gen, nubase=self.nubase, fy=self.fy, config_file=self.config, nubase_config=self.nubase_config) # Make a new nuclide instance with 
                                                                                                                                            # the next nuclide in the decay chain
                                                                                                                                            # all databases are passed in to reduce
                                                                                                                                            # RAM usage
            df = pd.DataFrame({"nuclide": [next_nuclide], "dm": [dm], 
                               "intensity": [self.convert_intensity_to_decimal(ratio)], 
                               'dintensity': [self.convert_intensity_to_decimal(dratio)],
                               "AZI": next_gen}, 
                               index=[path_num+str(new_path_num - path_num_offset)])  # Update the decay_path dataframe with the next nuclide. Path_num
            decay_paths = pd.concat([decay_paths, df])                                                      # is updated as a string representing the decay path of the decay
                                                                                                            # chain. For example, if path_num is 0, and new_path_num is 1, then
                                                                                                            # the index is 01

        self.decay_chain = pd.concat([self.decay_chain, decay_paths])   # add this decay path to self.decay_chain
                                                                        # Up to now we have found all of the daughter nuclides for a generation, now we will go through each daughter
                                                                        # nuclide and find the next generation, keeping track of the parents using path_num
        for index, path in decay_paths.iterrows():  
            path['nuclide'].make_decay_chain(include_isomer_decay=include_isomer_decay, path_num=index, generation=generation) # Go through the daughter nuclides of this generation and recursively make
            self.decay_chain = pd.concat([self.decay_chain, path['nuclide'].decay_chain])               # the next generation, passing in the parents path_num to keep track of 
                                                                                                        # the decay chain family. When the recursion is done, add this to self.decay_chain

    def break_decay_chain_branches(self) -> None:
        """
        This method seperates the decay chains generated using self.make_decay_chain into the individual 
        branches. This method updates self.decay_chain with the new broken branches.

        The dtype of self.decay_chain is originally a pd.DataFrame. After running this method, it is
        updated as a List[pd.DataFrame]. The dataframe is the linear decay chain, with columns 'nuclide'
        'dm', like that of the original self.decay_chain, and each element of the list is the unique
        decay chain.
        """
        self.decay_chain_made = True
        if isinstance(self.decay_chain, list):
            print("Branches already broken")
            return
        linear_chain = []
        replica_decay_chain = self.decay_chain.copy()               # Make a copy of the original decay chain
        replica_decay_chain['used'] = False                         # Make a new column called 'used' and set all values to False. This will keep track of
                                                                    # whether a nuclide has been used in a decay chain already, to avoid creating sub branches
        for path_num, row in replica_decay_chain[::-1].iterrows():  # Go through decay chain in reverse order.                                       
            # print(replica_decay_chain)
            # print(row)
            # print('___________________')
            if replica_decay_chain.loc[path_num, 'used']:                                            # Check if this row has been used, if it has then skip (it's apart of another decay chain)
               continue                                             # This works without having to go in order of longest chain to shortest because
                                                                    # self.decay_chain is generated sequentially, filling out a decay chain one at a time in
                                                                    # full before moving on.
            this_path = {"nuclide": [row.nuclide], 
                         "dm": [row.dm],
                         "intensity": [row.intensity],
                         "dintensity": [row.dintensity],
                         "path_index": [path_num],
                         "AZI": [row.AZI]}
            flag_used = [path_num]                                  # Record what paths have been used
            for follow_path in range(1, len(path_num)):             # Index through the path_num, excluding the last number. For example, if the path_num is 
                prev_nuclide_path = path_num[:-follow_path]         # 600, then these two lines make the previous path: 60, 6. Therefore, we are going from
                                                                    # the last nuclide in the decay chain and working backwared
                prev_nuclide = replica_decay_chain.loc[prev_nuclide_path]                       # get previous nuclide
                this_path["nuclide"].append(prev_nuclide.nuclide)                               # Update the decay path
                this_path['dm'].append(prev_nuclide.dm)
                this_path['intensity'].append(prev_nuclide.intensity)
                this_path['dintensity'].append(prev_nuclide.dintensity)
                this_path['path_index'].append(prev_nuclide_path)
                this_path['AZI'].append(prev_nuclide.AZI)
                flag_used.append(prev_nuclide_path)                                             # Add used flag            
            replica_decay_chain.loc[flag_used, 'used'] = True       # Update replica_decay_chain used column with used flag
            this_path['nuclide'].reverse()                          # reverse order of this_path so it's in order [Parent, child, grandchild, ...] for a branch
            this_path['dm'].reverse()
            this_path['intensity'].reverse()
            this_path['dintensity'].reverse()
            this_path['path_index'].reverse()
            this_path['AZI'].reverse()
            this_path = pd.DataFrame(this_path)
            linear_chain.append(this_path)                          # Add this path to the linear chain
        self.decay_chain = linear_chain                             # Update self.decay_chain with the linear chains.


    @staticmethod
    def _mermaid_flowchart_template(element: str, connect_to: str, note: str=''):
        return f'{element} --"{note}"--> {connect_to}'
    
    @staticmethod
    def _make_folder(path: str, **kwargs) -> None:
        """Makes a folder given the desired path, only if the path doesn't already exist"""
        cwd = kwargs.get('cwd', os.getcwd())
        if path[-1] == '/':
            path = path[:-1]
        if path not in os.listdir(cwd):
            os.makedirs(path)

    def display_decay_chain(self, save_path: str='decay_chains/', **kwargs) -> List[pd.DataFrame]:
        """
        Take the broken decay chains and return the same pandas dataframe but with the nuclide
        symbol instead of the nuclide object location. This is used for visually checking that
        the decay chain makes sense and is what was expected. This method also saves a mermaid
        flow chart representing this decay chain and saves it by default in 'decay_chains/' 
        with the name '{nuclide symbol} decay chain.md'. This can then be viewed on platforms
        such as github or Obsidian, or any markdown viewer with mermaid available. 

        ## Params
        save_path: str, default is decay_chains/, the path you want to save your mermaid.md
                        diagrams
        **kwargs: passed to _make_folder method
        
        ## Returns
        replica_decay_chain: List[pd.DataFrame], a list of the linear branches an element can
                                                 decay down, with the symbols shown in the
                                                 nuclide column for ease of reading.
        
        """
        if isinstance(self.decay_chain, pd.DataFrame):                      # Break the decay chain branches if not already done
            self.break_decay_chain_branches()
        
        mermaid_template = \
'''
```mermaid 
flowchart TD
'''                                                                 # Excuse the ugly syntax, this is just to get the
                                                                    # right spacing for nice looking .md files. This 
                                                                    # will change in the future to reference a variable
                                                                    # in a different file. 
        record_connections = []                                     # Record connections records what mermaid connections
                                                                    # have already been drawn (stops copies).
        replica_decay_chain = [df.copy() for df in self.decay_chain]# Make copy of self.decay_chain so we don't change 
                                                                    # self.decay_chain.
        for i, dc in enumerate(replica_decay_chain):                # go through each of the linear paths
            for index, row in dc.iterrows():                        # Go through the elements for each generation of a linear path
                cur_element_sym = row.nuclide.nuclide_nubase_info.loc[row.nuclide.AZI, 'A El'] # get the nuclide symbol of the current nuclide
                decay_mode = row['dm']                              # get the decay mode of the current nuclide
                if not index:                                       # If index==0, then this will be a connection from the original nuclide
                    if decay_mode is not None:

                        mermaid_flow = self._mermaid_flowchart_template(self.nuclide_nubase_info.loc[self.AZI, 'A El'], # Make the mermaid flow chart
                                                                    cur_element_sym,                                # for this connection
                                                                        note=decay_mode)
                        if mermaid_flow not in record_connections:      # Only add this to the mermaid template if it isn't already there
                            record_connections.append(mermaid_flow)
                            mermaid_template += \
f'''{mermaid_flow}
'''
                if index == dc.shape[0]-1:                          # Stop if this is the last element in the linear chain
                    break
                next_nuclide = dc.loc[index+1, 'nuclide']           # Otherwise, get the next nuclide in this chain
                this_decay_mode = dc.loc[index+1, 'dm']             # Get the decay mode of this nuclide (remember self.decay_mode is
                                                                    # nuclide, the decay mode that got to this nuclide).
                next_element_sym = next_nuclide.nuclide_nubase_info.loc[next_nuclide.AZI, 'A El']  # Get the next nuclides symbol
                
                mermaid_flow = self._mermaid_flowchart_template(cur_element_sym, next_element_sym, note=this_decay_mode) # Create the mermaid
                                                                    # flow chart
                if mermaid_flow in record_connections:              # Skip adding it if it's already there
                    continue
                record_connections.append(mermaid_flow)             # Add this mermaid flowchart if it hasn't already been used and record
                                                                    # that it's been used
                mermaid_template += \
f'''{mermaid_flow}
''' 
            dc['nuclide'] = dc['nuclide'].apply(lambda row: row.nuclide_nubase_info.loc[row.AZI, 'A El']) # Convert the whole nuclide column
                                                                                                          # to the symbols
        mermaid_template += '```'
        
        self._make_folder(save_path, **kwargs)                                 # make the folder (if needed)
        cwd = kwargs.get('cwd', os.getcwd())
        save_path = '\\'+save_path.replace('/', '\\')
        with open(f"{cwd}{save_path}{self.nuclide_nubase_info.loc[self.AZI, 'A El']} decay chain.md", 'w+') as f: # save the mermaid file
            f.writelines(mermaid_template)
            
        return replica_decay_chain                                  # return the updated readable decay chain  
    
    def convert_half_life(self, half_life, unit):
        if not isinstance(half_life[0], float):
            try:
                half_life = (float(half_life[0]), half_life[1])
            except ValueError:
                return half_life
        prefix = 1
        if len(unit)>1:
            prefix = self.config['prefix'][unit[0]]
        half_life = (half_life[0] * prefix * self.config['unit'][unit[-1]], half_life[1])
        return half_life
    
    @staticmethod
    def _check_half_life_flags(nubase_half_life, allow_theoretical):
        if not isinstance(nubase_half_life, tuple):
            return nubase_half_life, None
        if nubase_half_life[1] == '#' and allow_theoretical:
            return nubase_half_life[0], '#'
        return nubase_half_life

    def get_half_life_and_convert(self, chain_index, allow_theoretical):
        linear_chain = self.decay_chain[chain_index]
        linear_chain['half_life'] = linear_chain['nuclide'].apply(
            lambda row: self.convert_half_life(self._check_half_life_flags(row.nuclide_nubase_info.loc[row.AZI, 'T #'], allow_theoretical),
                                               row.nuclide_nubase_info.loc[row.AZI, 'unit T']))
        linear_chain['half_life_uncert'] = linear_chain['nuclide'].apply(
            lambda row: self.convert_half_life(self._check_half_life_flags(row.nuclide_nubase_info.loc[row.AZI, 'dT'], allow_theoretical),
                                               row.nuclide_nubase_info.loc[row.AZI, 'unit T']))

    @staticmethod
    def convert_half_lives_to_decay_constant(half_life):
        if half_life == '':
            return 0
        return np.log(2)/half_life

    def initialise_bateman_objects(self, allow_theoretical=True, *args, **kwargs):
        """
        calculate the concentration profiles of your decay chains for this nuclide. These will be saved
        in self.decay_chain in a new column called 'concentration_profile'. These concentrations are the
        bateman_equation methods of the Bateman object, where the Bateman objects are is stored in the
        'bateman_object' column.

        Parameters
        ---------------
        *args, passed to self.make_decay_chain
        **kwargs, passed to self.make_decay_chain

        self.make_decay_chain.__docs__ 
        ---------------
        Decay chains are generated using the decay mode information for the nuclide for each generation of the decay chain. 
        
        This method recursively generates the decay modes, only exiting the recursion if it meets a stable nuclide. 

        Parameters
        ---------------
        include_isomer_decay: bool=False, If the fission product is in a higher energy level, then you can choose whether to include 
                                            it's internal decay as part of the decay chain by setting include_isomer_decay to True, 
                                            default is False.
        allow_energetically_possible: bool=False, Some decay modes are energetically possible but not observed. You can choose to
                                                  allow or disallow these decay paths using this toggle. By default it's set to True.
        path_num: leave as '', it is passed on during the recursion.
        """
        
        if isinstance(self.decay_chain, pd.DataFrame):
            if self.decay_chain.shape[0] == 0:
                self.make_decay_chain(*args, **kwargs)
            self.break_decay_chain_branches()
        
        
        for chain_index, linear_chain in enumerate(self.decay_chain):
            self_chain_data = pd.DataFrame({"nuclide": [self], "dm": [None], "intensity": 1.0, "dintensity": 0.0, "path_index": [None], 'AZI': self.AZI})
            self.decay_chain[chain_index] = pd.concat([self_chain_data, linear_chain], ignore_index=True)
            self.get_half_life_and_convert(chain_index, allow_theoretical)
            half_lives = self.decay_chain[chain_index]['half_life'].values
            lambdai = [0 if hl[0] == 'stbl' else self.convert_half_lives_to_decay_constant(hl[0]) for hl in half_lives]
            self.decay_chain[chain_index]['bateman_object'] = self.decay_chain[chain_index].apply(lambda row: _Bateman(row.name, lambdai, self.AZI), axis=1)
    
    @staticmethod
    def concentration_profile(linear_chain: pd.DataFrame):
        return linear_chain.apply(lambda row: _BatemanProfile(row.bateman_object, 'bateman_equation').calculate, axis=1)

    def calculate_profiles(self, t, linear_chain):
        df = self.concentration_profile(linear_chain)
        return df.apply(lambda row: row(t))

    def _map_element_wise_multiplication(self, vec1, vec2):
        calc = 0
        for elem in vec1:
            calc += elem * vec2
        return calc

    def nuclide_concentration_profile(self, *args, **kwargs):
        """
        *args and **kwargs passed to self.make_decay_chain
        """
        if isinstance(self.decay_chain, pd.DataFrame,):
            if self.decay_chain.shape[0] == 0:
                self.make_decay_chain(*args, **kwargs)
            self.break_decay_chain_branches()

        concentration_profile = []
        for linear_chain in self.decay_chain:
            concentration_profile.append(linear_chain.apply(lambda row: _WeightedBatemanProfile(row.bateman_object, 'bateman_equation', row.intensity).calculate, axis=1))
        return concentration_profile
        
    
    def xe_poisoning(self, **neutron_kwargs):
        """
        Xe135 poisoning transmutates Xe135 to Xe136. This can be included in the concentration models
        
        neutron_kwargs are passed to _XeConcentrationModel. There are three kwargs, which are
        neutron_percent: Union[Callable, float] = 1.0, 
        neutron_flux: float = 1.6e9,
        cross_section: float = 2.805e-18
        """
        Xe135AZI = '1350540'
        focused_AZI = ['1350530', Xe135AZI, '1350550', '1350560']
        # go through decay chains in this nuclide, and check if Xe135 is present
        for i, linear_chain in enumerate(self.decay_chain):
            # Check if Xe135 is in linear_chain
            linear_chain['XePoisObject'] = linear_chain['bateman_object'].apply(lambda row: row.bateman_equation)
            linear_chain['XePois'] = linear_chain['XePoisObject'].copy()

            half_lives = linear_chain['half_life'].values

            half_lives = [hl for AZI, hl in zip(linear_chain['AZI'].values, half_lives) if AZI in focused_AZI]

            lambdai = [0 if hl[0] == 'stbl' else self.convert_half_lives_to_decay_constant(hl[0]) for hl in half_lives]
            
            precursor_concentration = None
            precursor_dc = None

            linear_chain_size = linear_chain.shape[0]

            if '1350530' in linear_chain['AZI'].values:
                precursor_concentration = linear_chain.loc[linear_chain_size-4,'bateman_object'].bateman_equation
                precursor_dc = lambdai[0]
                lambdai = lambdai[1:]
            linear_chain.loc[linear_chain_size-3:, 'XePoisObject'] = linear_chain[linear_chain_size-3:].apply(lambda row: _XePoisConc(row.name-linear_chain_size+3, lambdai, precursor_concentration=precursor_concentration, precursor_decay_constant=precursor_dc, **neutron_kwargs), axis=1)
            linear_chain.loc[linear_chain_size-3:, 'XePois'] = linear_chain.loc[linear_chain_size-3:, 'XePoisObject'].apply(lambda row: row.calculate)

            Xe136 = Nuclide(136, 54, nubase=self.nubase, fy=self.fy, nubase_config=self.nubase_config,
                            config_file = self.config)
            
            Xe136_concentration = _XePoisConc(3, lambdai, linear_chain.loc[linear_chain_size-4, 'XePois'], lambdai[-3], **neutron_kwargs)
            
            Xe136_contribution = pd.DataFrame([[Xe136, Xe136.AZI, ('stbl', None), Xe136_concentration.calculate, Xe136_concentration]], columns=['nuclide', 'AZI', 'half_life', 'XePois', 'XePoisObject'], index=[linear_chain.shape[0]])
            
            updated_linear_chain = pd.concat([linear_chain, Xe136_contribution], axis=0)

            self.decay_chain[i] = updated_linear_chain
    
    def weighted_xe_poisoning(self, **neutron_kwargs):
        self.xe_poisoning(**neutron_kwargs)
        for linear_chain in self.decay_chain:
            linear_chain['WeightXePois'] = linear_chain.apply(lambda row: _WeightedBatemanProfile(row.XePoisObject, 'calculate', row.intensity), axis=1)
            

class _XePoisConc:

    def __init__(self, chain_index: int, xe_chain_decay_constants: List[Nuclide], precursor_concentration: Optional[_Bateman]=None, precursor_decay_constant: Optional[Nuclide]=None,
                 neutron_percent: Union[Callable, float]=1.0,
                 neutron_flux: float=1.6e9,
                 neutron_cross_section: float = 2.805e-18
                 ):
        self.prior_concentration = precursor_concentration if precursor_concentration is not None else lambda t: 0
        self.xe_creation_rate = precursor_decay_constant if precursor_decay_constant is not None else 0
        self.xe_dc, self.cs_dc, self.ba_dc = xe_chain_decay_constants
        self.n = chain_index
        if callable(neutron_percent):
            self.neutron_absorption_term = lambda t: neutron_percent(t) * neutron_flux * neutron_cross_section
        else:
            self.neutron_absorption_term = lambda t: neutron_percent * neutron_flux * neutron_cross_section
        self.N0 = [0, 0, 0, 0] if precursor_concentration is not None else [1, 0, 0, 0]
    
    def model(self, N, t):
        Nx, NCs, NBa, Nx136 = N
        dNx = self.xe_creation_rate * self.prior_concentration(t)  - (self.xe_dc + self.neutron_absorption_term(t)) * Nx
        dNCs = self.xe_dc * Nx - self.cs_dc * NCs
        dNBa = self.cs_dc * NCs
        dNx136 = self.neutron_absorption_term(t) * Nx
        return [dNx, dNCs, dNBa, dNx136]
    
    def calculate(self, t):
        N = odeint(self.model, self.N0, t)[:, self.n]
        return N



if __name__ == "__main__":
    nuclide = Nuclide(AZI='0910420')

    nuclide.make_decay_chain()
    nuclide.display_decay_chain()

    t = np.logspace(0, 16, 100)
    # nuclide.bateman_equation_latex(2)
    # nuclide.concentration_profiles(t)
    # print(nuclide.decay_chain)
    # # print(nuclide.display_decay_chain())
    # import matplotlib.pyplot as plt
    # for linear_chain in nuclide.decay_chain:
    #     concentrations = linear_chain['concentration_profile']
    #     plt.stackplot(t, *concentrations, labels=[nuc.nuclide_nubase_info.loc[nuc.AZI, 'A El'] for nuc in linear_chain['nuclide'].values])
    # plt.xscale('log')
    # plt.legend()
    # plt.show()