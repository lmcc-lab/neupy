import pandas as pd
from databases import load_databases as ld
from pprint import pprint
import logging
import numpy as np
from typing import List
import os


class Nuclide:

    def __init__(self, A: int=None, Z: int=None, level: int = 0, AZI:int=None, nubase: pd.DataFrame=None, fy: dict=None, config_file: dict=None, nubase_config: dict=None, search_in_fy=False) -> None:
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
                nubase, config_file = ld.load_nubase(log_level=logging.ERROR)
                nubase_config = config_file['nubase2020']
            except FileNotFoundError:
                print("WARNING: nubase2020.txt not found, trying nubase2016.txt\n")
                try:
                    nubase, config_file = ld.load_nubase(filename='nubase2016.txt', log_level=logging.ERROR)
                    nubase_config = config_file['nubase2016']
                except FileNotFoundError:
                    raise FileNotFoundError("Nuclide requires nubase2020.txt or nubase2016.txt to be present in .databases/ directory. Please check that one of these files exists.")

        self.nubase = nubase
        self.config = config_file
        self.nubase_config = nubase_config

        if fy is None:
            fy = ld.load_all_fy_databases()
        
        self.fy = fy

        # Find AZI in fy
        self.nuclide_fission_info = {}
        if search_in_fy:
            for fy_type, db in self.fy.items():
                self.nuclide_fission_info[fy_type] = db.loc[db['AZI'] == self.AZI]

        # Find AZI in nubase
        self.nuclide_nubase_info = self.nubase.loc[self.nubase['AZI'] == self.AZI]

        # Initilise some params
        self.decay_chain = pd.DataFrame([], columns=['nuclide', 'dm'])

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

    def make_decay_chain(self, include_isomer_decay=False, path_num=''):
        """
        Decay chains are generated using the decay mode information for the nuclide for each generation of the decay chain. 
        
        If the fission product is in a higher energy level, then you can choose whether to include it's internal decay as part of the decay chain
        by setting include_isomer_decay to True, defaul is False.

        This method recursively generates the decay modes, only exiting the recursion if it meets a stable nuclide. 
        """
        if not include_isomer_decay:                                    # If this flag is False, then set nuclide level to 0 and assume no isomer
            self.level = 0                                              # decay. Will make True case in future
            self.AZI = self._azi_gen(self.A, self.Z, self.level)
        br = self.nuclide_nubase_info.loc[self.AZI, 'BR']               # Get the decay modes from self.nuclide_nubase_info
        half_life = self.nuclide_nubase_info.loc[self.AZI, 'T #']
        paths = br.split(';')                                           # Break them into the different decay modes (seperated by ';')
        decay_mode_key = self.nubase_config['decay_mode_key']           # Get the key from the config file
        es = self.nubase_config['equality_symbols']                     # Get the equality symbols from the config file
        decay_paths = pd.DataFrame([], columns=['nuclide', 'dm'])       # Initialise the decay_paths dataframe. Dataframe structure is | decay_path (Index) | nuclide | dm |  
        path_num_offset = 0
        for new_path_num, path in enumerate(paths):                     # Go through the decay modes of the nuclide, labeling the path num using new_path_num

            for sym in es:                                              # This section is just to split the decay mode from the branch ratio, where we split based on the
                if sym not in path:                                     # different possible equality symbols used. For example, B-=100 or A~90 etc.
                    continue                                            # 
                dm, ratio = path.split(sym, 1)                          #
                if half_life == 'stbl':                                 # We also check if the decay mode is IS, meaning the nuclide is stable and it's giving the abundance
                    return                                              # instead. We return out of this function if this is true.
                break                                                   #
            
            transform_key = decay_mode_key[dm]['transform']
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
            df = pd.DataFrame({"nuclide": [next_nuclide], "dm": [dm]}, index=[path_num+str(new_path_num - path_num_offset)])  # Update the decay_path dataframe with the next nuclide. Path_num
            decay_paths = pd.concat([decay_paths, df])                                                      # is updated as a string representing the decay path of the decay
                                                                                                            # chain. For example, if path_num is 0, and new_path_num is 1, then
                                                                                                            # the index is 01

        self.decay_chain = pd.concat([self.decay_chain, decay_paths])   # add this decay path to self.decay_chain
                                                                        # Up to now we have found all of the daughter nuclides for a generation, now we will go through each daughter
                                                                        # nuclide and find the next generation, keeping track of the parents using path_num
        for index, path in decay_paths.iterrows():  
            path['nuclide'].make_decay_chain(include_isomer_decay=include_isomer_decay, path_num=index) # Go through the daughter nuclides of this generation and recursively make
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
                         "path_index": [path_num]}
            flag_used = [path_num]                                  # Record what paths have been used
            for follow_path in range(1, len(path_num)):             # Index through the path_num, excluding the last number. For example, if the path_num is 
                prev_nuclide_path = path_num[:-follow_path]         # 600, then these two lines make the previous path: 60, 6. Therefore, we are going from
                                                                    # the last nuclide in the decay chain and working backwared
                prev_nuclide = replica_decay_chain.loc[prev_nuclide_path]                       # get previous nuclide
                this_path["nuclide"].append(prev_nuclide.nuclide)                               # Update the decay path
                this_path['dm'].append(prev_nuclide.dm)
                this_path['path_index'].append(prev_nuclide_path)
                flag_used.append(prev_nuclide_path)                                             # Add used flag
            replica_decay_chain.loc[flag_used, 'used'] = True       # Update replica_decay_chain used column with used flag
            this_path['nuclide'].reverse()                          # reverse order of this_path so it's in order [Parent, child, grandchild, ...] for a branch
            this_path['dm'].reverse()
            this_path['path_index'].reverse()
            this_path = pd.DataFrame(this_path)
            linear_chain.append(this_path)                          # Add this path to the linear chain
        self.decay_chain = linear_chain                             # Update self.decay_chain with the linear chains.


    @staticmethod
    def mermaid_flowchart_template(element: str, connect_to: str, note: str=''):
        return f'{element} --"{note}"--> {connect_to}'
    
    @staticmethod
    def _make_folder(path: str) -> None:
        """Makes a folder given the desired path, only if the path doesn't already exist"""
        if path[-1] == '/':
            path = path[:-2]
        if path not in os.listdir():
            os.makedirs(path)

    def display_decay_chain(self, save_path: str='decay_chains/') -> List[pd.DataFrame]:
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
                    mermaid_flow = self.mermaid_flowchart_template(self.nuclide_nubase_info.loc[self.AZI, 'A El'], # Make the mermaid flow chart
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
                
                mermaid_flow = self.mermaid_flowchart_template(cur_element_sym, next_element_sym, note=this_decay_mode) # Create the mermaid
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
        
        self._make_folder(save_path)                                 # make the folder (if needed)
        with open(f"{save_path}{self.nuclide_nubase_info.loc[self.AZI, 'A El']} decay chain.md", 'w+') as f: # save the mermaid file
            f.writelines(mermaid_template)
            
        return replica_decay_chain                                  # return the updated readable decay chain  
        
