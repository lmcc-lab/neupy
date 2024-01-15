from neupy.nuclide import Nuclide, pprint, pd, np, DecayChainDepthWarning, max_chain_depth, _Bateman
from neupy.databases.load_databases import load_all_fy_databases
from tqdm import tqdm
from typing import Union, Tuple, Dict, List
from neupy.config import *
import pickle
from neupy import path, sublist


class _CumulativeNeutrinos:

    def __init__(self, neut_scale, bateman_object: _Bateman, bateman_method: str) -> None:
        self.neut_scale = neut_scale
        self.function = bateman_object.__getattribute__(bateman_method)

    def _map_element_wise_multiplication(self, vec1, vec2):
        row = vec1.flatten()
        calc = np.zeros(len(row)).tolist()
        for i, elem in enumerate(row):
            calc[i] = elem * vec2
        calc = np.array(calc)
        calc = calc.reshape(*vec1.shape, len(vec2))
        return calc
    
    def calculate(self, t):
        return self._map_element_wise_multiplication(self.neut_scale, self.function(t))
    

class _WeightedCumNeutrinos:

    def __init__(self, fission_yield: float, total_branch_ratio_intensity: float, cum_neutrinos: _CumulativeNeutrinos):
        self.fission_yield = fission_yield
        self.branch_ratio_intensity = total_branch_ratio_intensity
        self.cum_neutrinos = cum_neutrinos
    
    def calculate(self, t):
        return self.fission_yield * self.branch_ratio_intensity * self.cum_neutrinos(t)

class _TotalCumNeutrinos:

    def __init__(self, nuclide: Nuclide, linear_chain_column: str):
        self.nuclide = nuclide
        self.column = linear_chain_column
    
    def calculate(self, t):
        total_neutrino_profile = 0
        for linear_chain in self.nuclide.decay_chain:
            total_neutrino_profile += np.sum(np.vstack([cum_function(t) for i, cum_function in enumerate(linear_chain[self.column].values) if linear_chain['AZI'].values[i] != '1360540']), axis=0)
        return total_neutrino_profile

class _FissionProductTotalNeutrinos:

    def __init__(self, fission_nuclide_array: List[Nuclide], nuclide_attribute: str) -> None:
        self.fission_nuclide_array = fission_nuclide_array
        self.attribute = nuclide_attribute
    
    def calculate(self, t):
        total_neutrinos = 0
        for nuclide in self.fission_nuclide_array:
            total_neutrinos += nuclide.__getattribute__(self.attribute)(t)
        return total_neutrinos
            

class NeutrinoEmission:

    def __init__(self) -> None:
        self.neutrino_table = {"B-": (1, -1, 'e'),
                               "B+": (1, 1, 'e'),
                               "2B-": (2, -1, 'e'),
                               "2B+": (2, 1, 'e'),
                               "e+": (1, 1, 'e'),
                               "EC": (1, 1, 'e'),
                               "EC+B+": (1, 1, 'e'),
                               "B-A": (1, -1, 'e'),
                               "B-n": (1, -1, 'e'),
                               "B-2n": (1, -1, 'e'),
                               "B-3n": (1, -1, 'e'),
                               "B-4n": (1, -1, 'e'),
                               "B-p": (1, -1, 'e'),
                               "B-d": (1, -1, 'e'),
                               "B-t": (1, -1, 'e'),
                               "B-SF": (1, -1, 'e'),
                               "B+p": (1, 1, 'e'),
                               "B+2p": (1, 1, 'e'),
                               "B+3p": (1, 1, 'e'),
                               "B+A": (1, 1, 'e'),
                               "B+pA": (1, 1, 'e'),
                               "B+SF": (1, 1, 'e')}
        self.neutrino_table_key = ("number of neutrinos", 
                                   "matter(antimatter) +(-)",
                                   "flavour (e, mu, tau)")

    def neutrino_vec(self, linear_chain: pd.DataFrame, ignore_matter_type=True, ignore_flavour=True, cumulative=True):
        shape = [linear_chain.shape[0]]
        if not ignore_matter_type:
            shape.append(2)
        
        if not ignore_flavour:
            shape.append(3)
        
        neut_vec = np.zeros(shape)

        dm = linear_chain.dm.values

        if dm[0] is None and len(dm)>1:
            dm = dm[1:]

        neut_data = [self.neutrino_table.get(d, (0, 0, None)) for d in dm]

        for i, neutrino_emission in enumerate(neut_data):
            if ignore_matter_type and ignore_flavour:
                neut_vec[i+1] = neutrino_emission[0]
                continue

            if not ignore_matter_type and ignore_flavour:
                neut_vec[i+1, int((neutrino_emission[1]+1)/2)] = neutrino_emission[0]
                continue

            if not ignore_matter_type and not ignore_flavour:
                neut_vec[i+1, int((neutrino_emission[1]+1)/2), ['e', 'mu', 'tau'].index(neutrino_emission[2])] = neutrino_emission[0]
                continue

            if ignore_matter_type and not ignore_flavour:
                neut_vec[i+1, ['e', 'mu', 'tau'].index(neutrino_emission[2])] = neutrino_emission[0]
                continue
        
        if cumulative:
            neut_vec = neut_vec.cumsum(0)
        return neut_vec
    
    def cumulative_neutrino_profile(self, linear_chain: pd.DataFrame, **neut_vec_kwargs):
        if 'bateman_object' not in linear_chain.columns:
            raise ValueError("bateman_object must have been calculated before running cumulative_neutrino_profile")
        neut_vec_kwargs.pop('cumulative', None)
        neut_vec = self.neutrino_vec(linear_chain, cumulative=True, **neut_vec_kwargs)
        linear_chain['cum_neut'] = linear_chain.apply(
            lambda row: _CumulativeNeutrinos(neut_vec[row.name], row.bateman_object, 'bateman_equation').calculate, axis=1)
        linear_chain['int_cum_neut'] = linear_chain.apply(
            lambda row: _CumulativeNeutrinos(neut_vec[row.name], row.bateman_object, 'integral_bateman_equation').calculate, axis=1)
        linear_chain['der_cum_neut'] = linear_chain.apply(
            lambda row: _CumulativeNeutrinos(neut_vec[row.name], row.bateman_object, 'derivative_bateman_equation').calculate, axis=1)
            
    def check_neutrino_profile_type(self, linear_chain):
        cum_neutrino_test = (linear_chain.shape[0], *linear_chain.loc[0, 'cum_neut'].shape)
        if len(cum_neutrino_test) == 2:
            print("Cumulative neutrinos were calculated without seperating neutrino flavour or matter-antimatter\n"
                  f"Slice index 0: Nuclide in decay chain, with {cum_neutrino_test[0]} nuclides"
                  f"Slice index 1: Time, with length {cum_neutrino_test[1]}")
        elif len(cum_neutrino_test) == 3 and cum_neutrino_test[1] == 2:
            print("Cummulative neutrinos were calculated with seperation of matter-antimatter, but not seperating flavour\n"
                  f"Slice index 0: Nuclide in decay chain, with {cum_neutrino_test[0]} nuclides\n"
                  "Slice index 1: Matter(antimatter) seperation index 0(1)\n"
                  f"Slide index 2: Time, with length {cum_neutrino_test[2]}")
        elif len(cum_neutrino_test) == 3 and cum_neutrino_test[1] == 3:
            print("Cummulative neutrinos were calculated with seperation of flavour, but not matter-antimatter\n"
                  f"Slice index 0: Nuclide in decay chain, with {cum_neutrino_test[0]} nuclides\n"
                  "Slice index 1: flavour, e, mu, tau seperation index with index 0, 1, 2 respectively\n"
                  f"Slice index 2: Time, with length {cum_neutrino_test[2]}")
        elif len(cum_neutrino_test) == 4:
            print("Cummulative neutrinos were calculated with seperation of both flavour and matter-antimatter\n"
                  f"Slice index 0: Nuclide in decay chain, with {cum_neutrino_test[0]} nuclides\n"
                  "Slice index 1: Matter(antimatter) seperation index 0(1)\n"
                  "Slice index 2: flavour, e, mu, tau seperation index with index 0, 1, 2 respectively\n",
                  f"Slice index 3: Time, with length {cum_neutrino_test[3]}")

    
class Neupy(NeutrinoEmission):
    
    def __init__(self, **fy_kwargs) -> None:
        """
        Neupy is a child class of Neutrino emission, which adds extra information about a decay chains
        neutrino profile, specifically the matter and flavour type. 

        Parameters
        -------------------
        - fy_kwargs, passed to databases.load_databases.load_all_fy_databases

        Docs from load_all_fy_databases
        -------------------
        Load any txt file of the form fy{}.txt file, add them to a dictionary with key {}
        and add a column of AZI, assuming all fy files are in the same format (sep = "   ", header=0 with
        columns "Z", "A", "Level", "YI", "YI uncert")

        @params
        path: str, path where fy{}.txt files is contained.

        @returns
        fiss_data: Dict[str: pd.DataFrame]
        """
        super().__init__()
        self.fy = load_all_fy_databases(**fy_kwargs)
        self.nuclide_template = Nuclide(1, 1, database_path=fy_kwargs.get('path', path+'\\databases\\'), search_in_fy=False)
        self.missing_nuclides = []
        self.no_contribution = []
        self._fission_induced_neutrinos_cache = None
        
    def weighted_cumulative_neutrinos(self, linear_chain: pd.DataFrame, fission_yield: float):
        total_branch_ratio_intensity = np.prod(linear_chain['intensity'])
        linear_chain['weight_cum_neut'] = linear_chain.apply(lambda row: _WeightedCumNeutrinos(fission_yield, total_branch_ratio_intensity, row['cum_neut']).calculate, axis=1)
        linear_chain['weight_int_cum_neut'] = linear_chain.apply(lambda row: _WeightedCumNeutrinos(fission_yield, total_branch_ratio_intensity, row['int_cum_neut']).calculate, axis=1)
        linear_chain['weight_der_cum_neut'] = linear_chain.apply(lambda row: _WeightedCumNeutrinos(fission_yield, total_branch_ratio_intensity, row['der_cum_neut']).calculate, axis=1)

    def total_weighted_neutrinos(self, nuclide: Nuclide):
        nuclide.total_neutrino_profile = _TotalCumNeutrinos(nuclide, 'weight_cum_neut').calculate
        nuclide.integral_total_neutrino_profile = _TotalCumNeutrinos(nuclide, 'weight_int_cum_neut').calculate
        nuclide.derivative_total_neutrino_profile = _TotalCumNeutrinos(nuclide, 'weight_der_cum_neut').calculate
            
    def fission_induced_neutrinos(self, nuclide: Nuclide, fy: float, **flag_kwargs):
        nuclide.initialise_bateman_objects(**flag_kwargs)
        for linear_chain in nuclide.decay_chain:
            self.cumulative_neutrino_profile(linear_chain, **flag_kwargs)
            self.weighted_cumulative_neutrinos(linear_chain, fy)
        self.total_weighted_neutrinos(nuclide)
    
    def all_fission_induced_neutrinos(self, 
                                      neutron_energy_range: Union[Tuple[float], float, None] = None, 
                                      use_elements: Union[str, List[str]] = ['U235', 'PU239', 'U233'],
                                      load_saved=False, 
                                      save = True,
                                      file_name = 'neutrino_results',
                                      file_path = f'',
                                      force_level_0 = True,
                                      **kwargs) -> dict:
        """
        
        Parameters
        -------------
        neutron_energy_range: Tuple[float] | float | None, if None it will include all database data without
        consideration of the neutron energy, if it's a float then it will only include databases where the
        neutron energy is what you gave it, otherwise it expects a tuple representing the start and end 
        of the neutron range in question (start, end) inclusive of both ends.

        use_elements: str | List[str], default is ['U235', 'PU239', 'U233']. This specifies which fission elements
        you want to focus on.

        load_saved: bool, default is False. If True, it loads the file under {file_path}{file_name} to reduce
        time spent recalculating everything.

        save: bool, default is True. If True, it saves the results in a file under {file_path}{file_name}.

        file_name: str, default is 'neutrino_results'.

        file_path: str, default is 'databases/' 

        force_level_0: bool, default is True. This is an argument passed to the Nuclide class, which forces an isomer of a nuclide to the ground state.

        **kwargs: passed on to fission_induced_neutrinos
        
        """
        if load_saved:
            print(f'loading fission data from {file_path}{file_name}.pickle')
            with open(f'{file_path}{file_name}.pickle', 'rb') as f:
                fission_product_neutrino_data = pickle.load(f)
            print('Fission data loaded')
            self._fission_induced_neutrinos_cache = fission_product_neutrino_data
            return fission_product_neutrino_data
        
        if isinstance(use_elements, str):
            use_elements = [use_elements]

        fission_product_neutrino_data = dict()
        for fission_element, ne_dict in self.fy.items():

            if fission_element not in use_elements:
                continue
            
            element_fpnd = dict()
            for ne, db in ne_dict.items():
                
                if ne != 'independent':
                    if isinstance(neutron_energy_range, float) and ne != neutron_energy_range:
                        continue
                    
                    if isinstance(neutron_energy_range, tuple) and not neutron_energy_range[0] <= ne <= neutron_energy_range[1]:
                        continue
                
                nuclide_neutrino_data = []
                for AZI, row in tqdm(db.iterrows(), total=db.shape[0], desc=f"{fission_element} {ne}"):
                    if row.Y == 0.0:
                        self.no_contribution.append(AZI)
                        continue
                    nuclide = Nuclide(AZI=AZI, nubase=self.nuclide_template.nubase, fy=self.fy, config_file=self.nuclide_template.config, nubase_config=self.nuclide_template.nubase_config, force_level_0=force_level_0)

                    if not nuclide.found:
                        self.missing_nuclides.append(AZI)
                        continue
                    try:
                        self.fission_induced_neutrinos(nuclide, row.Y, **kwargs)
                    except DecayChainDepthWarning:
                        print(f"Nuclide {AZI}'s decay chain depth exceeded maximum {max_chain_depth}")
                        continue
                    
                    nuclide_neutrino_data.append(nuclide)
                element_fpnd[ne] = {"total_neutrinos": _FissionProductTotalNeutrinos(nuclide_neutrino_data, 'total_neutrino_profile').calculate, 
                                    "int_total_neutrinos": _FissionProductTotalNeutrinos(nuclide_neutrino_data, 'integral_total_neutrino_profile').calculate,
                                    "der_total_neutrinos": _FissionProductTotalNeutrinos(nuclide_neutrino_data, 'derivative_total_neutrino_profile').calculate,
                                    'nuclide_specific': nuclide_neutrino_data}
            fission_product_neutrino_data[fission_element] = element_fpnd
        
        if save:
            with open(f'{file_path}{file_name}.pickle', 'wb+') as f:
                pickle.dump(fission_product_neutrino_data, f)

        self._fission_induced_neutrinos_cache = fission_product_neutrino_data

        return fission_product_neutrino_data
    
    def convert_thermal_to_burn_rate(self, thermal_power: np.ndarray) -> dict:
        # Basic thermal model, where U235 makes up 94% of the thermal power generation
        return {'U235': 0.94 * thermal_power / u235_energy_per_fission_MeV}
    
    def nuclide_cummulative_neutrinos_from_burn_rate(self, nuclide: Nuclide, burn_rate: Union[float, np.ndarray], 
                                                     time: Union[list, np.ndarray]):
        neutrinos = nuclide.total_neutrino_profile(time)
        dt = time[1] - time[0]
        cum_neutrinos = 0
        for i,_ in enumerate(time):
            neut_array = np.zeros(len(time))
            neut_array[i:] = neutrinos[i:]
            cum_neutrinos += burn_rate * neut_array * dt
        return cum_neutrinos

    def cummulative_neutrinos_from_burn_rate(self, burn_rate: dict, time: Union[list, np.ndarray],
                                             assumed_on_prior_s: float=0, 
                                             **fiss_induced_neut_kwargs):
        fy_neutrinos = self.all_fission_induced_neutrinos(**fiss_induced_neut_kwargs)
        dt = time[1] - time[0]

        if assumed_on_prior_s > 0:
            extended_time = np.linspace(0, assumed_on_prior_s, int(assumed_on_prior_s/dt))
            time = np.append(extended_time, time + assumed_on_prior_s)
            print("Extended time")

        cummulative_neutrinos = dict()
        for element, burn_profile in burn_rate.items():
            element_neutrinos = fy_neutrinos.get(element, None)

            if element_neutrinos == None:
                continue
            
            cum_neut_for_elem = dict()
            
            if assumed_on_prior_s > 0:
                extended_burn_rate = np.linspace(burn_profile[0], burn_profile[0], len(extended_time))
                burn_profile = np.append(extended_burn_rate, burn_profile)

            for ne, neutrino_data in element_neutrinos.items():
                print("Calculating cumulative neutrinos over time span")
                total_neutrinos = neutrino_data['total_neutrinos'](time)
                cum_neutrinos = 0
                print("Done")
                for i, instant_burn_rate in enumerate(tqdm(burn_profile, total=len(time), desc=f"Cum neut {element} {ne}")):
                    neut_array = np.zeros(len(time))
                    neut_array[i:] = total_neutrinos[i:]
                    cum_neutrinos += instant_burn_rate * neut_array * dt
                cum_neut_for_elem[ne] = cum_neutrinos[int(assumed_on_prior_s/dt):]

            cummulative_neutrinos[element] = cum_neut_for_elem

        return cummulative_neutrinos
    
    def slope(self, array, dx):
        """
        Find slope at each point, using central difference method for index's i, 0 < i < N, with forward
        and backward difference method for the ends.
        """
        central_diff = (array[2:]-array[:-2])/(2 * dx)
        start = (array[1] - array[0])/dx
        end = (array[-1] - array[-2])/dx
        slope = np.append(start, central_diff)
        slope = np.append(slope, end)
        return slope

    def neutrino_emission_rate(self, burn_rate: Union[dict, None]=None,
                               time: Union[list, np.ndarray, None]=None,
                               assumed_on_prior_s: float = 0,
                               cum_neutrinos: Union[None, Dict[str, Dict[str, np.ndarray]]]=None,
                               dt: float = 1.0,
                               **fiss_induced_neut_kwargs) -> Union[Dict[str, Dict[str, np.ndarray]], Tuple[Dict, float]]:
        
        if (burn_rate is None and time is None) and cum_neutrinos is None:
            raise ValueError("Burn rate and time have to be provided, or cum_neutrinos")
        
        if any([check is None for check in [burn_rate, time]]) and cum_neutrinos is None:
            raise ValueError("Must provide both burn_rate and time if cum_neutrinos are not provided.")

        if cum_neutrinos is None:
            cum_neutrinos = self.cummulative_neutrinos_from_burn_rate(burn_rate, time, assumed_on_prior_s=assumed_on_prior_s, 
                                                                    **fiss_induced_neut_kwargs)
            dt = time[1] - time[0]
        
        emission_rate = dict()
        for element, ne_dict in cum_neutrinos.items():
            element_emission_rate = dict()
            for ne, neutrino_data in ne_dict.items():
                element_emission_rate[ne] = np.gradient(neutrino_data, dt)
            emission_rate[element] = element_emission_rate
        
        return emission_rate, (cum_neutrinos, dt)
            
    def neutrino_flux_for_moving_observer(self, neutrino_emission_rate: Dict[str, Dict[str, np.ndarray]], 
                                          dt=1.0, 
                                          radial_velocity=0.0, 
                                          r0=0.0,
                                          r=None):

        ne_flux = dict()
        for element, ne_dict in neutrino_emission_rate.items():
            fuel_type_specific = dict()
            for ne, neutrino_data in ne_dict.items():
                # Calculate the radial position over time for this observer
                if r is None:
                    time = np.linspace(0, len(neutrino_data) * dt, len(neutrino_data))
                    r = r0 + radial_velocity * time
                flux = neutrino_data / (4 * np.pi * r**2)
                fuel_type_specific[ne] = flux
            
            ne_flux[element] = fuel_type_specific
        return ne_flux
    
    def cumulative_neutrino_profile_xe_pois(self, linear_chain: pd.DataFrame, **neut_vec_kwargs):
        neut_vec_kwargs.pop('cumulative', None)
        neut_vec = self.neutrino_vec(linear_chain, cumulative=True, **neut_vec_kwargs)
        linear_chain['cum_neut_xe_pois'] = linear_chain.apply(
            lambda row: _CumulativeNeutrinos(neut_vec[row.name], row.XePoisObject, 'calculate' if not isinstance(row.XePoisObject, _Bateman) else 'bateman_equation').calculate, axis=1)

    def weighted_cumulative_neutrinos_xe_pois(self, linear_chain: pd.DataFrame, fission_yield: float):
        total_branch_ratio_intensity = np.prod(linear_chain['intensity'])
        linear_chain['weight_cum_neut_xe_pois'] = linear_chain.apply(lambda row: _WeightedCumNeutrinos(fission_yield, total_branch_ratio_intensity, row['cum_neut_xe_pois']).calculate, axis=1)
    
    def total_weighted_neutrinos_xe_pois(self, nuclide: Nuclide):
        nuclide.total_neutrino_profile_with_xe_poisoning = _TotalCumNeutrinos(nuclide, 'weight_cum_neut_xe_pois').calculate

    def fission_induced_neutrinos_xe_pois(self, nuclide: Nuclide, fy: float, **kwargs):
        flag_kwargs = {'ignore_matter_type': kwargs.get('ignore_matter_type', True),
                       'ignore_flavour': kwargs.get('ignore_flavour', True),
                       'cumulative': kwargs.get('cumulative', True)}
        for key in flag_kwargs:
            kwargs.pop(key, None)
        nuclide.xe_poisoning(**kwargs)
        for linear_chain in nuclide.decay_chain:
            self.cumulative_neutrino_profile_xe_pois(linear_chain, **flag_kwargs)
            self.weighted_cumulative_neutrinos_xe_pois(linear_chain, fy)
        self.total_weighted_neutrinos_xe_pois(nuclide)

    def all_fission_induced_neutrinos_xe_pois(self, **kwargs) -> dict:
        """
        
        Parameters
        -------------
        neutron_energy_range: Tuple[float] | float | None, if None it will include all database data without
        consideration of the neutron energy, if it's a float then it will only include databases where the
        neutron energy is what you gave it, otherwise it expects a tuple representing the start and end 
        of the neutron range in question (start, end) inclusive of both ends.

        use_elements: str | List[str], default is ['U235', 'PU239', 'U233']. This specifies which fission elements
        you want to focus on.

        load_saved: bool, default is False. If True, it loads the file under {file_path}{file_name} to reduce
        time spent recalculating everything.

        save: bool, default is True. If True, it saves the results in a file under {file_path}{file_name}.

        file_name: str, default is 'neutrino_results'.

        file_path: str, default is 'databases/' 

        **kwargs: passed on to fission_induced_neutrinos_xe_pois
        
        """
        fission_product_neutrino_data = {}
        for fission_element, fiss_dict in self._fission_induced_neutrinos_cache.items():
            element_fpnd = {}
            for ne, ne_dict in fiss_dict.items():
                nuclide_neutrino_data = []
                for nuclide in tqdm(ne_dict['nuclide_specific'], total=len(ne_dict['nuclide_specific']), desc=f"{fission_element} {ne}"):
                    fy = nuclide.nuclide_fission_info[fission_element][ne]
                    self.fission_induced_neutrinos_xe_pois(nuclide, fy, **kwargs)
                    nuclide_neutrino_data.append(nuclide)
                element_fpnd[ne] = {"total_neutrinos": _FissionProductTotalNeutrinos(nuclide_neutrino_data, 'total_neutrino_profile_with_xe_poisoning').calculate,
                                    'nuclide_specific': nuclide_neutrino_data}
            fission_product_neutrino_data[fission_element] = element_fpnd
        return fission_product_neutrino_data


                    

        
        

        


if __name__ == "__main__":

    neupy = Neupy()
    time = np.linspace(0, 1000, 1000)
    # pprint(neupy.fy)
    thermal_power = np.zeros(len(time))
    thermal_power[100:200] = 20e6  # 20 MW
    burn_rate = neupy.convert_thermal_to_burn_rate(thermal_power)
    cum_neutrinos = neupy.neutrino_emission_rate(burn_rate, time)
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots()
    for i, element in enumerate(cum_neutrinos):
        for j, ne in enumerate(cum_neutrinos[element]):
            ax.plot(time, cum_neutrinos[element][ne], label=element)
    
    axtwin = ax.twinx()
    axtwin.plot(time, thermal_power, color='r')
    # plt.legend()
    # plt.xscale('log')
    plt.show()

    # import matplotlib.pyplot as plt
    # plt.stackplot(time, *[nuc.total_neutrino_profile if not isinstance(nuc.total_neutrino_profile, int) else [0]*len(time) for nuc in neupy.fiss_product_neutrino_data['u235thermal']['nuclide_specific']],
    #               labels=[nuc.nuclide_nubase_info.loc[nuc.AZI, 'A El'] for nuc in neupy.fiss_product_neutrino_data['u235thermal']['nuclide_specific']])

    # # nuclide = Nuclide(135, 52)
    # # neupy.fission_induced_neutrinos(nuclide, neupy.fy['u235thermal'].loc[nuclide.AZI, 'YI'], time)

    
    # # plt.plot(time, neupy.fiss_product_neutrino_data['u235thermal']['total_neutrinos'])
    # plt.xscale('log')
    # plt.legend()
    # plt.show()
