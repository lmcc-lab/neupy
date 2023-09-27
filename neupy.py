from nuclide import Nuclide, pprint, pd, np, DecayChainDepthWarning, max_chain_depth, _Bateman
from databases.load_databases import load_all_fy_databases
from tqdm import tqdm
from typing import Union, Tuple, Dict, List
from config import *
import pickle


class _CumulativeNeutrinos:

    def __init__(self, neut_scale, concentration_profile: _Bateman) -> None:
        self.neut_scale = neut_scale
        self.concentration_profile = concentration_profile

    def _map_element_wise_multiplication(self, vec1, vec2):
        row = vec1.flatten()
        calc = np.zeros(len(row)).tolist()
        for i, elem in enumerate(row):
            calc[i] = elem * vec2
        calc = np.array(calc)
        calc = calc.reshape(*vec1.shape, len(vec2))
        return calc
    
    def calculate(self, t):
        return self._map_element_wise_multiplication(self.neut_scale, self.concentration_profile(t))

class _WeightedCumNeutrinos:

    def __init__(self, fission_yield: float, total_branch_ratio_intensity: float, cum_neutrinos: _CumulativeNeutrinos):
        self.fission_yield = fission_yield
        self.branch_ratio_intensity = total_branch_ratio_intensity
        self.cum_neutrinos = cum_neutrinos
    
    def calculate(self, t):
        return self.fission_yield * self.branch_ratio_intensity * self.cum_neutrinos(t)

class _TotalCumNeutrinos:

    def __init__(self, nuclide: Nuclide):
        self.nuclide = nuclide
    
    def calculate(self, t):
        total_neutrino_profile = 0
        for linear_chain in self.nuclide.decay_chain:
            total_neutrino_profile += np.sum(np.vstack([cum_function(t) for cum_function in linear_chain['weight_cum_neut'].values]), axis=0)
        return total_neutrino_profile

class _FissionProductTotalNeutrinos:

    def __init__(self, fission_nuclide_array: List[Nuclide]) -> None:
        self.fission_nuclide_array = fission_nuclide_array
    
    def calculate(self, t):
        total_neutrinos = 0
        for nuclide in self.fission_nuclide_array:
            total_neutrinos += nuclide.total_neutrino_profile(t)
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
        if 'concentration_profile' not in linear_chain.columns:
            raise ValueError("concentration_profile must have been calculated before running cumulative_neutrino_profile")
        neut_vec_kwargs.pop('cumulative', None)
        neut_vec = self.neutrino_vec(linear_chain, cumulative=True, **neut_vec_kwargs)
        linear_chain['cum_neut'] = linear_chain.apply(
            lambda row: _CumulativeNeutrinos(neut_vec[row.name], row.concentration_profile).calculate, axis=1)
            
    def check_neutrino_profile_type(self, linear_chain):
        cum_neutrino_test = (linear_chain.shape[0], *linear_chain.loc[0, 'cum_neut'].shape)
        if len(cum_neutrino_test) == 2:
            print("Cumulative neutrinos were calculated without seperating neutrino flavour or matter-antimatter\n"
                  f"Slide index 0: Nuclide in decay chain, with {cum_neutrino_test[0]} nuclides"
                  f"Slide index 1: Time, with length {cum_neutrino_test[1]}")
        elif len(cum_neutrino_test) == 3 and cum_neutrino_test[1] == 2:
            print("Cummulative neutrinos were calculated with seperation of matter-antimatter, but not seperating flavour\n"
                  f"Slice index 0: Nuclide in decay chain, with {cum_neutrino_test[0]} nuclides\n"
                  "Slice index 1: Matter(antimatter) seperation index 0(1)\n"
                  f"Slide index 2: Time, with length {cum_neutrino_test[2]}")
        elif len(cum_neutrino_test) == 3 and cum_neutrino_test[1] == 3:
            print("Cummulative neutrinos were calculated with seperation of flavour, but not matter-antimatter\n"
                  f"Slice index 0: Nuclide in decay chain, with {cum_neutrino_test[0]} nuclides\n"
                  "Slice index 1: flavour, e, mu, tau seperation index with index 0, 1, 2 respectively\n"
                  f"Slide index 2: Time, with length {cum_neutrino_test[2]}")
        elif len(cum_neutrino_test) == 4:
            print("Cummulative neutrinos were calculated with seperation of both flavour and matter-antimatter\n"
                  f"Slice index 0: Nuclide in decay chain, with {cum_neutrino_test[0]} nuclides\n"
                  "Slice index 1: Matter(antimatter) seperation index 0(1)\n"
                  "Slice index 2: flavour, e, mu, tau seperation index with index 0, 1, 2 respectively\n",
                  f"Slide index 3: Time, with length {cum_neutrino_test[3]}")

    
class Neupy(NeutrinoEmission):
    
    def __init__(self, **fy_kwargs) -> None:
        """
        Neupy is a child class of Neutrino emission, which adds extra information about a decay chains
        neutrino profile, specifically the matter and flavour type. 

        Parameters
        -------------------
        - fy_kwargs, passed to databases.load_databases.load_all_fy_databases

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
        self.nuclide_template = Nuclide(1, 1)
        self.missing_nuclides = []
        
    def weighted_cumulative_neutrinos(self, linear_chain: pd.DataFrame, fission_yield: float):
        total_branch_ratio_intensity = np.prod(linear_chain['intensity'])
        linear_chain['weight_cum_neut'] = linear_chain.apply(lambda row: _WeightedCumNeutrinos(fission_yield, total_branch_ratio_intensity, row['cum_neut']).calculate, axis=1)
    
    def total_weighted_neutrinos(self, nuclide: Nuclide):
        nuclide.total_neutrino_profile = _TotalCumNeutrinos(nuclide).calculate
            
    def fission_induced_neutrinos(self, nuclide: Nuclide, fy: float, **flag_kwargs):
        nuclide.concentration_profiles(**flag_kwargs)
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
                                      file_path = 'databases/',
                                      **kwargs) -> dict:
        """
        neutron_energy_range: Tuple[float] | float | None, if None it will include all database data without
        consideration of the neutron energy, if it's a float then it will only include databases where the
        neutron energy is what you gave it, otherwise it expects a tuple representing the start and end 
        of the neutron range in question (start, end) inclusive of both ends.
        
        """
        if load_saved:
            print('loading fission data')
            with open(f'{file_path}{file_name}.pickle', 'rb') as f:
                fission_product_neutrino_data = pickle.load(f)
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
                for i, row in tqdm(db.iterrows(), total=db.shape[0], desc=f"{fission_element} {ne}"):
                    if row.Y == 0.0:
                        continue
                    AZI = i
                    nuclide = Nuclide(AZI=AZI, nubase=self.nuclide_template.nubase, fy=self.fy, config_file=self.nuclide_template.config, nubase_config=self.nuclide_template.nubase_config)
                    # print(AZI)
                    if not nuclide.found:
                        self.missing_nuclides.append(AZI)
                        continue
                    try:
                        self.fission_induced_neutrinos(nuclide, row.Y, **kwargs)
                    except DecayChainDepthWarning:
                        print(f"Nuclide {AZI}'s decay chain depth exceeded maximum {max_chain_depth}")
                        continue
                    nuclide_neutrino_data.append(nuclide)
                element_fpnd[ne] = {"total_neutrinos": _FissionProductTotalNeutrinos(nuclide_neutrino_data).calculate, 'nuclide_specific': nuclide_neutrino_data}
            fission_product_neutrino_data[fission_element] = element_fpnd
        
        if save:
            with open(f'{file_path}{file_name}.pickle', 'wb+') as f:
                pickle.dump(fission_product_neutrino_data, f)

        return fission_product_neutrino_data
    
    def convert_thermal_to_burn_rate(self, thermal_power: np.ndarray) -> dict:
        # Basic thermal model, where U235 makes up 94% of the thermal power generation
        return {'U235': 0.94 * thermal_power / u235_energy_per_fission_MeV}
    
    def cummulative_neutrinos_from_burn_rate(self, burn_rate: dict, time: Union[list, np.ndarray],
                                             assumed_on_prior_s: float=0, 
                                             **fiss_induced_neut_kwargs):
        fy_neutrinos = self.all_fission_induced_neutrinos(time, **fiss_induced_neut_kwargs)
        dt = time[1] - time[0]

        cummulative_neutrinos = dict()
        for element, burn_profile in burn_rate.items():
            element_neutrinos = fy_neutrinos.get(element, None)

            if element_neutrinos == None:
                continue
            
            cum_neut_for_elem = dict()

            for ne, neutrino_data in element_neutrinos.items():
                total_neutrinos = neutrino_data['total_neutrinos'](time)
                cum_neutrinos = 0
                for i, instant_burn_rate in enumerate(burn_profile):
                    neut_array = np.zeros(len(time))
                    neut_array[i:] = total_neutrinos[i:]
                    cum_neutrinos += instant_burn_rate * neut_array * dt
                cum_neut_for_elem[ne] = cum_neutrinos

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

    def neutrino_emission_rate(self, burn_rate: dict,
                               time: Union[list, np.ndarray],
                               assumed_on_prior_s: float = 0,
                               **fiss_induced_neut_kwargs) -> Dict[str, Dict[str, np.ndarray]]:
        
        cum_neutrinos = self.cummulative_neutrinos_from_burn_rate(burn_rate, time, assumed_on_prior_s=assumed_on_prior_s, 
                                                                  **fiss_induced_neut_kwargs)
        emission_rate = dict()
        dt = time[1] - time[0]
        for element, ne_dict in cum_neutrinos.items():
            element_emission_rate = dict()
            for ne, neutrino_data in ne_dict.items():
                element_emission_rate[ne] = np.gradient(neutrino_data, dt)
            emission_rate[element] = element_emission_rate
        
        return emission_rate
            
            


        
        

        


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
