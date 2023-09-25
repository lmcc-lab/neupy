from nuclide import Nuclide, pprint, pd, np, DecayChainDepthWarning, max_chain_depth
from databases.load_databases import load_all_fy_databases
from tqdm import tqdm
from typing import Union, Tuple
from config import *
import pickle

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
    
    def _map_element_wise_multiplication(self, vec1, vec2):
        row = vec1.flatten()
        calc = np.zeros(len(row)).tolist()
        for i, elem in enumerate(row):
            calc[i] = elem * vec2
        calc = np.array(calc)
        calc = calc.reshape(*vec1.shape, len(vec2))
        return calc
    
    def cumulative_neutrino_profile(self, linear_chain: pd.DataFrame, **neut_vec_kwargs):
        if 'concentration_profile' not in linear_chain.columns:
            raise ValueError("concentration_profile must have been calculated before running cumulative_neutrino_profile")
        neut_vec_kwargs.pop('cumulative', None)
        neut_vec = self.neutrino_vec(linear_chain, cumulative=True, **neut_vec_kwargs)
        linear_chain['cum_neut'] = linear_chain.apply(
            lambda row: self._map_element_wise_multiplication(neut_vec[row.name], row.concentration_profile), axis=1)
            
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
        super().__init__()
        self.fy = load_all_fy_databases(**fy_kwargs)
        self.nuclide_template = Nuclide(1, 1)
        self.missing_nuclides = []
        
    def weighted_cumulative_neutrinos(self, linear_chain: pd.DataFrame, fission_yield: float):
        linear_chain['weight_cum_neut'] = fission_yield * linear_chain['cum_neut'] * np.prod(linear_chain['intensity'])
    
    def total_weighted_neutrinos(self, nuclide: Nuclide):
        for linear_chain in nuclide.decay_chain:
            nuclide.total_neutrino_profile += np.sum(np.vstack(linear_chain['weight_cum_neut'].values), axis=0)
            
    def fission_induced_neutrinos(self, nuclide, fy, time, **flag_kwargs):
        nuclide.concentration_profiles(time, **flag_kwargs)
        for linear_chain in nuclide.decay_chain:
            self.cumulative_neutrino_profile(linear_chain, **flag_kwargs)
            self.weighted_cumulative_neutrinos(linear_chain, fy)
        self.total_weighted_neutrinos(nuclide)
    
    def all_fission_induced_neutrinos(self, time: Union[list, np.ndarray], 
                                      neutron_energy_range: Union[Tuple[float], float, None] = None, 
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
        
        fission_product_neutrino_data = dict()
        for fission_element, ne_dict in self.fy.items():
            element_fpnd = dict()
            for ne, db in ne_dict.items():
                
                if ne != 'independent':
                    if isinstance(neutron_energy_range, float) and ne != neutron_energy_range:
                        continue
                    
                    if isinstance(neutron_energy_range, tuple) and not neutron_energy_range[0] <= ne <= neutron_energy_range[1]:
                        continue
                
                total_neutrinos = 0
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
                        self.fission_induced_neutrinos(nuclide, row.Y, time, **kwargs)
                    except DecayChainDepthWarning:
                        print(f"Nuclide {AZI}'s decay chain depth exceeded maximum {max_chain_depth}")
                        continue
                    nuclide_neutrino_data.append(nuclide)
                    total_neutrinos += nuclide.total_neutrino_profile
                element_fpnd[ne] = {"total_neutrinos": total_neutrinos, 'nuclide_specific': nuclide_neutrino_data}
            fission_product_neutrino_data[fission_element] = element_fpnd
        
        if save:
            with open(f'{file_path}{file_name}.pickle', 'wb+') as f:
                pickle.dump(fission_product_neutrino_data, f)

        return fission_product_neutrino_data
    
    def convert_thermal_to_burn_rate(self, thermal_power: np.ndarray) -> dict:
        # Basic thermal model, where U235 makes up 94% of the thermal power generation
        return {'U235': 0.94 * thermal_power / u235_energy_per_fission_MeV}
    
    def cummulative_neutrinos_from_burn_rate(self, burn_rate: dict, time: Union[list, np.ndarray], 
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
                total_neutrinos = neutrino_data['total_neutrinos']
                cum_neutrinos = 0
                for i, instant_burn_rate in enumerate(burn_profile):
                    neut_array = np.zeros(len(time))
                    neut_array[i:] = total_neutrinos[i:]
                    cum_neutrinos += instant_burn_rate * neut_array * dt
                cum_neut_for_elem[ne] = cum_neutrinos

            cummulative_neutrinos[element] = cum_neut_for_elem

        return cummulative_neutrinos

            
            


        
        

        


if __name__ == "__main__":

    neupy = Neupy()
    time = np.linspace(0, 1000, 1000)
    # pprint(neupy.fy)
    thermal_power = np.zeros(len(time))
    thermal_power[100:200] = 20e6  # 20 MW

    cum_neutrinos = neupy.cummulative_neutrinos_from_burn_rate({'U235': thermal_power}, time, load_saved=True)
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