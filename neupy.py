from nuclide import Nuclide, pprint, pd, np
from databases.load_databases import load_all_fy_databases
from tqdm import tqdm

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
        self.fiss_product_neutrino_data = {fiss_elem: [] for fiss_elem in self.fy}
    
    def weighted_cumulative_neutrinos(self, linear_chain: pd.DataFrame, fission_yield: float):
        linear_chain['weight_cum_neut'] = fission_yield * linear_chain['cum_neut'] * (linear_chain['intensity']/100)
    
    def fission_induced_neutrinos(self, time, **flag_kwargs):
        for fission_element, db in self.fy.items():
            for i, row in tqdm(db.iterrows(), total=db.shape[0]):
                AZI = row.name

                nuclide = Nuclide(AZI=AZI, nubase=self.nuclide_template.nubase, fy=self.fy, config_file=self.nuclide_template.config, nubase_config=self.nuclide_template.nubase_config)
                if not nuclide.found:
                    self.missing_nuclides.append(AZI)
                    continue

                nuclide.concentration_profiles(time, **flag_kwargs)
                for linear_chain in nuclide.decay_chain:
                    self.cumulative_neutrino_profile(linear_chain, **flag_kwargs)
                    self.weighted_cumulative_neutrinos(linear_chain, row.YI)
                
                self.fiss_product_neutrino_data[fission_element].append(nuclide)


        


if __name__ == "__main__":

    neupy = Neupy()
    time = np.logspace(0, 16, 10)
    neupy.fission_induced_neutrinos(time)
    pprint(neupy.fiss_product_neutrino_data)
    # import matplotlib.pyplot as plt
    # neupy = NeutrinoEmission()
    # nuclide = Nuclide(135, 52)
    # nuclide.make_decay_chain()
    # nuclide.break_decay_chain_branches()
    # t = np.logspace(0, 16, 100)
    # nuclide.concentration_profiles(t)
    # for i, linear_chain in enumerate(nuclide.decay_chain):
    #     neupy.cumulative_neutrino_profile(linear_chain, ignore_matter_type=True, ignore_flavour=True)
    #     neupy.check_neutrino_profile_type(linear_chain)
    #     display_names = nuclide.display_decay_chain()[i]
    #     plt.stackplot(t, *linear_chain['cum_neut'].values, labels=display_names.nuclide.values)
    # print(nuclide.display_decay_chain())
    # plt.legend()
    # plt.xscale('log')
    # plt.show()