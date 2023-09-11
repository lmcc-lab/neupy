from nuclide import Nuclide, pprint, pd, np
from databases.load_databases import load_all_fy_databases


class NeutrinoEmission:

    def __init__(self, **fy_kwargs) -> None:
        self.fy = load_all_fy_databases(**fy_kwargs)
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
            print("Cumulative neutrinos were calculated without seperating neutrino flavour or matter-antimatter")
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
                  f"Slide index 2: Time, with length {cum_neutrino_test[2]}")

    
class Neupy:

    def __init__(self) -> None:
        pass


if __name__ == "__main__":
    import matplotlib.pyplot as plt
    neupy = NeutrinoEmission()
    nuclide = Nuclide(135, 52)
    nuclide.make_decay_chain()
    nuclide.break_decay_chain_branches()
    t = np.logspace(0, 16, 100)
    nuclide.concentration_profiles(t)
    for linear_chain in nuclide.decay_chain:
        neupy.cumulative_neutrino_profile(linear_chain, ignore_matter_type=True, ignore_flavour=True)
        neupy.check_neutrino_profile_type(linear_chain)
        plt.stackplot(t, *linear_chain['cum_neut'].values)
    plt.xscale('log')
    plt.show()