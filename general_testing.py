import numpy as np
from numpy import exp
import matplotlib.pyplot as plt
import pandas as pd
from sympy import *
import constants
from scipy.integrate import odeint
import dependencies as dp

data = pd.read_csv('./Contributing_chains/NotSimpleConstributingChainsMaxCDF.csv')

def check_for_Xe(father, Chain):
    if father == '135Xe-54':
        flag = True
    else:
        for daughter in Chain:
            if daughter == '135Xe-54':
                flag = True
                break
            else:
                flag = False
    return flag

class neutronAbsorbtion():

    def __init__(self, n_I, t, x0):
        self.time = t
        self.n_I = n_I
        self.x0 = x0

    def solve(self):
        def Xe_chain_numerical(x0, t):
            n_X = x0[0]
            n_Cs = x0[1]
            n_Ba = x0[2]
            n_Xe135 = x0[3]

            dn_Xdt = eval('- constants.XeDC * n_X + constants.IDC * ('+self.n_I+') - constants.sigma_x * n_X * constants.phi_0')
            dn_X136dt = constants.sigma_x * n_X * constants.phi_0
            dn_Csdt = - constants.CsDC * n_Cs + constants.XeDC * n_X
            dn_Badt = constants.CsDC * n_Cs

            return [dn_Xdt, dn_Csdt, dn_Badt, dn_X136dt]

        x = odeint(Xe_chain_numerical, self.x0, self.time)
        self.nx = x[:, 0]
        self.ncs = x[:, 1]
        self.nba = x[:, 2]
        self.nXe = x[:, 3]

t = np.logspace(-2, 30, 200)
counter = 0
IDCinBateman = '2.926153244511758e-05'
diff = 0
full_CDF = 0
max_contrib = 0
with open('./Contributing_chains/NotSimpleCDF_full.txt', 'r') as f:
    full = f.readline()
full_CDF_orig = np.array(list(map(eval('lambda t:' + full), t)))

for index, row in data.head(n=data.shape[0]).iterrows():
    Chain = dp.conv_str_to_list(row['Chain'])
    father = row['Father']
    bateman = row['Weighted CDF solution']
    max_contrib = max_contrib + row['Max CDF contrib']
    flag = check_for_Xe(father, Chain)
    print(index)
    if flag == False:
        CDF = np.array(list(map(eval('lambda t:' + bateman), t)))
        full_CDF = full_CDF + CDF
    else:
        if father == '135Xe-54':
            # Use the numerical solution for Xe, Ba, Cs chain
            numerical = neutronAbsorbtion('0', t, [1, 0, 0, 0])
            numerical.solve()
            weighting = float(bateman.split('*')[0])
            nx = 0*numerical.nx
            ncs = 1*numerical.ncs
            nba = 2*numerical.nba
            nXe136 = 0*numerical.nXe
            counter +=1
            CDF = weighting*(nx + ncs + nba + nXe136)

            orig_CDF = np.array(list(map(eval('lambda t:' + bateman), t)))

            full_CDF = full_CDF + CDF

            diff = diff + abs(orig_CDF[-1] - CDF[-1])
            print(father, Chain)
            # plt.plot(t, CDF, color='k')
            # plt.plot(t, orig_CDF, color='b')
            # plt.xscale('log')
            # plt.show()

        else:
            for i, daughter in enumerate(Chain):
                if daughter == '135Xe-54' and bateman != '0':
                    # Get the idodine solution
                    print(father, Chain)
                    # print(bateman, end='\n\n')
                    Ival = bateman.split('2.92615324451175',1)[1].split('e',1)[0]
                    weighting = float(bateman.split('*')[0])
                    Ival = '2.92615324451175'+Ival+'e-05'
                    chunks = bateman.split('exp(-'+Ival+'*t)/(', 1)
                    denom = chunks[1].split('))',1)[0]
                    bateman_to_I = chunks[0] + 'exp(-'+Ival+'*t)/(' + denom +'))'
                    extrap = bateman_to_I.split('*(')
                    Isolution = extrap[-1]
                    startingProd = extrap[-2].split('*')[-1]
                    neutrinoNumber = extrap[-2].split('*')[-2].split('+')
                    if len(neutrinoNumber) > 1:
                        neutrinoNumber = neutrinoNumber[1]
                    else:
                        neutrinoNumber = neutrinoNumber[0]

                    Isolution = startingProd + '*(' + Isolution
                    first_bateman = bateman.split(neutrinoNumber+'*'+Isolution, 1)[0]
                    first_bateman = first_bateman.split('*',1)[1]
                    first_bateman = first_bateman[:-1]

                    neutrinoVector = [float(neutrinoNumber) + i for i in range(4)]
                    # print(first_bateman)


                    numerical = neutronAbsorbtion(Isolution, t, [0, 0, 0, 0])
                    numerical.solve()
                    nx = neutrinoVector[1]*numerical.nx
                    ncs = neutrinoVector[2]*numerical.ncs
                    nba = neutrinoVector[3]*numerical.nba
                    nXe = neutrinoVector[1]*numerical.nXe
                    Isolution = neutrinoVector[0]*np.array(list(map(eval('lambda t:'+Isolution), t)))
                    if first_bateman != '':
                        first_bateman = list(map(eval('lambda t:'+first_bateman +')'), t))
                    else:
                        first_bateman = list(map(eval('lambda t: 0'), t))

                    CDF = first_bateman + Isolution + nx + ncs + nba + nXe
                    CDF = weighting * CDF
                    orig_CDF = np.array(list(map(eval('lambda t:' + bateman), t)))
                    full_CDF = full_CDF + CDF

                    diff = diff + abs(orig_CDF[-1] - CDF[-1])
                    counter += 1
                    # plt.plot(t, Isolution, label = 'I')
                    # plt.plot(t, nx, label='Xe')
                    # plt.plot(t, nXe, label='Xe-136')
                    # plt.plot(t, ncs, label='Cs')
                    # plt.plot(t, nba, label='Ba')
                    # plt.plot(t, CDF, color = 'k')
                    # plt.plot(t, orig_CDF, color ='b')
                    # plt.xscale('log')
                    # plt.show()
plt.plot(t, full_CDF_orig, color='k', label = 'Original CDF')
plt.plot(t, full_CDF, color='b', linestyle ='--', label = 'Xe poisoning')
plt.legend()
plt.grid()
plt.title('Xe modulation of cumulative neutrinos')
plt.xlabel('Time (s)')
plt.ylabel('Cumulative neutrinos')
plt.xscale('log')
plt.show()

print(diff, counter)
percent_change = diff/max_contrib
print(percent_change)