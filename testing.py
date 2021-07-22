import numpy as np
import constants
import breakBranches
import naming
import bateman_1 as bateman
import matplotlib.pyplot as plt
import findChain as fc
from dependencies import conv_str_to_list as ctl
from naming import readable

readableChain, isodata, yeild, fatherDecayConst, _ = fc.find_by_atomic(87, 34)
# readableChain, isodata, yeild, fatherDecayConst = fc.find_by_max_yield(1)
# print(readableChain)
# print(isodata)
# print(yeild)
# print(fatherDecayConst)
father = isodata['NU isomer']
decayChain = ctl(isodata['Decay Chain'])
branchIndex = ctl(isodata['Branch index'])
branchRatio = ctl(isodata['Branch Ratio'])
decayModes = ctl(isodata['Decay Mode'])
decayConstChain = ctl(isodata['Decay Constant'])
brokenBranches = breakBranches.branches(decayChain, branchIndex, branchRatio, decayModes, decayConstChain)
brokenBranches.breakBranches()
    
CDF = 0
PDF = 0
for i, chain in enumerate(brokenBranches.brokenDecayChainBranches):
    readableChain = np.zeros(len(chain)).tolist()
    for j, daught in enumerate(chain):
        readableChain[j] = readable(daught)
    print(readableChain)
    print(brokenBranches.brokenDecayModes[i])
    print(brokenBranches.brokenDecayConsts[i])
    print(brokenBranches.weighting)
    rateEquations = bateman.bateman(father, fatherDecayConst, yeild, brokenBranches.weighting[i], brokenBranches.brokenDecayConsts[i], brokenBranches.brokenDecayModes[i])
    rateEquations.neutrinoVector()
    rateEquations.generator()
    rateEquations.PDF_gen()
    CDF = CDF + rateEquations.weightedCDF
    PDF = PDF + rateEquations.weightedPDF
    plt.plot(constants.t, CDF)
    # plt.plot(constants.t, PDF, linestyle = ':')
    # plt.xscale('log')
    plt.grid()
    plt.show()
    plt.grid()
    plt.plot(constants.t[50000:], PDF[50000:], linestyle = ':')
    plt.show()
    # plt.pause(0.05)
    # plt.close

# readableBranches = np.zeros(len(brokenBranches.brokenDecayChainBranches)).tolist()
# for i, chain in enumerate(brokenBranches.brokenDecayChainBranches):
    # readableBranches[i] = [naming.readable(j) for j in chain]
    # print(readableBranches[i], brokenBranches.brokenBranchRatios[i], brokenBranches.weighting[i])

