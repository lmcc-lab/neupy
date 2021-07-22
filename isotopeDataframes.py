'''
Author: Liam McClelland

Last Edited: 31/05/21
'''

import decayChain
import nubaseSorter as nuSorter
import pandas as pd
import matplotlib.pyplot as plt
import preferences
import os
'''
isotopeDataframes generates the first generation of daughters from the fission products, 
along with their associated branch ratios, decay modes, decay constants etc. Data is 
saved in IgnoreTheoreticEnergy='True or False'_isomerData.csv. True of False based on 
preferences.ignore_energetically_possible. This gives all of the fission products their 
first generation of daughters.

'''
path = '.'
outdir=path+'/decay_chain_data/'
if (os.access(outdir,os.F_OK))==0:
    os.mkdir(outdir)

allData = {'NU isomer': [],'Yield': [],'Yield uncertainty': [],'Energy level': [],'Decay Constant':[],'Daughters':[],'Branch Ratios': [],'Decay modes':[]}
allData = pd.DataFrame(allData, columns = ['NU isomer','Yield','Yield uncertainty','Energy level','Decay Constant','Daughters','Branch Ratios','Decay modes'])

allDecayChains = {'NU isomer': [], 'Decay Chain':[], 'Branch index':[], 'Branch Ratio':[], 'Decay Constant': [], 'Decay Mode':[], 'Missing isomers': []}
allDecayChains = pd.DataFrame(allDecayChains, columns = ['NU isomer', 'Decay Chain', 'Branch index', 'Branch Ratio', 'Decay Constant', 'Decay Mode', 'Missing isomers'])


for i,a in enumerate(nuSorter.A):
    z = nuSorter.Z[i]
    level = nuSorter.Level[i]
    yeild = nuSorter.Yeild[i]
    yeildUncert = nuSorter.Yeild_uncert[i]
    fissProd = nuSorter.isomer(z, a, level, yeild, yeildUncert)
    try:
        fissProd.find_nubase()
        fissProd.excit_states()
        fissProd.half_life()
        fissProd.decay_info()
        fissProd.decay_constant()
        fissProd.daughter()
        chain = decayChain.chain(fissProd.daughters, fissProd.branchRatio, fissProd.decay_const, fissProd.decayModes)
        chain.chain_gen()
        addData = pd.DataFrame([[int(fissProd.NU), fissProd.yeild, fissProd.yeild_uncert, fissProd.level, fissProd.decay_const, fissProd.daughters, fissProd.branchRatio, fissProd.decayModes]],columns = ['NU isomer','Yield','Yield uncertainty','Energy level','Decay Constant','Daughters','Branch Ratios','Decay modes'])
        allData = allData.append(addData, ignore_index = True)
        addChainData = pd.DataFrame([[int(fissProd.NU), chain.isomerChain, chain.chainIndex, chain.branchRatioChain, chain.decayConstChain, chain.decayModeChain, chain.Missing]], columns = ['NU isomer', 'Decay Chain', 'Branch index', 'Branch Ratio', 'Decay Constant', 'Decay Mode', 'Missing isomers'])
        allDecayChains = allDecayChains.append(addChainData, ignore_index = True)
    except ValueError:
        fissProd.missing_isomer()


allData.to_csv(outdir+'IgnoreTheoreticEnergy='+str(preferences.ignore_energetically_possible)+'_isomerData.csv',sep = ',')
allDecayChains.to_csv(outdir+'IgnoreTheoreticEnergy='+str(preferences.ignore_energetically_possible)+'_decayChainData.csv',sep = ',')