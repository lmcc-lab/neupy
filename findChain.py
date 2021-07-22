import pandas as pd
import numpy as np
from naming import readable
from dependencies import conv_str_to_list, nMax
import preferences

def find_by_max_yield(nthMax):
    '''
    Find the decay chain data from maxmium yield.

    @param nthMax; int, nth from the maximum.

    returns readableDecayChain, isoDecayData, yeild, fatherDecayConst
    '''
    decayChainData = pd.read_csv('./decay_chain_data/IgnoreTheoreticEnergy='+str(preferences.ignore_energetically_possible)+'_decayChainData.csv')
    allData = pd.read_csv('./decay_chain_data/IgnoreTheoreticEnergy='+str(preferences.ignore_energetically_possible)+'_isomerData.csv')
    allyeild =  allData['Yield'].tolist()
    fatherIndex = allyeild.index(nMax(nthMax, allyeild))
    yeild = allyeild[fatherIndex]
    fatherDecayConst = allData['Decay Constant'].loc[fatherIndex]
    isoDecayData = decayChainData.loc[fatherIndex]
    decayChain = isoDecayData['Decay Chain']
    decayChain = conv_str_to_list(decayChain)
    readableDecayChain = np.zeros(len(decayChain)).tolist()
    for i, gen in enumerate(decayChain):
        if gen != ['']:
            readableGen = [readable(j) for j in gen]
            readableDecayChain[i] = readableGen
    return readableDecayChain, isoDecayData, yeild, fatherDecayConst

def find_by_atomic(A, Z, level=0):
    '''
    Find decay chain data using atomic and proton numbers. Isomer level is also
    an option but defaults to 0.

    @param A; int, atomic number
    @param Z; int, proton number
    @param level=0; int, isomer level

    returns readableDecayChain, isoDecayData, yeild, fatherDecayConst, fatherIndex
    '''
    NU = Z*10000000+A*10000+level
    decayChainData = pd.read_csv('./decay_chain_data/IgnoreTheoreticEnergy='+str(preferences.ignore_energetically_possible)+'_decayChainData.csv')
    allData = pd.read_csv('./decay_chain_data/IgnoreTheoreticEnergy='+str(preferences.ignore_energetically_possible)+'_isomerData.csv')
    allyeild =  allData['Yield'].tolist()
    NUdata = allData['NU isomer'].tolist()
    try:
        fatherIndex = NUdata.index(NU)
    except:
        raise ValueError("Isomer input not in database")
    yeild = allyeild[fatherIndex]
    fatherDecayConst = allData['Decay Constant'].loc[fatherIndex]
    isoDecayData = decayChainData.loc[fatherIndex]
    decayChain = isoDecayData['Decay Chain']
    decayChain = conv_str_to_list(decayChain)
    readableDecayChain = np.zeros(len(decayChain)).tolist()
    for i, gen in enumerate(decayChain):
        if gen != ['']:
            readableGen = [readable(j) for j in gen]
            readableDecayChain[i] = readableGen
    return readableDecayChain, isoDecayData, yeild, fatherDecayConst, fatherIndex

