import numpy as np
import matplotlib.pyplot as plt
from numpy.lib.utils import deprecate
import pandas as pd
import dependencies as dp
import constants
import bateman as bt


def thermal(time, tOn, tOff, power):
    thermal_output = np.zeros(len(time))
    thermal_output[dp.find_index(time, tOn) : dp.find_index(time, tOff)] = power
    return thermal_output

def fission(thermal_power):
    return thermal_power/constants.U235_fiss

def linear():
    graph_off = 10**5
    time = np.linspace(0,int(graph_off), int(graph_off/2+1))
    tOn = 0
    tOff = 10**4

    #Toy model 
    #One fission per second
    thermal_output = thermal(time, tOn, tOff, constants.U235_fiss)
    fissionRate = fission(thermal_output)
    #Yields
    Ayeild = 1
    #Decay Chain
    #Seconds
    AFatherDC = np.log(2)/(1)
    AdecayChain = ['A1']
    AdecayModes = ['b-']
    AbranchRatio = [100]
    AbranchIndex = [0]

    ANeutrinoCDF = bt.bateman(270660000, AFatherDC, Ayeild, 1, [0], AdecayModes, time).neutrinoVector().CDF_gen().PDF_gen() #show_stacked=True, decayChain=[270660000]

    #Decay Chain
    #Seconds
    Byeild = 1
    BFatherDC = np.log(2)/(100)
    BdecayChain = ['B1', 'B2']
    BdecayModes = ['b-', 'b-']
    BbranchRatio = [100, 100]
    BbranchIndex = [0, 00]
    BdecayConsts = [np.log(2)/(10000), 0]

    BNeutrinoCDF = bt.bateman(270660000, BFatherDC, Byeild, 1, BdecayConsts, BdecayModes, time).neutrinoVector().CDF_gen() #show_stacked=True, decayChain=[270660000, 270660000]

    #This is the total neutrino CDF for a single thermal event
    totalCDF = ANeutrinoCDF.weightedCDF + BNeutrinoCDF.weightedCDF


    #Total neutrino CDF for thermal reactor
    CummulitiveneutrinoFlux = 0
    for i in range(len(time)):
        if i>0:
            newNeutrinoCDF = np.roll(totalCDF, i)
            newNeutrinoCDF[:i] = 0
            CummulitiveneutrinoFlux = CummulitiveneutrinoFlux + (time[i]-time[i-1])*fissionRate[i]*newNeutrinoCDF
        else:
            CummulitiveneutrinoFlux = CummulitiveneutrinoFlux + fissionRate[i]*totalCDF

    neutrinoFlux = dp.derivative(CummulitiveneutrinoFlux, time)

    offset_time = 0

    startPlot = dp.find_index(time, offset_time)
    newTime = time-offset_time
    plt.plot(newTime[startPlot :], neutrinoFlux[startPlot :])
    plt.vlines(tOff-offset_time, 0, max(neutrinoFlux[startPlot :]))
    plt.xscale('log')
    plt.ylabel('Neutrino Flux ('+r'$\bar{\nu}s^{-1}$)')
    plt.xlabel('Time (s)')
    plt.show()

def log():
    graph_off = 5
    graph_on = -3
    time = np.logspace(graph_on,graph_off, int(100))
    tOn = 10**-2
    tOff = 10**4

    #Toy model 
    #One fission per second
    thermal_output = thermal(time, tOn, tOff, constants.U235_fiss)
    fissionRate = fission(thermal_output)
    #Yields
    Ayeild = 1
    #Decay Chain
    #Seconds
    AFatherDC = np.log(2)/(1)
    AdecayChain = ['A1']
    AdecayModes = ['b-']
    AbranchRatio = [100]
    AbranchIndex = [0]

    ANeutrinoCDF = bt.bateman(270660000, AFatherDC, Ayeild, 1, [0], AdecayModes, time).neutrinoVector().CDF_gen().PDF_gen() #show_stacked=True, decayChain=[270660000]

    #Decay Chain
    #Seconds
    Byeild = 1
    BFatherDC = np.log(2)/(100)
    BdecayChain = ['B1', 'B2']
    BdecayModes = ['b-', 'b-']
    BbranchRatio = [100, 100]
    BbranchIndex = [0, 00]
    BdecayConsts = [np.log(2)/(10000), 0]

    BNeutrinoCDF = bt.bateman(270660000, BFatherDC, Byeild, 1, BdecayConsts, BdecayModes, time).neutrinoVector().CDF_gen() #show_stacked=True, decayChain=[270660000, 270660000]

    #This is the total neutrino CDF for a single thermal event
    totalCDF = ANeutrinoCDF.weightedCDF + BNeutrinoCDF.weightedCDF


    #Total neutrino CDF for thermal reactor
    CummulitiveneutrinoFlux = 0
    for i in range(len(time)):
        if i > 0:
            newNeutrinoCDF = np.roll(totalCDF, i)
            newNeutrinoCDF[:i] = 0
            plt.plot(time, newNeutrinoCDF)
            plt.xscale('log')
            plt.show()
            CummulitiveneutrinoFlux = CummulitiveneutrinoFlux + (time[i]-time[i-1])*fissionRate[i]*newNeutrinoCDF

    neutrinoFlux = dp.derivative(CummulitiveneutrinoFlux, time)

    plt.plot(time, totalCDF)
    plt.vlines(tOff, 0, max(neutrinoFlux))
    plt.xscale('log')
    plt.ylabel('Neutrino Flux ('+r'$\bar{\nu}s^{-1}$)')
    plt.xlabel('Time (s)')
    plt.show()
linear()