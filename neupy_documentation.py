'''
Author: Liam McClelland
Last Edited: 31/05/2021

Neupy is a collection of scripts for use in sorting nuclear data and generating new data. 
It's main goal is to generate neutrino flux and spectrum for use in toy reactor models and for 
experimental data. There is however other things neupy can do.

Neupy can be seperated into 2 parts, the parsing of nuclear data and generation of sorted data
and then the neutrino emission part. These two parts should be run seperately to save computing
time. 

This script also has a few dependencies including pandas, numpy, matplotlib, math, os, scipy.

                                                nubaseSorter.py
__________________________________________________________________________________________________
Starting with it's parsing of nuclear data, you can access the raw nubase and fission yeild data
using

fission yeild data = nubaseSorter.fiss_data
nuclear data = nubaseSorter.nubase

The original databases should be kept in a folder ./databases/ with file headings fyu235thermal.txt 
and Parsed_nubase2016.xlsx for the fission data and nubase data respectively. 
'''

import nubaseSorter as ns

# fissionData = ns.fiss_data            #Get all raw fission data in dataframe object.
# nubaseData = ns.nubase                #Get all raw nubase data in dataframe object.

# fiss_column_Z = ns.Z                  #Get column Z from fission data, list object.
# print(len(fiss_column_Z))
# ns.A, ns.Level, ns.Yeild, ns.Yeild_uncert all accessable

# nubase_column_Z = ns.nuZ              #Get column A from nubase data, list object.
# ns.nuA is also accessable.

'''
Once nubaseSorter is imported, the databases are loaded into these objects. 
You can then call the isomer() class for sorting through nubase and gathering relavent data.

isomer class has 7 methods, including
        find_nubase()
        excit_states()
        half_life()
        missing_isomer()
        decay_info()
        decay_constant()
        daughter()

These objects should be called in this order.
One can only initialise the isomer class with a single isotope.
'''
# Example isotope
# Z = 4
# A = 7
# level = 1
# yeild = 1
# yeild_uncert = 0
# sortingData = ns.isomer(Z, A, level, yeild, yeild_uncert)
# sortingData.find_nubase()
# print(sortingData.nubase_data)
# sortingData.excit_states()
# print(sortingData.isomer_states)

# help(ns.isomer)        # Read the documentation using.


'''
                                        IsotopeDataframe.py
_________________________________________________________________________________________________
The isotope dataframes are then generated using the isotopeDataframes module. It simply runs
through all of the fission yield isotopes and finds all the accociated information from
nubaseSorter and saves the data into 2 .csv files named 
IgnoreTheoreticEnergy='True or False'_isomerData.csv
IgnoreTheoreticEnergy='True or False'_decayChainData.csv

These are saved in the folder ./decay_chain_data/. True or False are set by changing 
the ignore_energetically_possible toggle in preferences.py. This is used in the chain generation
part of this loop. If it is True then it will ignore branches that have nan values, if False
it will calculate what the branch ratio would be based off of other branches. If the sum
of the branch ratios is > 100 then it sets this branch to 0.

This module should be run to generate the parsed databases, but once run it does not need to be
run again.
'''


'''
                                                decayChain.py
________________________________________________________________________________________________
The decayChain module deals with generating the next generation of daughters, including grabbing
all relevant information for the current generation. 

The decayChain module has a chain class, which is initialised with a list of daughters, branch
ratios, decay constants and decay modes. This list should be the last generation of the current
decay chain.

Once initialised, the chain_gen function can be called, which then generates the next generation
and appends it to the list of generation information. There is also a preference here to set 
a maximum decay chain depth. This is done through preferences.max_decay_chain_depth. Once exceeded
the next generation will not be generated, a message will be printed explaining that this cap 
has been reached and a chainDepthExceeded counter will count up. You can access the number of
times the chain depth exceeds the maximum preference by printing self.chainDepthExceeded.

See the example below
'''

import decayChain as dc
import naming

#View all documentation of this module here.
# help(dc.chain)

# Z = 24
# A = 66
# level = 0
# isotopeNU = naming.NU(Z, A, level)

# print(isotopeNU)

#Get all data from you isotope
# isotopeData = ns.isomer(Z, A, level, 1, 0).find_nubase().excit_states().half_life().missing_isomer().decay_info().decay_constant().daughter()

#Call any of the data from the istopeData object
# branchRatio = isotopeData.branchRatio
# print(branchRatios)

#Get chain data object
# chainData = dc.chain(isotopeData.daughters, isotopeData.branchRatio, isotopeData.decay_const, isotopeData.ground_decay_modes).chain_gen()

#For example, see the full decay chain using
# print(chainData.isomerChain)

# Make it easily readable using
# readableChain = naming.readableChain(chainData.isomerChain)
# for _, generation in enumerate(readableChain):
#    print(generation)

'''
                                                naming.py
______________________________________________________________________________________________
This is a good time to introduce one of the side modules, naming.py. This module handles
changing naming types of the isotopes. The different types include, NU format (int), seperated 
form (list) and readable format (string). There are a few functions for these different
conversions, which include

readable
seperate
short_readable
NU
readableChain

The functions are pretty self explanitary, but for a quick description, readable takes a NU 
formatted isotope name and makes it easily readable in the usual Z-ID-A format. ID is the 
isotope symbol, i.e. H for hydrogen. Short readable does the same as readable except it 
returns Z-ID. Seperate function takes a NU format and seperates it into a list with 
[Z, A, level]. NU does the opposite, taking a Z, A and level and returning the NU format.
readableChain takes a full decayChain all in NU format and converts them to readables.

'''

'''
                                                breakBranches.py
______________________________________________________________________________________________
Up to this point, most modules (except naming.py) are exclusively used in the isotopeDataframes
module for the first half of the code and all data is saved in the decay_chain_data folder. 
This next part then grabs this data and seperates the branches into seperate linear branches.

This module has a branches class, which is initialised using a decayChain, chainIndex and so
on (run help(breakBranches.branches)) then has a breakBranches function. This function
seperates each branch and returns new lists with each chain and there associated weighting
and other data. See the example below.
'''

import breakBranches as bb

#setup using isotope data
# Z = 24
# A = 66
# level = 0
# isotopeNU = naming.NU(Z, A, level)
# daughters = [[isotopeNU]]
# isotopeData = ns.isomer(Z, A, level, 1, 0).find_nubase().excit_states().half_life().missing_isomer().decay_info().decay_constant().daughter()
# chainData = dc.chain(isotopeData.daughters, isotopeData.branchRatio, isotopeData.decay_const, isotopeData.ground_decay_modes).chain_gen()

#Break branches, data contained in brokenBranches
# brokenBranches = bb.branches(chainData.isomerChain, chainData.chainIndex, chainData.branchRatioChain, chainData.decayModeChain , chainData.decayConstChain).breakBranches()

# Make these broken branches readable
# brokenChains = naming.readableChain(brokenBranches.brokenDecayChainBranches)

#print each branch with its weighting on its own line.
# for i,gen in enumerate(brokenChains):
#    print(gen, brokenBranches.weighting[i])

'''
                                                bateman.py
______________________________________________________________________________________________
The bateman module deals with generating neutrino vectors, CDFs and PDFs for neutrino emissions.
It works off of the bateman equation of the form
                        n           n-1       n         n
                Nn(t) = Σ [Ni(0) x ( ∏ λj) x (Σ e^-λjt/(∏ (λp-λj)))] 
                        i=1          j=i      j=i     p=i, p≠j
Where Nn(t) is the concentration of the next isotope based off of the current and next isotopes
decay constants. Assuming the first isotopes concentration is 1 (1 element) allows us to scale
this by the neutrino vector to get neutrino yeild instead of isotope concentrations. 

For example, a chain with decay modes
                        [b-,  b-,  2b-]
Would have a neutrino vector of
                        [1,   2,   4]
Since a b- decay mode releases 1 electron antineutrino, and 2b- releases 2 electron antineutrinos.
By calculating the concentration of the nth isotope in the broken chain, and then multiplying
it by the nth neutrino vector component, you are left with the contribution to the neutrino CDF
by the nth isotope. This is equivilent to a dot product of the isotope CDF vector and the
neutrino vector
                        solution = [CDF1, CDF2, CDF3] dot [1, 2, 4]
This solution is then multiplied by a weighting which comes from the weighting of the broken chain
and the yeild of the father isotope. Summing all of these solutions for all yeild products and
all chains gives a total CDF for the neutrino yeild. Literature has this approach ~6, which is
true also for this script. 

This module has 3 methods

neutrinoVector
CDF_gen
PDF_gen

The PDF is the derivative of the CDF and is calculated using the derivative of the bateman
equation. 
                             n           n-1           n         n
                d/dt Nn(t) = Σ [Ni(0) x ( ∏ λj) x -λj*(Σ e^-λjt/(∏ (λp-λj)))] 
                             i=1          j=i         j=i      p=i, p≠j

It gives the probability distribution function for neutrino emission for the nth isotope 
in a chain. 

See the example below. Note that the example CDF goes past 6 neutrinos. This is because the
yeild is set to 1, which isn't true.

CDF_gen has the optional parameter to have show_stacked, which displays a stacked chart of
all of the concentrations of all of the isotopes in a stacked plot. Notice that the stack
remains at 1, therefore preserving concentration.
'''

import bateman
import matplotlib.pyplot as plt
import numpy as np

#setup using isotope data
Z = 53
A = 135
level = 0
isotopeNU = naming.NU(Z, A, level)
isotopeData = ns.isomer(Z, A, level, 1, 0).find_nubase().excit_states().half_life().missing_isomer().decay_info().decay_constant().daughter()
chainData = dc.chain(isotopeData.daughters, isotopeData.branchRatio, isotopeData.ground_decay_modes).chain_gen()
print(chainData.isomerChain, chainData.chainIndex, chainData.branchRatioChain, chainData.decayModeChain, chainData.decayConstChain)
brokenBranches = bb.branches(chainData.isomerChain, chainData.chainIndex, chainData.branchRatioChain, chainData.decayModeChain , chainData.decayConstChain).breakBranches()
brokenChains = naming.readableChain(brokenBranches.brokenDecayChainBranches)

#Setup a time axis
time = np.logspace(0,15,100)
#Use bateman module to generate CDF
# fullCDF = 0
for i, dec in enumerate(brokenBranches.brokenDecayConsts):
    print(brokenChains[i])
    batemanData = bateman.bateman(isotopeNU, isotopeData.decay_const, 1, brokenBranches.weighting[i], dec, brokenBranches.brokenDecayModes[i]).neutrinoVector().CDF_gen(xlog = True, show_stacked=True, decayChain = brokenBranches.brokenDecayChainBranches[i], time=time)
    #Plot each chains CDF
    # plt.plot(time, batemanData.weightedCDF)
    # plt.xlabel('time (s)')
    # plt.ylabel('Cumulitive neutrinos')
    # plt.xscale('log')
    # plt.show()
    # fullCDF = fullCDF + batemanData.weightedCDF

#Plot the full CDF
# plt.plot(time,fullCDF)
# plt.xscale('log')
# plt.xlabel('time (s)')
# plt.ylabel('Cumulitive neutrinos')
# plt.show()

'''
                                    save_updated_data.py
______________________________________________________________________________________________
This module is designed to help speed up analysis by saving all known data into .csv files. 
It begins with the export_data class, which is initialised with a time parameter. This class
automatically grabs the saved data from the isotopeDataframes module which is kept in the
decay_chain_data folder. It also creates a new folder called Contrituting_chains for storing
it's own CDF and PDF data. 

This class has can be run faster by chaning the preferences.simple toggle. If True, then the
script will only consider full chains that are linear, i.e. Xe135 - I135 - C135 and only b-
decay modes. So a chain like [b-, b-, b-, b-] would be considered simple. If false, then 
this script will consider all chains, no matter the complexcity.

This class has 2 methods, contributers and maxContribtions. The first deals with generating
all contributing CDF and PDF data for each chain individually. The second then goes through
these contributions and adds them to a total CDF and PDF based on the top contribtions,
where the script will stop if the contributions meet a certain number of chains or a
certain percent of the full. See help(save_updated_data) for more notes.
'''

import save_updated_data as sud
import preferences
import pandas as pd
from scipy.integrate import quad
import constants as ct
from dependencies import find_index as fi

#Setup a time axis
starting_time = -2
ending_time = 5
time = np.logspace(starting_time, ending_time, 100) #int(ending_time-starting_time)+1
thermal = np.zeros(len(time))
tOn = 10**1
tOff = 10**4

thermal[fi(time, tOn):fi(time,tOff)] = ct.Reactor_output
#
# #Export the CDF and PDF contributions by maximum contribution.
# CDF_output = sud.export_data().contributers().maxContributions(full=True)
# with open('./Contributing_chains/'+preferences.simpleTitle+'CDF_full.txt', 'r') as f:
#     CDF = f.readline()
# #
# #
# # print(CDF)
# CDF_lam = eval('lambda t: '+CDF)
# print(CDF_lam)
# CDF = list(map(CDF_lam, time))
# #
# plt.plot(time, CDF)
# plt.xscale('log')
# plt.show()
# print(CDF[-1])
'''
                                    ContributingChartGenerator.py
______________________________________________________________________________________________
The contributingChartGenerator module is dedicated solely to producing a svg chat with all of
the fission yield products in the standard A vs Z format. They are then colour coded to show
much of a contribution each one makes to the CDF. It is initialised using the 
running_contributions and running_percent from maxContributions method in the save_updated_data module. This svg
file is saved in the Contributing_chains folder and is titled 
contributing_Isomers_*running percent*.svg
'''

import ContributingChartGenerator as ccg

# time = np.logspace(-2,10,100)
# _, _, _, running_contribtions, running_percent = sud.export_data(time).maxContributions()
# ccg.gen_svg(running_contribtions, running_percent)

'''
                                    findChain.py
______________________________________________________________________________________________
More of an independent module, obviously to find a decay chain. It has 2 functions,
find_by_max_yield and find_by_atomic. These two are fairly obvious in there functions. See
help(findChain.function) for notes.
'''
import findChain as fc

# nthMax = 4
# readableDecayChain, isoDecayData, yeild, fatherDecayConst = fc.find_by_max_yield(nthMax)

#Find the Iodine 135 decay chain data, useful for Xenon poisoning analysis.
# readableDecayChain, isoDecayData, yeild, fatherDecayConst, fatherIndex = fc.find_by_atomic(135,53, 0)
# print(readableDecayChain)
# print(isoDecayData)

'''
                                    neutrinoFlux.py
______________________________________________________________________________________________
This module is designed to output the neutrinoFlux given with a U235 burn rate function. This
can be calculated from the thermal output of the reactor.
'''
import neutrinoFlux as nf
import constants as ct
from dependencies import conv_str_to_list as ctl
from dependencies import find_index as fi

# time = np.logspace(0,6,100)
# # sud.export_data(time).contributers().maxContributions(full=True)
# thermal = np.zeros(len(time))
# tOn = 10**1
# tOff = 10**4
# #
# thermal[fi(time, tOn):fi(time,tOff)] = ct.Reactor_output
# #
# CDF_data = pd.read_csv('./Contributing_chains/'+preferences.simpleTitle+'CDF_full.csv')
# CDF = CDF_data.iloc[0,1]
# flux = nf.neuFlux(time, thermal, CDF).flux().neutrinoFlux
# plt.plot(time, thermal)
# plt.plot(time, flux)
# time_offset = 10**4
# start_plot = fi(time, time_offset)
# newTime = time-time_offset
# plt.plot(newTime[start_plot:], flux[start_plot:])
# plt.xscale('log')
# plt.show()
# plt.xscale('log')
# plt.yscale('log')
# plt.vlines(tOff,0,max(flux))
# plt.show()
# plt.plot(time, text.burnt)
# plt.xscale('log')
# plt.show()
# plt.plot(time, text.Uburn)
# plt.plot(time, text.returnUburn)
# plt.xscale('log')
# plt.show()



'''
                                    plotResults.py
______________________________________________________________________________________________
This script is also an independent module, with the obvious function of plotting results.
See the help(plotResults) for details on the different functions.
'''
import plotResults as pr

# time = np.logspace(-2,20,200)

# Plot the fission yeild plot.
# pr.plotfissionYield(ylog=True)
# Plot full CDF/PDF
# pr.plotFull(time, CDFPDF='CDF', xlog=True, ylog=False)

# CDF = np.genfromtxt('CDF.csv', delimiter=',')

# plt.plot(CDF[:,0], CDF[:,1])
# plt.xscale('log')
# plt.show()
# Plot each individual chain
# pr.plotEachContributingChain(time, CDFPDF = 'CDF', xlog=True)

# Plot stacked chain concentrations for a single decay chain of 1 atom.
# Get needed data (in this example it is Xenon 135)
# atomicNumbers = naming.readableToSeperate('135I-53')
# readableDecayChain, isoDecayData, yeild, fatherDecayConst, fatherIndex = fc.find_by_atomic(atomicNumbers[1], atomicNumbers[0], atomicNumbers[2])
# pr.plotChainConcentrationStacked(int(isoDecayData['NU isomer']), float(fatherDecayConst), float(yeild), isoDecayData, time, xlog=True)

