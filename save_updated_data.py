from typing import final
import pandas as pd
import numpy as np
from dependencies import conv_str_to_list as ctl, find_index
import preferences
import naming
import constants
import breakBranches as bb
import bateman
import os
import sys

path = '.'
outdir=path+'/Contributing_chains/'
if (os.access(outdir,os.F_OK))==0:
    os.mkdir(outdir)

if (os.access('./cache/', os.F_OK))==0:
    os.mkdir('./cache/')


def checkSimple(unBrokendecayModes):
    '''
    Checks if the current chain is a simple, linear chain
    '''
    for _, gen in enumerate(unBrokendecayModes):
        if len(gen)>1:
            linearFlag = False
            break
        else:
            linearFlag = True
    for _, gen in enumerate(unBrokendecayModes):
        if linearFlag == True:
            if len(gen)>0:
                if gen[0] == 'b-':
                    simpleFlag = True
                else:
                    simpleFlag = False
                    break
            else:
                simpleFlag = True
        else:
            simpleFlag = False
    return simpleFlag

class export_data():
    '''
    Export data class deals with saving all known CDF, PDF and contribution
    data for speeding up calculation times.
    '''      
    def __init__(self, time):
        '''
        Gets isotope data and chain data from decy_chain_data folder.

        parameters
        ------------------
        time - time array
        '''
        self.time = time
        self.isotopeData = pd.read_csv('./decay_chain_data/IgnoreTheoreticEnergy='+str(preferences.ignore_energetically_possible)+'_isomerData.csv')
        self.chainData = pd.read_csv('./decay_chain_data/IgnoreTheoreticEnergy='+str(preferences.ignore_energetically_possible)+'_decayChainData.csv')
        self.chainCounter = 0
        self.father_counter = 0
        # Cache previous time values
        with open('./cache/timeValues.txt', 'w+') as f:
            f.write(str(time[0])+'\n'+str(time[-1])+'\n'+str(len(time)))

    def contributers(self):
        #initilalising the progress bar
        toolbar_width = 40
        sys.stdout.write("[%s]" % (" "*toolbar_width))
        sys.stdout.flush()
        sys.stdout.write("\b"*(toolbar_width+1))
        '''
        Output all individual chain contributions including max CDF contribution for
        the purpose of reducing the ~2500 chains for more direct analysis without
        loosing the majority of the information

        exports
        ----------------
        *simpleTitle*ConstributingChainsMaxCDF.csv
        Where simpleTitle is either Simple or NotSimple depending on if preferences.simple
        is equal to True or False respectively.
        '''
        progression = np.linspace(0,self.isotopeData.shape[0],toolbar_width)
        progression = [int(np.round(progression[i])) for i in range(toolbar_width)]
        for index, _ in self.isotopeData.head(n=self.isotopeData.shape[0]).iterrows():
            distribution_contributions = {'Father':[], 'Chain':[], 'Chain weigthing':[], 'Weighted CDF solution':[],'Max CDF contrib':[], 'Weighted PDF solution':[], 'PDF Area contrib':[]}
            distribution_contributions = pd.DataFrame(distribution_contributions, columns= ['Father','Chain', 'Chain weigthing','Weighted CDF solution','Max CDF contrib','Weighted PDF solution','PDF Area contrib'])
            if index in progression:
                sys.stdout.write("-")
            self.father_counter += 1
            # print(run_bateman.counter,'/',numEntries)
            fatherIsomer = int(self.isotopeData['NU isomer'][index])
            fatherDecayConst = float(self.isotopeData['Decay Constant'][index])
            fatherYeild = float(self.isotopeData['Yield'][index])
            decayChain = ctl(self.chainData['Decay Chain'][index])
            branchRatios = ctl(self.chainData['Branch Ratio'][index])
            chainIndex = ctl(self.chainData['Branch index'][index])
            decayConstChain = ctl(self.chainData['Decay Constant'][index])
            decayModeChain = ctl(self.chainData['Decay Mode'][index])
            simpleFlag = checkSimple(decayModeChain)
            if preferences.simple == True and simpleFlag == True:
                brokenBranches = bb.branches(decayChain, chainIndex, branchRatios, decayModeChain, decayConstChain)
                brokenBranches.breakBranches()
                for i, chain in enumerate(brokenBranches.brokenDecayModes):
                    singleCDF = bateman.bateman(fatherIsomer, fatherDecayConst, fatherYeild, brokenBranches.weighting[i], brokenBranches.brokenDecayConsts[i], brokenBranches.brokenDecayModes[i], self.time)
                    singleCDF.neutrinoVector()
                    singleCDF.CDF_gen()
                    singleCDF.PDF_gen()
                    CDF_max = singleCDF.weightedCDF[-1]
                    PDF_area = np.trapz(singleCDF.weightedPDF, self.time)
                    self.chainCounter +=1
                    readableChain = np.zeros(len(chain)).tolist()
                    for j, daught in enumerate(brokenBranches.brokenDecayChainBranches[i]):
                        readableChain[j] = naming.readable(daught)
                    distributions = pd.DataFrame([[naming.readable(fatherIsomer),readableChain, brokenBranches.weighting[i]*fatherYeild, list(singleCDF.weightedCDF), CDF_max, list(singleCDF.weightedPDF), PDF_area]],columns= ['Father','Chain','Chain weigthing','Weighted CDF solution','Max CDF contrib','Weighted PDF solution','PDF Area contrib'])
                    distribution_contributions = distribution_contributions.append(distributions)
            elif preferences.simple == True and simpleFlag == False:
                pass
            elif preferences.simple == False:
                brokenBranches = bb.branches(decayChain, chainIndex, branchRatios, decayModeChain, decayConstChain)
                brokenBranches.breakBranches()
                for i, chain in enumerate(brokenBranches.brokenDecayModes):
                    singleCDF = bateman.bateman(fatherIsomer, fatherDecayConst, fatherYeild, brokenBranches.weighting[i], brokenBranches.brokenDecayConsts[i], chain, self.time)
                    singleCDF.neutrinoVector()
                    singleCDF.CDF_gen()
                    singleCDF.PDF_gen()
                    CDF_max = singleCDF.weightedCDF[-1]
                    PDF_area = np.trapz(singleCDF.weightedPDF, self.time)
                    self.chainCounter +=1
                    readableChain = np.zeros(len(chain)).tolist()
                    for j, daught in enumerate(brokenBranches.brokenDecayChainBranches[i]):
                        readableChain[j] = naming.readable(daught)
                    distributions = pd.DataFrame([[naming.readable(fatherIsomer),readableChain,brokenBranches.weighting[i]*fatherYeild, list(singleCDF.weightedCDF), CDF_max, list(singleCDF.weightedPDF), PDF_area]],columns= ['Father','Chain','Chain weigthing','Weighted CDF solution','Max CDF contrib','Weighted PDF solution','PDF Area contrib'])
                    distribution_contributions = distribution_contributions.append(distributions)
            sys.stdout.flush()
            print('Ready to output')
            if index == 0:
                distribution_contributions.to_csv(outdir+preferences.simpleTitle+'ConstributingChainsMaxCDF.csv', sep=',')
                print('Output for i=0')
            else:
                distribution_contributions.to_csv(outdir+preferences.simpleTitle+'ConstributingChainsMaxCDF.csv', mode='a', header=False, sep=',')
                print('Output for i=', index)
        print('\nExport complete')
        return self
    

    def maxContributions(self, full=False):
        '''
        This method deals with returning and adding the CDF contributions based on
        the top contributers. This will stop when the percentage of the CDF meets
        preferences.percent_of_final, if preferences.pref_percent is True, or it
        will stop if the number of chains used is greater than preferences.top_num
        if preferences.pref_percent is False.

        Parameters
        ----------------------------
        full - optional, default is False. If True then it will run through all
        contributions and output a CDF_PDF_full.csv file. The returned data will
        all be False.

        Generates
        ----------------------------
        CDF_PDF_full.csv 
        *simpleTitle*SelectContributingChains.csv

        Returns
        ----------------------------
        running_CDF - The CDF before it cuts off at maximum contribtions
        running_PDF - The PDF before it cuts off at maximum contribtions
        num_contrib - The number of contributing chains for these outputs
        running_contributions - The raw chain data for these contributions
        running_percent - The percentage of the contributing CDF to the full CDF.
        '''
        final_CDF = 0
        final_PDF = 0
        Contributing_data = pd.read_csv('./Contributing_chains/'+preferences.simpleTitle+'ConstributingChainsMaxCDF.csv')
        if full==True:
            for i, _ in Contributing_data.head(n=Contributing_data.shape[0]).iterrows():
                final_CDF = final_CDF + np.array(ctl(Contributing_data.iloc[i,4]))
                final_PDF = final_PDF + np.array(ctl(Contributing_data.iloc[i,6]))
            output_CDF = pd.DataFrame({'CDF':[], 'PDF':[]}, columns = ['CDF','PDF'])
            final = pd.DataFrame([[list(final_CDF),list(final_PDF)]], columns = ['CDF','PDF'])
            output_CDF = output_CDF.append(final)
            output_CDF.to_csv(outdir+preferences.simpleTitle+'CDF_PDF_full.csv',sep=',')
            running_CDF = False
            running_PDF = False
            num_contrib = False
            running_contributions = False
            running_percent = False
        else:
            Contributing_data = Contributing_data.sort_values(by='Max CDF contrib', ascending=False)
            Contributing_data = Contributing_data.reset_index()
            sumMax = sum(Contributing_data['Max CDF contrib'])
            running_max = 0
            running_CDF, running_PDF = 0,0
            for i, _ in Contributing_data.head(n=Contributing_data.shape[0]).iterrows():
                running_max += Contributing_data.iloc[i,5]
                running_percent = running_max/sumMax*100
                running_CDF = running_CDF + np.array(ctl(Contributing_data.iloc[i,4]))
                running_PDF = running_PDF + np.array(ctl(Contributing_data.iloc[i,6]))
                if running_percent >= preferences.percent_of_final and preferences.pref_percent == True:
                    running_contributions = Contributing_data[:i]
                    num_contrib = i
                    break
                elif i+1 == preferences.top_num and preferences.pref_percent == False:
                    running_contributions = Contributing_data[:i]
                    num_contrib = i
                    break
            running_contributions.to_csv(outdir+preferences.simpleTitle+'SelectContributingChains.csv',sep=',')
        print('Data saving complete')
        return running_CDF, running_PDF, num_contrib, running_contributions, running_percent

    def full_CDF(self, save = False):

        toolbar_width = 40
        sys.stdout.write("[%s]" % (" " * toolbar_width))
        sys.stdout.flush()
        sys.stdout.write("\b" * (toolbar_width + 1))
        progression = np.linspace(0, self.isotopeData.shape[0], toolbar_width)
        progression = [int(np.round(progression[i])) for i in range(toolbar_width)]
        final_CDF = 0
        for index, _ in self.isotopeData.head(n=self.isotopeData.shape[0]).iterrows():
            if index in progression:
                sys.stdout.write("-")

            fatherIsomer = int(self.isotopeData['NU isomer'][index])
            fatherDecayConst = float(self.isotopeData['Decay Constant'][index])
            fatherYeild = float(self.isotopeData['Yield'][index])
            decayChain = ctl(self.chainData['Decay Chain'][index])
            branchRatios = ctl(self.chainData['Branch Ratio'][index])
            chainIndex = ctl(self.chainData['Branch index'][index])
            decayConstChain = ctl(self.chainData['Decay Constant'][index])
            decayModeChain = ctl(self.chainData['Decay Mode'][index])
            simpleFlag = checkSimple(decayModeChain)
            if preferences.simple == True and simpleFlag == True:
                brokenBranches = bb.branches(decayChain, chainIndex, branchRatios, decayModeChain, decayConstChain)
                brokenBranches.breakBranches()
                for i, chain in enumerate(brokenBranches.brokenDecayModes):
                    singleCDF = bateman.bateman(fatherIsomer, fatherDecayConst, fatherYeild, brokenBranches.weighting[i], brokenBranches.brokenDecayConsts[i], brokenBranches.brokenDecayModes[i], self.time)
                    singleCDF.neutrinoVector()
                    singleCDF.CDF_gen()
                    final_CDF = final_CDF + singleCDF.weightedCDF
            elif preferences.simple == True and simpleFlag == False:
                pass
            elif preferences.simple == False:
                brokenBranches = bb.branches(decayChain, chainIndex, branchRatios, decayModeChain, decayConstChain)
                brokenBranches.breakBranches()
                for i, chain in enumerate(brokenBranches.brokenDecayModes):
                    singleCDF = bateman.bateman(fatherIsomer, fatherDecayConst, fatherYeild, brokenBranches.weighting[i], brokenBranches.brokenDecayConsts[i], chain, self.time)
                    singleCDF.neutrinoVector()
                    singleCDF.CDF_gen()
                    final_CDF = final_CDF + singleCDF.weightedCDF
            sys.stdout.flush()
        CDF_output = np.hstack((np.array([self.time]).T, np.array([final_CDF]).T))
        if save == True:
            np.savetxt('CDF.csv',CDF_output, delimiter=',')
            print('\nExport complete')
            
        return CDF_output
