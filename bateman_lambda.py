'''
bateman.py deals with denerating neutrino vectors, CDF's and PDF's of neutrino emissions for a single fission event.

Author: Liam McClelland

Last Edited: 31/05/21

'''


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from dependencies import split
import preferences
import naming
import constants
import math

class bateman_function():
    '''
    bateman class is initialised using the father isomer (this is the fission product), the fathers decay constant,
    the fathers yeild,

    '''
    def __init__(self, fatherIsomer, fatherDecayConstant, fatherYeild, weighting, decayConsts, decayModes):
        '''
        parameters
        -------------
        fatherIsomer, int - NU formated isotope
        fatherDecayConstant, float - father isotopes decay constant
        fatherYeild, float - father isotopes fission yield
        decayChain, list - broken decay chain NU format
        weighting, float - broken chains weighting
        decayConsts, list - float elements, broken chains decay constants
        decayModes, list - string elements, broken chains decay modes
        time, array - time span.
        '''
        self.fatherIsomer = naming.readable(fatherIsomer)
        self.fatherDecayConstant = fatherDecayConstant
        self.yeild = fatherYeild
        self.rateWeight = weighting*self.yeild
        self.decayConsts = decayConsts
        self.decayConsts.insert(0, self.fatherDecayConstant)
        self.decayModes = decayModes

    def neutrinoVector(self):
        '''
        Neutrino vector is a list that encapulates how many neutrinos should be present
        in a decay chain due to there decay mode.

        generates
        --------------

        self.nVector
        '''
        self.nVector = np.zeros(len(self.decayModes)+1).tolist()
        for i,dm in enumerate(self.decayModes):
            if i < len(self.nVector)-1:
                if dm == 'b-':
                    self.nVector[i+1] = self.nVector[i]+1
                elif dm == '2b-':
                    self.nVector[i+1] = self.nVector[i]+2
                elif dm == 'b+':
                    self.nVector[i+1] = self.nVector[i]+1
                elif dm == '2b+':
                    self.nVector[i+1] = self.nVector[i]+2
                elif dm == 'ec':
                    self.nVector[i+1] = self.nVector[i]+1
                else:
                    self.nVector[i+1] = self.nVector[i]
        return self

    def CDF_gen(self, show_stacked = False, xlog = False, ylog = False, decayChain = []):
        '''
        Generating the cumulative distribution
        function for neutrinos.

        optional parameters
        -------------
        @param show_stacked  : booleen
                        shows stacked plot of chain concentrations over time, default to False
        @param xlog          : booleen
                        for the stacked plot, default is False
        @param ylog          : booleen
                        for the stacked plot, default is False
        @param decayChain    : list
                        for the stacked plot legend. Default is empty list []


        generates
        ------------
        self.CDF
        self.weightedCDF

        '''
        #         n           n-1       n         n
        # Nn(t) = Σ [Ni(0) x ( ∏ λj) x (Σ e^-λjt/(∏ (λp-λj)))]
        #        i=1          j=i      j=i     p=i, p≠j
        decayChain = naming.readableChain(decayChain)
        decayChain.insert(0, self.fatherIsomer)
        initCond = np.zeros(len(self.decayConsts))
        initCond[0] = 1

        def expSum(dcList):
            # Generate exp code
            # exp template
            print('Input dcList', dcList)
            denom = np.product([dcList[m] - dcList[0] for m in range(len(dcList)) if m != 0])
            print('Denominator', denom)
            template = 'np.exp(-' + str(dcList[0]) + '*t)/(' + str(denom) + ')'
            for i in range(len(dcList) - 1):
                denom = np.product([dcList[m] - dcList[i + 1] for m in range(len(dcList)) if m != i + 1])
                template = template.join(['', '+np.exp(-' + str(dcList[i + 1]) + '*t)/(' + str(denom) + ')'])
            starting_prod = np.product([dcList[i] for i in range(len(dcList) - 1)])
            # template = ''.join([str(starting_prod)+'*(', template, ')'])
            template = str(starting_prod) + '*(' + template + ')'
            print(template)
            return template


        CDF_contrib = str(self.nVector[0]) + '*' + expSum(self.decayConsts[:1])
        for i, dc in enumerate(self.decayConsts):
            if i > 0:
                CDF_contrib = CDF_contrib + '+' + str(self.nVector[i]) + '*' + expSum(self.decayConsts[:i + 1])
        self.CDF_contrib = str(self.rateWeight) + '*(' + CDF_contrib + ')'

        return self