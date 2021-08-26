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

class bateman():
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
        '''
        self.fatherIsomer = naming.readable(fatherIsomer)
        self.fatherDecayConstant = fatherDecayConstant
        self.yeild = fatherYeild
        self.rateWeight = weighting*self.yeild
        self.decayConsts = decayConsts
        self.decayConsts.insert(0, self.fatherDecayConstant)
        self.decayModes = decayModes
        # self.decayChain = decayChain
    
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

    def CDF_gen(self, show_stacked = False, xlog = False, ylog = False, decayChain = [], time = []):
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
        @param time          : array
                        Add a time array for plotting the stacked concentration plots

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
        # InitCond not used
        initCond = np.zeros(len(self.decayConsts))
        initCond[0] = 1
        #This will store the plots for each species. Not in use in current build
        species = np.zeros(len(self.decayConsts)).tolist()
        solutionPlot = np.zeros(len(species)).tolist()

        def expSum(dcList):
            """
            ExpSum generates a string representation of the bateman
            solution for nuclide j in a linear decay chain. This
            string can be evaluated to generate the output or it
            can be appended to futher strings to contribute to the
            full CDF.

            @param dcList: decay chain list from nuclide 1 to j.

            @return template: string of the executable function.
            """

            denom = np.product([dcList[m] - dcList[0] for m in range(len(dcList)) if m != 0])   #Calculate the first
                                                                                                # denominator
            template = 'np.exp(-' + str(dcList[0]) + '*t)/(' + str(denom) + ')'                 # Begin the template

            for i in range(len(dcList) - 1):    # Work through each decay constant.

                denom = np.product([dcList[m] - dcList[i + 1] for m in range(len(dcList)) if m != i + 1]) # Denominator

                template = template.join(['', '+np.exp(-' + str(dcList[i + 1]) + '*t)/(' + str(denom) + ')']) # append
                                                                                                            # Contribution

            starting_prod = np.product([dcList[i] for i in range(len(dcList) - 1)]) # Calculate the fist product

            template = str(starting_prod) + '*(' + template + ')'   # Append this to the template

            return template


        # Main loop
        CDF_contrib = str(self.nVector[0]) + '*' + expSum(self.decayConsts[:1]) # Generate the first contribution with neutrino scaling

        for i, dc in enumerate(self.decayConsts):

            if i > 0:

                CDF_contrib = CDF_contrib + '+' + str(self.nVector[i]) + '*' + expSum(self.decayConsts[:i + 1])

        self.CDF_contrib = str(self.rateWeight) + '*(' + CDF_contrib + ')'


        return self

    def CDF_gen(self, show_stacked = False, xlog = False, ylog = False, decayChain = [], time = []):
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
        @param time          : array
                        Add a time array for plotting the stacked concentration plots

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
        # InitCond not used
        initCond = np.zeros(len(self.decayConsts))
        initCond[0] = 1
        #This will store the plots for each species. Not in use in current build
        species = np.zeros(len(self.decayConsts)).tolist()
        solutionPlot = np.zeros(len(species)).tolist()

        def expSum(dcList):
            """
            ExpSum generates a string representation of the bateman
            solution for nuclide j in a linear decay chain. This
            string can be evaluated to generate the output or it
            can be appended to futher strings to contribute to the
            full CDF.

            @param dcList: decay chain list from nuclide 1 to j.

            @return template: string of the executable function.
            """

            denom = np.product([dcList[m] - dcList[0] for m in range(len(dcList)) if m != 0])   #Calculate the first
                                                                                                # denominator
            template = 'np.exp(-' + str(dcList[0]) + '*t)/(' + str(denom) + ')'                 # Begin the template

            if dcList[0] != 0:

                int_1_temp = str(-1/dcList[0]) + '*np.exp(-' + str(dcList[0]) + '*t)/(' + str(denom) + ')' # Used for the
                                                                                                    # first integration of G

            else:
                int_1_temp = 't/(' + str(denom) + ')'

            for i in range(len(dcList) - 1):    # Work through each decay constant.

                denom = np.product([dcList[m] - dcList[i + 1] for m in range(len(dcList)) if m != i + 1]) # Denominator

                template = template.join(['', '+np.exp(-' + str(dcList[i + 1]) + '*t)/(' + str(denom) + ')']) # append
                                                                                                            # Contribution
                if dcList[i+1] != 0:

                    int_1_temp = int_1_temp + '+'+str(-1/dcList[i + 1])+'*np.exp(-' + str(dcList[i + 1]) + '*t)/(' + str(denom) + ')'


                else:

                    int_1_temp = int_1_temp + '+t/(' + str(denom) + ')'


            starting_prod = np.product([dcList[i] for i in range(len(dcList) - 1)]) # Calculate the fist product

            template = str(starting_prod) + '*(' + template + ')'   # Append this to the template

            int_1_temp = str(starting_prod) + '*(' + int_1_temp + ')'

            #Finding integration constant

            # G_int_t0 = eval('lambda t:' + int_1_temp)(0)
            #
            # A = -G_int_t0
            #
            # int_1_temp = int_1_temp + '+' + str(A)

            return template, int_1_temp


        # Main loop
        CDF_contrib = '0' # Generate the first contribution with neutrino scaling

        int_1_CDF_contrib = '0'

        for i, dc in enumerate(self.decayConsts):

            if i > 0 and self.nVector[i] != 0:

                CDF_contrib = CDF_contrib + '+' + str(self.nVector[i]) + '*' + expSum(self.decayConsts[:i + 1])[0]

                int_1_CDF_contrib = int_1_CDF_contrib + '+' + str(self.nVector[i]) + '*' + \
                                    expSum(self.decayConsts[:i + 1])[1]

        if self.rateWeight != 0:
            self.CDF_contrib = str(self.rateWeight) + '*(' + CDF_contrib + ')'

            self.int_1_CDF_contrib = str(self.rateWeight) + '*(' + int_1_CDF_contrib + ')'

        else:

            self.CDF_contrib = '0'

            self.int_1_CDF_contrib = '0'


        return self



    def PDF_gen(self):
        '''
        Generating the probility distribution function for neutrinos.

        generates
        ------------
        self.PDF
        self.weightedPDF

        '''
        #         n           n-1       n         n
        # Nn(t) = Σ [Ni(0) x ( ∏ λj) x (Σ e^-λjt/(∏ (λp-λj)))] 
        #        i=1          j=i      j=i     p=i, p≠j
        initCond = np.zeros(len(self.decayConsts))
        initCond[0] = 1
        #This will store the plots for each species
        species = np.zeros(len(self.decayConsts)).tolist()
        self.PDF = 0
        solutionPlot = np.zeros(len(species)).tolist()
        for n,_ in enumerate(species):
            terms = np.zeros(n+1).tolist()
            num_of_exp = [len(terms)-i for i in range(len(terms))]
            for j in range(len(terms)):
                dc = [self.decayConsts[k+j] for k in range(num_of_exp[j])]
                firstProd = np.product(dc[:len(terms)-j-1])
                expVal = [(-dc[k]*np.exp(-dc[k]*self.time)) for k in range(num_of_exp[j])]
                denom = [[] for k in range(num_of_exp[j])]
                denom_dc = [[dc[m] for m in range(k)]+[dc[m] for m in range(k+1,num_of_exp[j])] for k in range(num_of_exp[j])]
                for k in range(len(dc)):
                    if len(denom_dc[k])>0:
                        denom[k] = np.array(denom_dc[k])-dc[k]
                    denom[k] = np.product(denom[k])
                allVals = [expVal[k]/denom[k] for k in range(len(expVal))]
                secondSum = 0
                for k in range(len(allVals)):
                    secondSum = secondSum+allVals[k]
                terms[j] = secondSum*firstProd*initCond[j]
            finalSolution = 0
            for j in range(len(terms)):
                finalSolution = finalSolution + terms[j]
            
            solutionPlot[n] = finalSolution
            
            finalSolution = np.array(finalSolution)*self.nVector[n]
            self.PDF = self.PDF + finalSolution
        self.weightedPDF = np.array(self.PDF)*self.rateWeight
        return self

