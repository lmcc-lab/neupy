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
    def __init__(self, fatherIsomer, fatherDecayConstant, fatherYeild, weighting, decayConsts, decayModes, time):
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
        self.time = time
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
        show_stacked  : booleen
                        shows stacked plot of chain concentrations over time, default to False
        xlog          : booleen
                        for the stacked plot, default is False
        ylog          : booleen
                        for the stacked plot, default is False
        decayChain    : list
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
        #This will store the plots for each species
        species = np.zeros(len(self.decayConsts)).tolist()
        self.CDF = 0
        solutionPlot = np.zeros(len(species)).tolist()
        for n,_ in enumerate(species):
            terms = np.zeros(n+1).tolist()
            num_of_exp = [len(terms)-i for i in range(len(terms))]
            for j in range(len(terms)):
                dc = [self.decayConsts[k+j] for k in range(num_of_exp[j])]
                firstProd = np.product(dc[:len(terms)-j-1])
                expVal = [np.exp(-dc[k]*self.time) for k in range(num_of_exp[j])]
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
            self.CDF = self.CDF + finalSolution
        self.weightedCDF = np.array(self.CDF)*self.rateWeight
        if show_stacked == True:
            if len(solutionPlot) == 6:
                plt.stackplot(self.time, solutionPlot[0], solutionPlot[1], solutionPlot[2], solutionPlot[3], solutionPlot[4], solutionPlot[5], labels=decayChain)
            elif len(solutionPlot) == 5:
                plt.stackplot(self.time, solutionPlot[0], solutionPlot[1], solutionPlot[2], solutionPlot[3], solutionPlot[4], labels=decayChain)
            elif len(solutionPlot) == 4:
                plt.stackplot(self.time, solutionPlot[0], solutionPlot[1], solutionPlot[2], solutionPlot[3], labels=decayChain)
            elif len(solutionPlot) == 3:
                plt.stackplot(self.time, solutionPlot[0], solutionPlot[1], solutionPlot[2], labels=decayChain)
            elif len(solutionPlot) == 2:
                plt.stackplot(self.time, solutionPlot[0], solutionPlot[1], labels=decayChain)
            elif len(solutionPlot) == 1:
                plt.stackplot(self.time, solutionPlot[0], labels=decayChain)
            else:
                plt.plot(self.time, np.zeros(len(self.time)))

            if xlog == True:
                plt.xscale('log')
            if ylog == True:
                plt.yscale('log')
            plt.xlabel('time (s)')
            plt.legend(loc='upper right')
            plt.ylabel('Concetration')
            plt.title('Fission isotope '+str(self.fatherIsomer))
            plt.show()
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

