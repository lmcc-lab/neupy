import nubaseSorter as sort
import numpy as np
import math
import preferences
import time

class chain(): 
    '''
    Chain class is used for generating decay chains based on decay modes. All decay chain data including branch
    ratios and decay constants are saved as part of the object.

    '''
    def __init__(self, daughters, branchRatios, decayModes, branchRatiosUncerts):
        '''
        @param daughters, list - All daughters must be in NU format.
        @param branchRatios, list - Containing float values
        @param decayModes, list - Contains strings.

        Chain is of the form [[Gen1], [Gen2], [Gen3]] such that
        len(Chain)=number chains and each Gen contains
        [daughter1, daughter2, daughter3...]

        Keeps track of chain depth and compares to preferences.max_chain_depth to see if it is maxed out.

        generates:  self.isomerChain - list with embedded daughters
                    self.branchRatioChain - list with embedded daughter branch ratios
                    self.decayModeChain - list with embedded daughter deacay modes
                    self.chainIndex - list, indexes all daughters for navagating branches
                    self.decayConstChain - list with embedded decay constants.
                    self.Missing - missing isomer list
                    self.daughters - input daughters
                    self.branches - length of input daughters
                    self.chainDepthExceeded - int, keeps track of chain depth.

                    Uncertainties
                    self.branchRatioUncerts
                    self.decayConstUncerts
        '''
        self.isomerChain = []
        self.branchRatioChain = []
        self.decayModeChain = []
        self.chainIndex = []
        self.decayConstChain = []
        self.Missing = []
        # Uncertainties
        self.branchRatioUncertsChain = []
        self.decayConstUncertsChain = []
        self.daughters = daughters
        self.branches = len(self.daughters)
        self.isomerChain.append(daughters) #Chain with generations imbedded. [[]]
        self.decayModeChain.append(decayModes)

        for j, br in enumerate(branchRatios):
            if preferences.ignore_energetically_possible == False:
                if math.isnan(br) == True:
                    if sum([sum(branchRatios[:j]),sum(branchRatios[j+1:])])>100:
                        branchRatios[j] = 0
                    else:
                        branchRatios[j] = 100-sum([sum(branchRatios[:j]),sum(branchRatios[j+1:])])

        self.branchRatioChain.append(branchRatios)
        self.branchRatioUncertsChain.append(branchRatiosUncerts)
        branchIndex = np.zeros(len(self.isomerChain[-1])).tolist()
        for i,_ in enumerate(self.isomerChain[-1]):
            branchIndex[i] = str(i)
        self.chainIndex.append(branchIndex)
        self.chainDepthExceeded = 0


    
    def chain_gen(self):
        '''
        Generates the next generation of daughters by going through the list of last generation daughters, finding their decay modes
        from nubase and adding it to the daughter object. It also saves all relevant data such as branch ratios, indexes, decay modes
        etc. 

        To save processing time, one can change the max_decay_chain_depth in preferences.py. If the chain reaches the max decay chain
        depth then the process of finding the next generation will stop and a message suggesting an increase of preference will be
        made. The number of times this maximum chain is reached is stored in self.chainDepthExceeded.
        '''
        currentGen = len(self.isomerChain)
        if currentGen < preferences.max_decay_chain_depth:
            self.isomerChain.append([])
            self.branchRatioChain.append([])
            self.decayModeChain.append([])
            self.chainIndex.append([])
            self.decayConstChain.append([])
            self.decayConstUncertsChain.append([])
            self.branchRatioUncertsChain.append([])
            for i, iso in enumerate(self.isomerChain[currentGen-1]):
                branchCode = self.chainIndex[currentGen-1][i]
                daughterIsomer = sort.isomer(iso//10000000,(iso//10000)%1000,iso%10000, 0, 0)
                try:
                    daughterIsomer.find_nubase()
                    daughterIsomer.excit_states()
                    daughterIsomer.half_life()
                    daughterIsomer.decay_info()
                    daughterIsomer.decay_constant()
                    daughterIsomer.daughter()
                    self.decayConstChain[currentGen-1].append(daughterIsomer.decay_const)
                    self.decayConstUncertsChain[currentGen - 1].append(daughterIsomer.decay_const_uncert)
                    for j, D in enumerate(daughterIsomer.daughters):
                        if preferences.ignore_energetically_possible == True:
                            if daughterIsomer.branchRatio[j] > preferences.min_branch_ratio and math.isnan(daughterIsomer.branchRatio[j]) == False:
                                self.isomerChain[currentGen].append(D)
                                self.branchRatioChain[currentGen].append(daughterIsomer.branchRatio[j])
                                self.branchRatioUncertsChain[currentGen].append(daughterIsomer.branchRatioUncert[j])
                                self.decayModeChain[currentGen].append(daughterIsomer.decayModes[j])
                                self.chainIndex[currentGen].append(str(branchCode)+str(j))
                        else:
                            if daughterIsomer.branchRatio[j] > preferences.min_branch_ratio or math.isnan(daughterIsomer.branchRatio[j]) == True:
                                if math.isnan(daughterIsomer.branchRatio[j]) == True:
                                    if sum([sum(daughterIsomer.branchRatio[:j]),sum(daughterIsomer.branchRatio[j+1:])])>100:
                                        daughterIsomer.branchRatio[j] = 0
                                    else:
                                        daughterIsomer.branchRatio[j] = 100-sum([sum(daughterIsomer.branchRatio[:j]),sum(daughterIsomer.branchRatio[j+1:])])
                                self.isomerChain[currentGen].append(D)
                                self.branchRatioChain[currentGen].append(daughterIsomer.branchRatio[j])
                                self.branchRatioUncertsChain[currentGen].append(daughterIsomer.branchRatioUncert[j])
                                self.decayModeChain[currentGen].append(daughterIsomer.decayModes[j])
                                self.chainIndex[currentGen].append(str(branchCode)+str(j))
                except ValueError:
                    daughterIsomer.missing_isomer()
                    self.Missing.append(daughterIsomer.NU)
        try:
            if currentGen < preferences.max_decay_chain_depth:
                if len(self.isomerChain[currentGen]) > 0:
                    self.chain_gen()
        except AttributeError:
            pass
        
        self.decayDepth = len(self.isomerChain)
        if self.decayDepth == preferences.max_decay_chain_depth:
            #raise 'Maxmimum decay chain depth reached. Consider increasing maximum decay chain depth.'
            self.chainDepthExceeded += 1
            print('Maxmimum decay chain depth reached. Consider increasing maximum decay chain depth.')
        return self
                  

    