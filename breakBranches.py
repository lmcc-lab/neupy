'''
Author: Liam McClelland

Last edited: 31/05/2021
'''


from dependencies import split, counter
import pandas as pd
import numpy as np
import preferences
import math

class branches():
	'''
	breakBranches module deals with seperating all of the decay branches into 
	single, easy to manage branches with associated weightings.

								   father
								   	|
								   daughter
									/   \
								   d2	 d1
								  / |	  |
								d3 d4     d5

	branch1 = [daughter, d2, d3]  with associated weighting = fatheryeild*br(daughter)*br(d2)*br(d3)
	branch2 = [daughter, d2, d4]  with associated weighting = fatheryeild*br(daughter)*br(d2)*br(d4)
	...
	These branchesa are embedded as [[branch1],[branch 2],[branch3]...]
	'''
	def __init__(self, decayChain, chainIndex, branchRatios, decayModes, decayConsts, branchRatioUncerts, decayConstsUncerts):
		'''
		parameters
		--------------
		decayChain, list - full decayChain with NU formatted elements.
		chainIndex, list - associated indexes for tracking individual branches with int elements
		branchRatio, list - associated branch ratios
		decayConsts, list - associated decay constants.
		'''
		self.decayChain = decayChain
		self.chainIndex = chainIndex
		self.branchRatios = branchRatios
		self.decayModes = decayModes
		self.decayConsts = decayConsts
		self.branchRatioUncerts = branchRatioUncerts
		self.decayConstsUncerts = decayConstsUncerts
		
	def breakBranches(self):
		'''
		Seperate each branch into seperate chain contained within a single list object.
		Data is kept in seperate but compareable lists.

		generates
		-------------
		self.brokenDecayChainBranches
		self.brokenBranchRatios
		self.brokenDecayModes
		self.brokenDecayConsts
		self.weighting
		'''
		i = -2
		duplicateChain = self.decayChain
		duplicateIndex = self.chainIndex
		duplicateBranchRatio = self.branchRatios
		duplicateDecayModes = self.decayModes
		duplicateDecayConsts = self.decayConsts
		duplicateBranchRatioUncerts = self.branchRatioUncerts
		duplicateDecayConstsUncerts = self.decayConstsUncerts
		self.brokenDecayChainBranches = []
		self.brokenBranchRatios = []
		self.brokenDecayModes = []
		self.brokenDecayConsts = []
		self.brokenBranchRatiosUncerts = []
		self.brokenDecayConstsUncerts = []
		while i*(-1) <= len(duplicateChain[:-1]):
			for j, index in enumerate(duplicateIndex[i]):
				linearChain = [duplicateChain[i][j]]
				linearBranchRatio = [float(duplicateBranchRatio[i][j])/100]
				linearDecayModes = [duplicateDecayModes[i][j]]
				linearDecayConsts = [duplicateDecayConsts[i+1][j]]
				linearDecayConstsUncerts = [duplicateDecayConstsUncerts[i + 1][j]]
				linearBranchRatioUncerts = [float(duplicateBranchRatioUncerts[i][j])/100]
				val = split(index)
				for k,_ in enumerate(val[:-1]):
					prev_gen = len(val[:-1-k])-1
					prev_index = ''.join(val[:-1-k])
					linearChain.append(duplicateChain[prev_gen][duplicateIndex[prev_gen].index(prev_index)])
					linearBranchRatio.append(duplicateBranchRatio[prev_gen][duplicateIndex[prev_gen].index(prev_index)]/100)
					linearDecayModes.append(duplicateDecayModes[prev_gen][duplicateIndex[prev_gen].index(prev_index)])
					linearDecayConsts.append(duplicateDecayConsts[prev_gen][duplicateIndex[prev_gen].index(prev_index)])
					linearDecayConstsUncerts.append(duplicateDecayConstsUncerts[prev_gen][duplicateIndex[prev_gen].index(prev_index)])
					linearBranchRatioUncerts.append(duplicateBranchRatioUncerts[prev_gen][duplicateIndex[prev_gen].index(prev_index)]/100)
				linearChain.reverse()
				linearBranchRatio.reverse()
				linearDecayModes.reverse()
				linearDecayConsts.reverse()
				linearDecayConstsUncerts.reverse()
				linearBranchRatioUncerts.reverse()
				checkNan = [math.isnan(linearBranchRatio[i]) for i in range(len(linearBranchRatio))]
				checkNanUncert = [math.isnan(linearBranchRatioUncerts[i]) for i in range(len(linearBranchRatioUncerts))]
				if True in checkNan:
					linearBranchRatio[checkNan.index(True)] = 0
				if True in checkNanUncert:
					linearBranchRatioUncerts[checkNanUncert.index(True)] = 0
				if len(self.brokenDecayChainBranches)>0:
					for m in range(len(self.brokenDecayChainBranches)):
						str1 = split(str(linearChain))
						str1.remove('[')
						str1.remove(']')
						str1 = ''.join(str1)
						str2 = split(str(self.brokenDecayChainBranches[m]))
						str2.remove('[')
						str2.remove(']')
						str2 = ''.join(str2)
						if str1 in str2:
							break
						elif m == len(self.brokenDecayChainBranches)-1:
							self.brokenDecayChainBranches.append(linearChain)
							self.brokenBranchRatios.append(linearBranchRatio)
							self.brokenDecayModes.append(linearDecayModes)
							self.brokenDecayConsts.append(linearDecayConsts)
							self.brokenDecayConstsUncerts.append(linearDecayConstsUncerts)
							self.brokenBranchRatiosUncerts.append(linearBranchRatioUncerts)
				else:
					self.brokenDecayChainBranches.append(linearChain)
					self.brokenBranchRatios.append(linearBranchRatio)
					self.brokenDecayModes.append(linearDecayModes)
					self.brokenDecayConsts.append(linearDecayConsts)
					self.brokenDecayConstsUncerts.append(linearDecayConstsUncerts)
					self.brokenBranchRatiosUncerts.append(linearBranchRatioUncerts)
			i = i-1
		self.weighting = [np.prod(i) for i in self.brokenBranchRatios]
		self.weightingUncert = [np.sqrt(np.sum([(i[j]/self.weighting[i][j])**2 for j in len(i)])) for i in self.brokenBranchRatiosUncerts]
		return self