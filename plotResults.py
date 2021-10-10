from nubaseSorter import A, Yeild, Yeild_uncert
import matplotlib.pyplot as plt
import save_updated_data as sud
import preferences
import pandas as pd
from dependencies import conv_str_to_list as ctl
import bateman
import breakBranches as bb

def plotfissionYield(ylog = True):
	allData = pd.read_csv('./decay_chain_data/IgnoreTheoreticEnergy='+str(preferences.ignore_energetically_possible)+'_isomerData.csv')
	A = [(allData['NU isomer'][i]//10000)%1000 for i in range(len(allData['NU isomer']))]
	Yeild = allData['Yield']
	Yeild_uncert = allData['Yield uncertainty']
	plt.errorbar(A, Yeild, Yeild_uncert,fmt = 'o', markersize = 2, elinewidth = 1, ecolor = 'k', color='k')
	if ylog == True:
		plt.yscale('log')
	plt.xlabel('A')
	plt.ylabel('Fission Yield')
	plt.grid()
	plt.show()

def plotFull(time, CDFPDF = 'CDF', xlog = False, ylog = False):
	try:
		CDFPDF_data = pd.read_csv('./Contributing_chains/'+preferences.simpleTitle+'CDF_PDF_full.csv')
		if CDFPDF == 'CDF':
			ploting = ctl(CDFPDF_data.iloc[0,1])
		else:
			ploting = ctl(CDFPDF_data.iloc[0,2])
	except FileNotFoundError:
		sud.export_data(time).contributers().maxContributions(full=True)
		CDFPDF_data = pd.read_csv('./Contributing_chains/'+preferences.simpleTitle+'CDF_PDF_full.csv')
		if CDFPDF == 'CDF':
			ploting = ctl(CDFPDF_data.iloc[0,1])
		else:
			ploting = ctl(CDFPDF_data.iloc[0,2])
	with open('./cache/timeValues.txt', 'r') as f:
		lines = f.readlines()
		start = float(lines[0])
		stop = float(lines[1])
		length = int(lines[2])
	if float(time[0]) != start or float(time[-1]) != stop or len(time) != length:
		print("Time has changed, updating data with new time array.")
		sud.export_data(time).contributers().maxContributions(full=True)
		CDFPDF_data = pd.read_csv('./Contributing_chains/'+preferences.simpleTitle+'CDF_PDF_full.csv')
		if CDFPDF == 'CDF':
			ploting = ctl(CDFPDF_data.iloc[0,1])
		else:
			ploting = ctl(CDFPDF_data.iloc[0,2])
	plt.plot(time, ploting)
	if xlog == True:
		plt.xscale('log')
	if ylog == True:
		plt.yscale('log')
	if CDFPDF == 'CDF':
		plt.ylabel('Cummulitive neutrino emissions')
	else:
		plt.ylabel('Probability of neutrino emission')
	plt.xlabel('Time (s)')
	plt.grid()
	plt.show()


def plotEachContributingChain(time, CDFPDF = 'CDF', xlog = False, ylog = False):
	try:
		contributions = pd.read_csv('./Contributing_chains/'+preferences.simpleTitle+'ConstributingChainsMaxCDF.csv')
		if CDFPDF == 'CDF':
			ploting = contributions.iloc[:,4]
		else:
			ploting = contributions.iloc[:,6]
	except FileNotFoundError:
		sud.export_data(time).contributers()
		contributions = pd.read_csv('./Contributing_chains/'+preferences.simpleTitle+'ConstributingChainsMaxCDF.csv')
		if CDFPDF == 'CDF':
			ploting = contributions.iloc[:,4]
		else:
			ploting = contributions.iloc[:,6]
	with open('./cache/timeValues.txt', 'r') as f:
		lines = f.readlines()
		start = float(lines[0])
		stop = float(lines[1])
		length = int(lines[2])
	if float(time[0]) != start or float(time[-1]) != stop or len(time) != length:
		print("Time has changed, updating data with new time array.")
		sud.export_data(time).contributers()
		contributions = pd.read_csv('./Contributing_chains/'+preferences.simpleTitle+'ConstributingChainsMaxCDF.csv')
		if CDFPDF == 'CDF':
			ploting = contributions.iloc[:,4]
		else:
			ploting = contributions.iloc[:,6]
	for row in ploting:
		plt.plot(time, ctl(row))
	if xlog == True:
			plt.xscale('log')
	if ylog == True:
		plt.yscale('log')
	if CDFPDF == 'CDF':
		plt.ylabel('Cummulitive neutrino emissions')
	else:
		plt.ylabel('Probability of neutrino emission')
	plt.xlabel('Time (s)')
	plt.grid()
	plt.show()

def plotChainConcentrationStacked(fatherIsomer, fatherDecayConstant, fatherYeild, isoDecayData, time, xlog = False, ylog = False):
	brokenBranches = bb.branches(ctl(isoDecayData['Decay Chain']), ctl(isoDecayData['Branch index']), ctl(isoDecayData['Branch Ratio']), ctl(isoDecayData['Decay Mode']), ctl(isoDecayData['Decay Constant'])).breakBranches()
	for i, dc in enumerate(brokenBranches.brokenDecayConsts):
		bateman.bateman(fatherIsomer, fatherDecayConstant, fatherYeild, brokenBranches.weighting[i], dc, brokenBranches.brokenDecayModes[i], time).neutrinoVector().CDF_gen(show_stacked=True, xlog=xlog, ylog=ylog, decayChain = brokenBranches.brokenDecayChainBranches[i])
	