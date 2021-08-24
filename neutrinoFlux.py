import numpy as np
from numpy.lib.index_tricks import ndenumerate
import constants
import matplotlib.pyplot as plt
import time
import scipy.integrate as integrate
import dependencies
import pandas as pd
import preferences
from matplotlib.animation import FuncAnimation
import save_updated_data as sud
from dependencies import conv_str_to_list as ctl
import dependencies as dp

class neuFlux():

	def __init__(self, time, thermal, CDF):
		'''
		parameters
		----------
		time    :   array
					time axis - must be linearly scalled and have t=0 at reactor turn on.
		thermal :   array
					thermal output function, J/s vs time
		CDF		:	string
					Singular U235 cummulitive neutrino emission
		
		returns
		--------
		neutrinoFlux :  array
						Full neutrino flux for reactor operations.
		'''
		self.time = time
		self.Uburn = 0.94*thermal/constants.U235_fiss
		self.CDF = lambda t: eval(CDF)
	
	def flux(self):
		CummulitiveneutrinoFlux = np.zeros(len(self.time))
		for i, t in enumerate(self.time):
			if i>0:
				# Time since last check
				DeltaTau = t- self.time[i-1]
				BurnRate = self.Uburn[i-1]
				# deltat = 1/BurnRate
				print(DeltaTau*BurnRate)

				times = np.linspace(0, DeltaTau, int(DeltaTau*BurnRate))
				CNF = 0
				for j in range(len(times)):
					CNF += self.CDF(DeltaTau-times[j])
				CummulitiveneutrinoFlux[i] = CummulitiveneutrinoFlux[i-1]+CNF

		
		self.neutrinoFlux = dp.derivative(CummulitiveneutrinoFlux, self.time)
		return self


	def convolute(self):
		'''
		Move CDF to each point in time on the burnt function
		and calculate the neutrinos emitted at that time
		'''
		self.neutrinoFlux = 0
		flux = {'timeframe':[], 'CDF':[], 'neuFlux':[], 'neutrinoFlux':[]}
		flux = pd.DataFrame(flux, columns= ['timeframe','CDF', 'neuFlux','neutrinoFlux'])
		for i in range(len(self.time)):
			#CDF is 0 before time = 0, shifting t0
			#define time[i] = t0 = 0, therefore CDF = 0 for t<t0.
			# Take rest of time t[i:]-t[i] to set to offset to 0 
			# and calculate CDF for new time 
			t0 = self.time[i]
			timeFrame = self.time[i:]-t0
			#Calculate CDF
			sud.export_data(timeFrame).contributers().maxContributions(full=True)
			CDFPDF_data = pd.read_csv('./Contributing_chains/'+preferences.simpleTitle+'CDF_PDF_full.csv')
			CDF = np.array(ctl(CDFPDF_data.iloc[0,1]))
			# PDF = ctl(CDFPDF_data.iloc[0,2])
			#Add in extra 0's at start
			if i>0:
				CDFstart = np.zeros(i)
				CDF = np.insert(CDF,0,CDFstart)
			neuFlux = np.array(CDF)*self.burnt[i]
			self.neutrinoFlux = self.neutrinoFlux + neuFlux
			saved_flux = pd.DataFrame([[timeFrame, list(CDF), list(neuFlux), list(self.neutrinoFlux)]], columns= ['timeframe','CDF', 'neuFlux','neutrinoFlux'])
			flux = flux.append(saved_flux)
		flux.to_csv('./Flux data.csv', sep=',')
		return self







	

































# def convolute(x, pdf, reactor_on_period, beginning_buffer):
#     uranium_burn_function = np.zeros(len(x))
#     beginning_buffer_index = find_index(x, beginning_buffer)
#     max_fiss = constants.Reactor_output/constants.U235_fiss
#     convoluted = np.zeros(len(x))
#     for i in range(beginning_buffer_index, len(x)):
#         on_period = x[i]-x[0]
#         if on_period>reactor_on_period:
#             off_time = x[i]-reactor_on_period
#             reactor_off_index = find_index(x, off_time)
#             uranium_burn_function[reactor_off_index:i] = max_fiss
#             fg = np.array(PDF)*np.array(uranium_burn_function)
#             convoluted[i] = np.trapz(fg, x)
#         else:
#             uranium_burn_function[:i] = max_fiss
#             fg = np.array(PDF)*np.array(uranium_burn_function)
#             convoluted[i]= np.trapz(fg, x)
#     plt.plot(x,uranium_burn_function)
#     plt.xscale('log')
#     plt.show()
#     plt.plot(x, convoluted)
#     plt.xscale('log')
#     plt.show()

# convolute(constants.t, PDF, reactor_on, beginning_buffer)







# Space = tOff-tOn
# on_index = find_index(constants.t, tOn)
# off_index = find_index(constants.t, tOff)

# Thermal_power[on_index:off_index] = constants.Reactor_output
# uranium_burn_function = Thermal_power/constants.U235_fiss

# def shift(startingx, x,f):
#     shiftedf = np.zeros(len(x))
#     on_index = find_index(constants.t, startingx)
#     off_index = find_index(constants.t, startingx+Space)
#     print(on_index, off_index)
#     space = constants.t[off_index]-constants.t[on_index]
#     shiftedf[on_index:off_index]=1
#     return shiftedf, space


# #SHow animation
# # fig = plt.figure()
# # ax = plt.axes(xlim=(constants.t[0], constants.t[-1]),ylim=(-0.5,1.5))
# # line, = ax.plot([],[])
# # def init():
# #     line.set_data([],[])
# #     return line,

# # xdata, ydata = [],[]

# # def animate(i):
# #     x = constants.t
# #     y,_ = shift(x[i],x,uranium_burn_function)
# #     line.set_data(x,y)
# #     return line,

# # anim = FuncAnimation(fig, animate, init_func=init, frames = len(constants.t), interval = 20, blit=True)
# # plt.xscale('log')
# # plt.show()

# # spaceing = np.zeros(len(constants.t))
# # for i,_ in enumerate(constants.t):
# #     shiftedf, space = shift(constants.t[i], constants.t, uranium_burn_function)
# #     spaceing[i]=space
# # plt.plot(constants.t, spaceing)
# # plt.xscale('log')
# # plt.show()

# fig = plt.figure()

# def convolute(x, f, g):
#     L = len(f)
#     convoluted = np.zeros(L)
#     for i in range(0, L):
#         shifted,_ = list(shift(x[i],x,f))
#         # plt.plot(constants.t, shifted)
#         # plt.xscale('log')
#         # plt.plot(constants.t, g)
#         # plt.show()
#         fg = [shifted[m]*g[m] for m in range(L)]
#         plt.plot(constants.t, fg)
#         plt.xscale('log')
#         plt.show()
#         convoluted[i] = np.trapz(constants.t, fg)
#     return convoluted

# convoluted = convolute(constants.t, uranium_burn_function, PDF)
# plt.plot(constants.t, convoluted)
# plt.xscale('log')
# plt.show()



# # x = np.logspace(0,2,50)
# # f = [x[-i]/max(x) for i in range(len(x))]
# # g = np.zeros(len(x))
# # start = 10**(0.75)
# # end = 10**(1.5)
# # on_index = find_index(x, start)
# # off_index = find_index(x, end)
# # g[on_index:off_index] = 1
# # plt.plot(x,f)
# # plt.plot(x,g)
# # plt.xscale('log')
# # plt.show()
# # convoluted = convolute(x, g, f)












# # x = np.linspace(0,10,100)
# # def neutrinoFlux(x, f, g):
# #     #Convolute uranium burn function with PDF
# #     #g*f -- inverse f and shift accross g
# #     L = len(f)
# #     convoluted = np.zeros(L)
# #     #shifting inverse_Uburn
# #     for i in range(1, (L+1)*2):
# #         if i<=L:
# #             shift = list(f[-i:])
# #             extraPart = list(np.zeros(L-i))
# #             full = shift+extraPart
# #         elif L*2>=i>L:
# #             shift = list(f[:(L-i)])
# #             extraPart = list(np.zeros(i-L))
# #             full = extraPart+shift
# #         if i%2 == 0:
# #             j = int(i/2)
# #             # fg = [full[m]*g[m] for m in range(L)]
# #             convoluted[j-1] = np.trapz(full, g)
# #             plt.plot(x, full)
# #             plt.show()
# #     return convoluted

	
# # f = np.zeros(len(x))
# # f[25:75] = 1
# # g = np.zeros(len(x))
# # g[25:75] = x[0:50]/max(x[0:50])
# # g = [g[-i] for i in range(len(g))]
# # convolute = neutrinoFlux(constants.t, uranium_burn_function, PDF)
# # fig, ((ax1,ax3)) = plt.subplots(2,1, figsize=(5,7))

# # colour = 'tab:red'
# # ax1.set_ylabel('Uranium burn rate ($U^{235}$/s)')
# # ax1.set_title('$U^{235}$ function and PDF')
# # ax1.plot(constants.t, uranium_burn_function, color= 'r')
# # ax1.tick_params(axis='y', labelcolor=colour)
# # ax1.set_xscale('log')
# # ax1.grid()
# # ax2 = ax1.twinx()

# # colour = 'tab:blue'
# # ax2.set_ylabel('Neutrino PDF')
# # ax2.plot(constants.t,PDF, color='b')
# # ax2.set_yscale('log')
# # # ax2.('log')
# # ax2.tick_params(axis='y', labelcolor=colour)
# # yLab = r"$\phi_{\bar{\nu}} (t)$"
# # ax3.set_xlabel('time (s)')
# # ax3.set_ylabel(yLab)
# # ax3.plot(constants.t, convolute)
# # ax3.set_xscale('log')
# # title = r"$\bar{\nu}$"
# # ax3.set_title('Convoluted '+title+' flux')
# # fig.subplots_adjust(hspace = 17)
# # txt = 'Assuming $U^{235}$ is providing 100% thermal power. Thermal power converted to fiss/s using average $U^{235}$ energy/fiss of $200x10^6$ eV. Max thermal power begins at t=$10^{'+str(int(np.log10(tOn)))+'}$ s, with nominal power assumed to be $20x10^6$ J/s, and ending at t=$10^{'+str(int(np.log10(tOff)))+'}$ s.'
# # plt.figtext(0.5,0.435, txt, wrap=True, horizontalalignment='center', fontsize =8)
# # fig.tight_layout()
# # ax3.grid()
# # plt.show()
