import pandas as pd
from pandas import ExcelWriter
from pandas import ExcelFile
import matplotlib.pyplot as plt
import numpy as np
import constants
import preferences
from dependencies import conv_str_to_list as ctl
from numpy import exp as exp
import neutrinoFlux as nf

shutdown_df = pd.read_excel('./experimental_data/ANSTO Reactor Shutdown & Startup.xlsx', sheet_name='Shutdown Stats')
startup_df = pd.read_excel('./experimental_data/ANSTO Reactor Shutdown & Startup.xlsx', sheet_name='Startup')

shutdown_time = pd.to_timedelta(shutdown_df['Date']-shutdown_df['Date'][0]).astype('timedelta64[s]').astype(int)
startup_time = pd.to_timedelta(startup_df['Date']-startup_df['Date'][0]).astype('timedelta64[s]').astype(int)

t = np.array(shutdown_time)

uranium_burn_function = np.array(shutdown_df['Reactor Power'])/constants.U235_fiss
with open('./Contributing_chains/NotSimpleCDF_full.txt', 'r') as f:
    CDF = f.readline()


epsilon = nf.neuFlux(t, uranium_burn_function, CDF, transient = 'off').flux().neutrinoEmissionRate

fig, ((ax1,ax3)) = plt.subplots(2,1, figsize=(5,7))
colour = 'tab:red'
ax1.set_ylabel('Uranium burn rate ($U^{235}$/s)')
ax1.set_title('$U^{235}$ function and PDF')
ax1.plot(t, uranium_burn_function, color= 'r')
ax1.tick_params(axis='y', labelcolor=colour)
# ax1.set_xscale('log')
ax1.grid()

colour = 'tab:blue'
yLab = r"$\phi_{\bar{\nu}} (t)$"
ax3.set_xlabel('time (s)')
ax3.set_ylabel(yLab)
ax3.plot(t, epsilon[1:])
# ax3.set_xscale('log')
title = r"$\bar{\nu}$"
ax3.set_title('Convoluted '+title+' Emission rate')
fig.subplots_adjust(hspace = 17)
# txt = 'Assuming $U^{235}$ is providing 100% thermal power. Thermal power converted to fiss/s using average $U^{235}$ energy/fiss of $200x10^6$ eV. Max thermal power begins at t=$10^{'+str(int(np.log10(tOn)))+'}$ s, with nominal power assumed to be $20x10^6$ J/s, and ending at t=$10^{'+str(int(np.log10(tOff)))+'}$ s.'
# plt.figtext(0.5,0.435, txt, wrap=True, horizontalalignment='center', fontsize =8)
fig.tight_layout()
ax3.grid()
plt.show()