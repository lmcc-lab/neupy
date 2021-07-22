import pandas as pd
from pandas import ExcelWriter
from pandas import ExcelFile
import matplotlib.pyplot as plt
import numpy as np
import constants
import preferences
from dependencies import conv_str_to_list as ctl
shutdown_df = pd.read_excel('./experimental data/ANSTO Reactor Shutdown & Startup.xlsx', sheet_name='Shutdown Stats')
startup_df = pd.read_excel('./experimental data/ANSTO Reactor Shutdown & Startup.xlsx', sheet_name='Startup')

shutdown_time = pd.to_timedelta(shutdown_df['Date']-shutdown_df['Date'][0]).astype('timedelta64[s]').astype(int)
startup_time = pd.to_timedelta(startup_df['Date']-startup_df['Date'][0]).astype('timedelta64[s]').astype(int)

constants.t = np.array(shutdown_time)
preferences.full = False
import neutpy
import neutrinoFlux as nf
# neutpy.single_fission_event(constants.t, update=True)
#Generating PDF

uranium_burn_function = np.array(shutdown_df['Reactor Power'])/constants.U235_fiss
CDFPDF_data = pd.read_csv('CDF_PDF_full.csv')
CDF = ctl(CDFPDF_data.iloc[0,1])
PDF = ctl(CDFPDF_data.iloc[0,2])
convolute = nf.neutrinoFlux(constants.t, uranium_burn_function, PDF)

fig, ((ax1,ax3)) = plt.subplots(2,1, figsize=(5,7))
colour = 'tab:red'
ax1.set_ylabel('Uranium burn rate ($U^{235}$/s)')
ax1.set_title('$U^{235}$ function and PDF')
ax1.plot(constants.t, uranium_burn_function, color= 'r')
ax1.tick_params(axis='y', labelcolor=colour)
ax1.set_xscale('log')
ax1.grid()
ax2 = ax1.twinx()

colour = 'tab:blue'
ax2.set_ylabel('Neutrino PDF')
ax2.plot(constants.t,PDF, color='b')
ax2.set_yscale('log')
# ax2.('log')
ax2.tick_params(axis='y', labelcolor=colour)
yLab = r"$\phi_{\bar{\nu}} (t)$"
ax3.set_xlabel('time (s)')
ax3.set_ylabel(yLab)
ax3.plot(constants.t, convolute)
ax3.set_xscale('log')
title = r"$\bar{\nu}$"
ax3.set_title('Convoluted '+title+' flux')
fig.subplots_adjust(hspace = 17)
# txt = 'Assuming $U^{235}$ is providing 100% thermal power. Thermal power converted to fiss/s using average $U^{235}$ energy/fiss of $200x10^6$ eV. Max thermal power begins at t=$10^{'+str(int(np.log10(tOn)))+'}$ s, with nominal power assumed to be $20x10^6$ J/s, and ending at t=$10^{'+str(int(np.log10(tOff)))+'}$ s.'
# plt.figtext(0.5,0.435, txt, wrap=True, horizontalalignment='center', fontsize =8)
fig.tight_layout()
ax3.grid()
plt.show()