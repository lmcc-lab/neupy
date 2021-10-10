import pandas as pd
from pandas import ExcelWriter
from pandas import ExcelFile
import matplotlib.pyplot as plt
import numpy as np
import constants
import preferences
from dependencies import conv_str_to_list as ctl
from numpy import exp as exp
import neutrinoEmissions as nf
from dependencies import derivative

shutdown_df = pd.read_excel('./experimental_data/ANSTO Reactor Shutdown & Startup.xlsx', sheet_name='Shutdown Stats')
startup_df = pd.read_excel('./experimental_data/ANSTO Reactor Shutdown & Startup.xlsx', sheet_name='Startup')

shutdown_time = pd.to_timedelta(shutdown_df['Date']-shutdown_df['Date'][0]).astype('timedelta64[s]').astype(int)
startup_time = pd.to_timedelta(startup_df['Date']-startup_df['Date'][0]).astype('timedelta64[s]').astype(int)

time_end = shutdown_df['Date'][0]
time_start = startup_df['Date'][len(startup_time)-1]


Thermal_on = np.array(startup_df['Reactor Power'])
t_on = np.array(startup_time)
Thermal_off = np.array(shutdown_df['Reactor Power'])
t_off = np.array(shutdown_time)


#Convert MW to W
Thermal_on = Thermal_on * 1e6
Thermal_off = Thermal_off * 1e6
Thermal_off = np.append(Thermal_off, np.zeros(int(6*24*3600/20)))
t_off = np.append(t_off, np.linspace(t_off[-1], t_off[-1]+6*24*3600, int(7*24*3600/20)))

print('starting')
import test_delete
data_off = nf.neuFlux(t_off, Thermal_off, state = 'off').emissions()
data_on = nf.neuFlux(t_on, Thermal_on, state = 'on').emissions()
uburn_off = data_off.Uburn
uburn_on = data_on.Uburn
epsilon_on = data_on.neutrinoEmissionRate
epsilon_off = data_off.neutrinoEmissionRate
print(epsilon_off[-1])

# epsilon = data.neutrinoEmissionRate
# cumalitive = data.cummulitiveNeutrinos
# epsilon[0] = epsilon[1]
# epsilon = derivative(cumalitive, t)
# uburn = data.Uburn

fig, ((ax1,ax3),(ax2, ax4)) = plt.subplots(2,2, figsize=(7,7))
colour = 'tab:red'
ax1.set_ylabel(r'$\zeta (t)$ $(U^{235}/s)$')
ax1.set_title('Turn on')
ax1.plot(t_on, uburn_on, color= 'r')
ax1.tick_params(axis='y', labelcolor=colour)
# ax1.set_xscale('log')
ax1.grid()

colour = 'tab:blue'
yLab = r"$\varepsilon (t)$"
# yLab = "Cumulative neutrinos"
ax2.set_xlabel('time (s)')
ax2.set_ylabel(yLab)
ax2.plot(t_on, epsilon_on)
# ax3.set_xscale('log')
ax2.grid()

colour = 'tab:red'
ax3.set_ylabel(r'$\zeta (t)$ $(U^{235}/s)$')
ax3.set_title('Turn off')
ax3.plot(t_off, uburn_off, color= 'r')
ax3.tick_params(axis='y', labelcolor=colour)
# ax1.set_xscale('log')
ax3.grid()

colour = 'tab:blue'
yLab = r"$\varepsilon (t)$"
# yLab = "Cumulative neutrinos"
ax4.set_xlabel('time (s)')
ax4.set_ylabel(yLab)
ax4.plot(t_off, epsilon_off)
# ax3.set_xscale('log')
ax4.grid()

fig.subplots_adjust(hspace = 17)
# txt = 'Assuming $U^{235}$ is providing 100% thermal power. Thermal power converted to fiss/s using average $U^{235}$ energy/fiss of $200x10^6$ eV. Max thermal power begins at t=$10^{'+str(int(np.log10(tOn)))+'}$ s, with nominal power assumed to be $20x10^6$ J/s, and ending at t=$10^{'+str(int(np.log10(tOff)))+'}$ s.'
# plt.figtext(0.5,0.435, txt, wrap=True, horizontalalignment='center', fontsize =8)
fig.tight_layout()



plt.show()

# for val in cumalitive:
#     with open('./experimental_data/experimental_neutrinos', 'a') as f:
#         f.write(str(val) + '\n')