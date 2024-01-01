from neupy import Neupy
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


neupy = Neupy()
# df = pd.read_excel('ANSTO_reactor/ANSTO Reactor Shutdown & Startup.xlsx')
# old_results = pd.read_csv('ANSTO_reactor/off_experimental_neutrinos')
df = pd.read_excel('ANSTO_reactor/ANSTO reactor on.xlsx')
old_results = pd.read_csv('ANSTO_reactor/on_experimental_neutrinos')

start_date = df.loc[0, 'Date']
end_date = df.loc[df.shape[0]-1, 'Date']
delta_time = (end_date - start_date).seconds  # seconds

thermal_power = df['Reactor Power'].values * 1e6


time = np.linspace(0, delta_time, len(thermal_power))

burn_rate = neupy.convert_thermal_to_burn_rate(thermal_power)
cum_neutrinos, emission_rate = neupy.neutrino_emission_rate(burn_rate, time, load_saved=True, assumed_on_prior_s=0)

fig, ax = plt.subplots()
thermal_line, =ax.plot(time, thermal_power, color='r', label='thermal power')
# old_neutrino_line, = ax.plot(time, old_results['neutrinos'].values, label='Old results', color='r')

axtwin = ax.twinx()
for element, ne_dict in emission_rate.items():
    for ne, neutrino_emission in ne_dict.items():
        new_neutrino_line, =axtwin.plot(time, neutrino_emission, label='{0}, Neutron energy {1:g}'.format(element, ne))

old_slope = neupy.slope(old_results['neutrinos'].values, time[1]-time[0])
old_neutrino_line, = axtwin.plot(time, old_slope, color='g')

ax.set_xlabel('time (s)')
# ax.set_ylabel('Reactor thermal power (MW)')
ax.set_ylabel(r'Old cumulative neutrinos ($\nu$)')
# axtwin.set_ylabel(r'neutrino emission rate ($\nu/s$)')
axtwin.set_ylabel(r'New cumulative neutrinos ($\nu$)')
ax.legend([thermal_line, new_neutrino_line, old_neutrino_line], ['Thermal power', 'New neutrino results', 'Old neutrino results'])

plt.show()
