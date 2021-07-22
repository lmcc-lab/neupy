import numpy as np
import preferences

#Graph time conversions

if preferences.logspace == False:
    graph_start = preferences.graph_start*3600                      #s
    graph_end = preferences.graph_end*3600                          #s
    offset_t = preferences.offset_t*3600                            #s
    t = np.linspace(graph_start,graph_end,preferences.resolution)
elif preferences.offset_t == 0:
    t = np.logspace(preferences.graph_start, preferences.graph_end, preferences.resolution)
else:
    t = np.logspace(preferences.graph_start+preferences.offset_t, preferences.graph_end, preferences.resolution)

#Constants
sigma_U_thermal = 584.3                             #Microscopic cross section of U 235 for thermal neutrons Barns
sigma_U_fast = 1                                    #Microscopic cross section of U 235 for fast neutrons Barns
barns_cm = 1e-24                                    #cm^2/barn
sigma_U = sigma_U_fast*barns_cm                     #Microscopic cross section of U235 in cm^2
amu = 235.0439299                                   #g/mol
amu_g = 1.6605402e-24                               #g/mol to g
U235_mass_g = amu*amu_g                             #g
Na = 6.0221409e23                                   #Avogadros number
U235_den = 19.1                                     #g/cm^3
atom_den = U235_den*Na/amu                          #cm^-3
Sum_f = atom_den*sigma_U                            #Macroscopic cross section of U 235 cm^-1
U235_fiss = 200e6                                   #Average U235 output energy eV
eV_joul = 1.60217662e-19                            #J/eV
U235_fiss = U235_fiss*eV_joul                       #J

#Reactor dependent constants

N_op = 1                                            #Neutron flux when operating (as a percentage of normal operating flux)
Reactor_output = 20e6                               #Watts (J/s)
V_reactor = 745402.44                               #cm^3

#Constants again

fiss_psec = Reactor_output/U235_fiss                #Fission events of U235/s
N_psec = fiss_psec*1.3                              #Assuming average of 1.3 neutrons emmited s^-1
phi_0 = Reactor_output/(Sum_f*V_reactor*U235_fiss)  #Neutron flux at nominal reactor operation cm^-2s^-1
