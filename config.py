from neupy.databases.load_databases import load_config
from neupy import path

ANSTO_config = load_config(f'{path}\\ANSTO_reactor', 'ansto_specifications.yaml')
ANSTO_thermal_neutron_flux = float(ANSTO_config['general_description_of_the_opal_facility']['max_neutron_flux']['thermal_irradiation_position_in_reflector_cm-2s-1'])

max_chain_depth = 150
u235_energy_per_fission_MeV = 3.20435313e-11

# https://wwwndc.jaea.go.jp/cgi-bin/Tab80WWW.cgi?lib=J40&iso=U235
u233_thermal_cross_section_barns = 531.3
u234_thermal_cross_section_barns = 67.02e-3
u235_thermal_cross_section_barns = 585.1
u236_thermal_cross_section_barns = 259.4e-6
u238_thermal_cross_section_barns = 16.80e-6

u233_average_cross_section_barns = 1.908
u234_average_cross_section_barns = 1.222
u235_average_cross_section_barns = 1.218
u236_average_cross_section_barns = 579.2e-3
u238_average_cross_section_barns = 306.4e-3

xe_neutron_cross_section_cm2 = 2.805e-18
amu_to_g = 1.660538921e-24
proton_amu = 1
neutron_amu = 1
kev_to_amu = 1.07354466425775e-6
barnes_to_cm2 = 1e-24

fissile_fuel_cross_sections_barns = {"U233": {0.0253: u233_thermal_cross_section_barns,
                                            "spectrum_average": u233_average_cross_section_barns},
                                   "U234": {0.0253: u234_thermal_cross_section_barns,
                                            "spectrum_average": u234_average_cross_section_barns},
                                   "U235": {0.0253: u235_thermal_cross_section_barns,
                                            "spectrum_average": u235_average_cross_section_barns},
                                   "U236": {0.0253: u236_thermal_cross_section_barns,
                                            "spectrum_average": u236_average_cross_section_barns},
                                   "U238": {0.0253: u238_thermal_cross_section_barns,
                                            "spectrum_average": u238_average_cross_section_barns}}