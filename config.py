from neupy.databases.load_databases import load_config
from neupy import path

ANSTO_config = load_config(f'{path}\\ANSTO_reactor', 'ansto_specifications.yaml')
ANSTO_thermal_neutron_flux = float(ANSTO_config['general_description_of_the_opal_facility']['max_neutron_flux']['thermal_irradiation_position_in_reflector_cm-2s-1'])

max_chain_depth = 150
u235_energy_per_fission_MeV = 3.20435313e-11
xe_neutron_cross_section = 2.805e-18
