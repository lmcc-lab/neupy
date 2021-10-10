import numpy as np
from dependencies import derivative

def correct_emissions(state, size):
    with open('./experimental_data/'+state+'_experimental_neutrinos', 'r') as f:
        lines = f.readlines()
    lines = [float(line.split('\n',1)[0]) for line in lines]

    time = np.linspace(20, len(lines)*20, len(lines))
    emissionRate = derivative(lines[-size:], time[-size:])
    return emissionRate