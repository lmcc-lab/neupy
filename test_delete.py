import ast

import numpy as np
from numpy import exp as exp
import matplotlib.pyplot as plt
from dependencies import derivative
test = [1,2,3,4]

# function = map(lambda t, dc: np.exp(-dc*t)/np.product([test[m]-dc for m in range(len(test)) if test[m] != dc]))
def correct_emissions(state):
    with open('./experimental_data/'+state+'_experimental_neutrnios', 'r') as f:
        lines = f.readlines()
    lines = [float(line.split('\n',1)[0]) for line in lines]

    time = np.linspace(20, len(lines)*20, len(lines))
    emissionRate = derivative(lines[-264:], time[-264:])
    return emissionRate

def expSum(dcList):
    # Generate exp code
    # exp template
    print('Input dcList',dcList)
    denom = np.product([dcList[m]-dcList[0] for m in range(len(dcList)) if m != 0])
    print('Denominator',denom)
    template = 'np.exp(-'+str(dcList[0])+'*t)/('+str(denom)+')'
    for i in range(len(dcList)-1):
        denom = np.product([dcList[m]-dcList[i+1] for m in range(len(dcList)) if m != i+1])
        template = template.join(['','+np.exp(-'+str(dcList[i+1])+'*t)/('+str(denom)+')'])
    starting_prod = np.product([dcList[i] for i in range(len(dcList)-1)])
    # template = ''.join([str(starting_prod)+'*(', template, ')'])
    template = str(starting_prod) + '*(' + template + ')'
    return template

def neutrino_emission_CDF(weighting, nVector, decayConstants):
    CDF_contrib = str(nVector[0])+'*'+expSum(decayConstants[:1])
    for i, dc in enumerate(decayConstants):
        if i>0:
            CDF_contrib = CDF_contrib +'+'+  str(nVector[i])+'*'+expSum(decayConstants[:i+1])
    CDF_contrib = str(weighting)+'*(' + CDF_contrib + ')'

    return CDF_contrib