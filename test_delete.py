import ast

import numpy as np
import matplotlib.pyplot as plt
test = [1,2,3,4]

# function = map(lambda t, dc: np.exp(-dc*t)/np.product([test[m]-dc for m in range(len(test)) if test[m] != dc]))

beforeCDF = list(np.zeros(1).tolist())
print(beforeCDF)
list2 = [1,2,3,4,6]
print(list2)
beforeCDF.extend(list2)
print(beforeCDF)


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