from sympy import *
import pandas as pd
import matplotlib.pyplot as plt
import constants as ct
import time
import numpy as np

#Import experimental_data
file_name = './experimental_data/ANSTO Reactor Shutdown & Startup.xlsx'
sheet = 'Shutdown Stats'
data = pd.read_excel(io = file_name, sheet_name = sheet)

# Convert the thermal power to burn rate
zeta = 0.94*data['Reactor Power']/ct.U235_fiss
T = data['Date']-data['Date'][0]
T = [v.total_seconds() for v in T]

# print(T)

# Set up symbols
tau = Symbol("tau")

# Make zeta continuous
arg = ''
for i, val in enumerate(T):
    if i < len(T)-1:
        arg = arg + '('+str(zeta[i])+', And(tau>='+str(T[i])+', tau<'+str(T[i+1])+')), '
    else:
        arg = arg + '('+str(zeta[i])+', tau >='+str(T[i])+')'

zeta_cont = eval('Piecewise('+arg+')')
