import nubaseSorter as ns
import naming
import findChain
import numpy as np
import preferences
import pandas as pd
import matplotlib.pyplot as plt
import constants as cst
import time as ti
from ast import literal_eval

with open('./databases/'+preferences.simpleTitle+'CDF_full.txt', 'r') as f:
    lines = f.readline()


#Just do every 1 second
time = np.linspace(0, 6e6, 300000)
time_step = 20
CDF = eval('lambda t:' + lines)
time_start = ti.time()
count = len(open('./databases/'+str(time_step)+'s_CDF_high_density.txt').readlines(  ))
print(count)
for i,t in enumerate(time):
    if i > count:
        Output = CDF(t)
        current_time = ti.time()
        diff_t = current_time - time_start
        Average_rate = diff_t / (i + 1)
        extrap_time = Average_rate*(len(time)-(i+1))
        print(extrap_time/3600, 'Hours', end='\r')
        with open('./databases/'+str(time_step)+'s_CDF_high_density.txt', 'a') as f:
            f.write(str(Output)+'\n')
