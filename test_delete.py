import numpy as np
import matplotlib.pyplot as plt
from dependencies import derivative
from functools import reduce

# def _channel_to_hex(color_val: int):
#     raw: str = hex(color_val)[2:]
#     return raw.zfill(2)


# def rgb_to_hex(red: int, green: int, blue: int):
#     return "#" + _channel_to_hex(red) + _channel_to_hex(green) + _channel_to_hex(blue)

# print(rgb_to_hex(1,1,1))

# x = np.linspace(0,15*np.pi,100)
# y = np.sin(x)
# plt.plot(x,y*x)
# plt.show()

A=6
Z=3
A1 = A/2+A%2
Z1 = Z/2 +Z%2
print(A1, Z1)

'''
N0 = init[0]*np.exp(-dc[0]*t)
N1 = init[0]*dc[0]*(np.exp(-dc[0]*t)/(dc[1]-dc[0])+np.exp(-dc[1]*t)/(dc[0]-dc[1])) + init[1]*np.exp(-dc[1]*t)
N2 = init[0]*dc[1]*dc[0]*(np.exp(-dc[0]*t)/((dc[1]-dc[0])*(dc[2]-dc[0]))+np.exp(-dc[1]*t)/((dc[0]-dc[1])*(dc[2]-dc[1]))+np.exp(-dc[2]*t)/((dc[0]-dc[2])*(dc[1]-dc[2]))) + init[1]*dc[1]*(np.exp(-dc[1]*t)/(dc[2]-dc[1])+np.exp(-dc[2]*t)/(dc[1]-dc[2])) + init[2]*np.exp(-dc[2]*t)

dN0 = init[0]*-dc[0]*np.exp(-dc[0]*t)
dN1 = init[0]*dc[0]*(-dc[0]*np.exp(-dc[0]*t)/(dc[1]-dc[0])-dc[1]*np.exp(-dc[1]*t)/(dc[0]-dc[1])) + init[1]*-dc[1]*np.exp(-dc[1]*t)
dN2 = init[0]*dc[1]*dc[0]*(-dc[0]*np.exp(-dc[0]*t)/((dc[1]-dc[0])*(dc[2]-dc[0]))-dc[1]*np.exp(-dc[1]*t)/((dc[0]-dc[1])*(dc[2]-dc[1]))-dc[2]*np.exp(-dc[2]*t)/((dc[0]-dc[2])*(dc[1]-dc[2]))) + init[1]*dc[1]*(-dc[1]*np.exp(-dc[1]*t)/(dc[2]-dc[1])-dc[2]*np.exp(-dc[2]*t)/(dc[1]-dc[2])) -dc[2]*init[2]*np.exp(-dc[2]*t)
CDF = N0*0+N1*1+N2*2
PDF = dN0*0+dN1*1+dN2*2
plt.plot(t, N0, label='N0')
plt.plot(t,dN0, linestyle = ':',label='dN0')
plt.plot(t,N1,label='N1')
plt.plot(t,dN1, linestyle = ':',label='dN1')
plt.plot(t,N2,label='N2')
plt.plot(t,dN2, linestyle = ':',label='dN2')
plt.legend()
plt.xscale('log')
plt.show()
plt.plot(t, CDF/max(CDF), label = 'CDF')
plt.plot(t, PDF/max(PDF), label = 'PDF')
plt.xscale('log')
plt.legend()
plt.show()
'''