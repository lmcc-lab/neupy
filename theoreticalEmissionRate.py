import matplotlib.pyplot as plt
import numpy as np
import neutrinoEmissions as nf
import random

endtime = 9         #hrs
endtime = endtime * 3600    #s
stepsize = 20               #s

time = np.linspace(0, 5000, int(5000/stepsize) + 1)

Thermal = np.zeros(len(time))

for i in range(50, len(Thermal)):
    Thermal[i] = 19.93592183+random.randrange(-5, 5)/10

# for i in range(50, 60):
#     Thermal[i] = 15+random.randrange(-25, 25)/10
#
# for i in range(60, 70):
#     Thermal[i] = 5+random.randrange(-25, 25)/10

# Convert Thermal MW to W
Thermal = Thermal*1e6

data = nf.neuFlux(time, Thermal, state='on').flux()
uburn = data.Uburn

epsilon = data.neutrinoEmissionRate
epsilon[0] = epsilon[1]

cum = data.cummulitiveNeutrinos

fig, ((ax1, ax2)) = plt.subplots(2,1, figsize=(5,7))
ax1.plot(time, uburn, color='r')
ax1.grid()
ax1.set_ylabel(r'$\zeta(t)$')
ax1.set_title('Reactor burn rate')

ax2.plot(time, epsilon)
ax2.grid()
ax2.set_ylabel(r'$\varepsilon(t)$')
ax2.set_title('Neutrino emission rate')
ax2.set_xlabel('time (s)')
plt.show()