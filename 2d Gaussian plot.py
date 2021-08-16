import numpy as np

import matplotlib.pyplot as plt

x = np.linspace(-7, 7, 100)
y = np.linspace(-7, 7, 100)
X, Y = np.meshgrid(x,y)
sigma_x = 1
sigma_y = 1

def Gaussian(sigma, x_axis):
    return 1/(np.sqrt(2*np.pi)*sigma) * np.exp(-1/2 * x_axis**2/sigma**2)

gaussian = Gaussian(sigma_x, X)*Gaussian(sigma_y, Y)

fig = plt.figure()
ax = plt.axes(projection='3d')
ax.contour3D(X, Y, gaussian, 50, cmap='binary')
ax.set_xlabel('x')
ax.set_ylabel('t')
ax.set_zlabel('$\psi$')
ax.set_zticks([])
ax.set_xticks([])
ax.set_yticks([])
plt.show()