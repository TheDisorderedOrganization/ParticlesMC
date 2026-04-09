import numpy as np
import matplotlib.pyplot as plt

##### Energy #####

N = 3*(4**3)

data     = np.loadtxt('energy.dat')
time     = data[:,0]
energies = data[:,1]/N

plt.plot(time,energies)
plt.xlabel('MC steps')
plt.ylabel('energy per particle')
plt.savefig('energy.png')