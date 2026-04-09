import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt('chains/1/energy.dat')

energy = data[:,1]
time = data[:,0]

plt.plot(time,energy)
plt.xlabel('t, MC step')
plt.ylabel('energy')
plt.savefig('energy.png')