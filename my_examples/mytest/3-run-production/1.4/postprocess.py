import atooms.postprocessing as pp
import numpy as np
import pandas as pd
from atooms.trajectory import TrajectoryXYZ
from matplotlib import pyplot as pl

path = 'chains/1/trajectory.xyz'
pathE = 'chains/1/energy.dat'

th = TrajectoryXYZ(path)
data = np.loadtxt(pathE)
E = data[:,1]/192
t = data[:,0]

th = TrajectoryXYZ(path)

##### energy #####

pl.figure()
pl.plot(t,E)
pl.xlabel('t')
pl.ylabel('energy per particle')
pl.savefig('energy.png')
pl.close()

##### MSD #####

msd = pp.MeanSquareDisplacement(th)#, rmax=3.0)
msd.compute()
pl.figure()
pl.plot(msd.grid, msd.value, '-o')  # on enlève t=0
pl.xlabel("$t$")
pl.ylabel("MSD")
pl.savefig('msd.png')
pl.close()
