import atooms.postprocessing as pp
import numpy as np
import pandas as pd
from atooms.trajectory import TrajectoryXYZ
from matplotlib import pyplot as pl

path = 'trajectory.xyz'

th = TrajectoryXYZ(path)

#print(th[-1])  ## last frame
#print(len(th)) ## number of frame

##### test tutorial #####

# # We let postprocessing choose a reasonable number of time origins for the average
# gr = pp.RadialDistributionFunction(th, dr=0.03)
# gr.compute()
# # We average over all the frames
# gr_all = pp.RadialDistributionFunction(th, dr=0.03, norigins=len(th))
# gr_all.compute()

# pl.plot(gr.grid, gr.value, label='Default')
# pl.plot(gr_all.grid, gr_all.value, label='All time origins')
# pl.legend()
# pl.xlabel("$r$")
# pl.ylabel("$g(r)$")
# pl.show()

##### MSD #####

msd = pp.MeanSquareDisplacement(th, rmax=3.0)
msd.compute()
pl.plot(msd.grid, msd.value, '-o')
pl.xlabel("$t$")
pl.ylabel("MSD")

pl.show()