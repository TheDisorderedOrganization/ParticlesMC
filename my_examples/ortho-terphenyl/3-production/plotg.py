import atooms.postprocessing as pp
import numpy as np
import pandas as pd
from atooms.trajectory import TrajectoryXYZ
from matplotlib import pyplot as plt

path = 'chains/1/trajectory.xyz'
th = TrajectoryXYZ(path)

gr_all = pp.RadialDistributionFunction(th, dr=0.03, norigins=len(th))
gr_all.compute()

plt.plot(gr_all.grid, gr_all.value)
plt.xlabel("$r$")
plt.ylabel("$g(r)$")
plt.savefig('g(r).png')