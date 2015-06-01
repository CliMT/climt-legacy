#!/usr/bin/env python

from numpy import *
from climt import *
from pylab import *
from myCmap import *

#Parameters
kwargs = {}
kwargs['afc'] = 0.4
kwargs['dt'] = 600.

nhours = 1000
nsteps = int((nhours*3600)/kwargs['dt'])
print 'Num Steps: ', nsteps

dycore = dynamics(scheme='gfs', **kwargs)

oldU = 0
newU = 0
incU = 0

for i in range(nsteps):
    oldU = dycore['U'].copy()

    dycore.step()

    incU = dycore.Inc['U'].copy()
    newU = dycore['U'].copy()

contourf(average(oldU - newU + incU, axis=0))
colorbar()
figure()
zonalmean = average(newU, axis=0).transpose()
maxLev = amax(abs(zonalmean))
levs = linspace(-maxLev, maxLev, 16)
contourf(zonalmean ,levels = levs, cmap=joyDivCmapRdBl, origin='lower')
colorbar()
show()

dycore.Extension.shutDownModel()
