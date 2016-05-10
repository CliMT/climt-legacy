#!/usr/bin/env python

from numpy import *
from dynamics import *
from pylab import *
from myCmap import *

#Parameters
kwargs = {}
kwargs['afc'] = 0.4
kwargs['dt'] = 600.

nhours = 100
nsteps = int((nhours*3600)/kwargs['dt'])
print 'Num Steps: ', nsteps

dycore = dynamics(scheme='gfs', **kwargs)

print '==========================='
oldU = 0
newU = 0
incU = 0


for i in range(nsteps):

    
    print '==========================='
    dycore.step()



dycore.Extension.shutDownModel()

