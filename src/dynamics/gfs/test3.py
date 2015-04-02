from _gfs_dynamics import _gfs_dynamics
import numpy as np;


test = _gfs_dynamics(192,94);

test.initModel();

#Initial conditions

nug,nvg,nvtg,ntrg,npsg,npg = test.getResult();

print nug.shape

for i in range(10):
    test.oneStepForward();

test.shutDownModel();
