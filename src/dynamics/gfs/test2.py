from _gfs_dynamics import _gfs_dynamics
import numpy as np;
import scipy as sp;


test = _gfs_dynamics(192,94);

test.initModel();

#Initial conditions
for i in range(5):
    test.oneStepForward();

nug,nvg,nvtg,ntrg,npsg,npg = test.getResult();

nugt = nug + sp.random.normal(scale=.0000001,size=nug.shape);

test.oneStepForward();

#Result of one forward integration
nug1,nvg1,nvtg1,ntrg1,npsg1,npg1 = test.getResult();

#Result of sending ICs in and stepping forward
#These are only increments (fnew - fold)
nug2,nvg2,nvtg2,ntrg2,npsg2,npg2 = \
        test.driver(nugt,nvg,nvtg,ntrg,npsg,npg,simTime=0);

#Should be very close to zero
print 'udiff: ', np.amax(abs(nug+nug2-nug1));
print 'gdiff: ', np.amax(abs(nvg+nvg2-nvg1));
print 'pdiff: ', np.amax(abs(npg+npg2-npg1));
print 'vtdiff: ', np.amax(abs(nvtg+nvtg2-nvtg1));

#test.shutDownModel();
