cimport numpy as cnp
import numpy as np

def driver(*args):

    beta, sigma, g, dt, Cpd, tau_inf, alpha,\
        T, p, ps, solin, Ts, pint = args

    new_ps = ps[:,:,np.newaxis]

    opt_depth = (tau_inf/beta)*(1. - (pint/new_ps)**alpha)


    lwflux, upwelling, downwelling = compute_fluxes(beta,sigma,T,opt_depth,Ts)

    Tdot = (g/Cpd) * (lwflux[:,:,0:-1] - lwflux[:,:,1::])/(pint[:,:,0:-1] - pint[:,:,1::])

    Tinc = Tdot*dt

    return Tinc


cdef compute_fluxes(double beta,
                    double sigma,
                    cnp.ndarray[dtype=cnp.double_t,ndim=3] T,
                    cnp.ndarray[dtype=cnp.double_t,ndim=3] tau,
                    cnp.ndarray[dtype=cnp.double_t,ndim=2] Ts):

    nlons = tau.shape[0]
    nlats = tau.shape[1]
    nlevs = tau.shape[2]

    up = np.zeros([nlons,nlats,nlevs], order='F')

    up[:,:,0] = sigma*(Ts**4)

    dtau = tau[:,:,1::] - tau[:,:,0:-1]


    for lev in range(1,nlevs):

        up[:,:,lev] = up[:,:,lev-1]*np.exp(-dtau[:,:,lev-1])\
            + sigma*(T[:,:,lev-1]**4)*(1. - np.exp(-dtau[:,:,lev-1]))



    down = np.zeros([nlons,nlats,nlevs], order='F')

    for lev in range(nlevs-2,-1,-1):

        down[:,:,lev] = down[:,:,lev+1]*np.exp(-dtau[:,:,lev])\
            - sigma*(T[:,:,lev]**4)*(1. - np.exp(-dtau[:,:,lev]))

    up = beta*up + (1.-beta)*sigma*(Ts[:,:,np.newaxis]**4)
    down = beta*down
    lwflux = up + down

    return lwflux, up, down
