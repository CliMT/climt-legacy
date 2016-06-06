import numpy as np
from component import Component

from pylab import *

class held_suarez(Component):
    """
    Interface to Held-Suarez forcing for 3-D dynamical cores.

    * Instantiation

    hs = climt.held_suarez( <args> )

    where <args> include the following REQUIRED arguments:
    Name           Dims         Meaning                        Units     Default                Notes

    latitudes       1           Latitudes in the model          degrees


    and the following OPTIONAL arguments:
    Name           Dims         Meaning                        Units     Default                Notes
    
    sigbot          scalar      Height of the boundary         no units    0.7
                                layer in sigma coordinates


    lapse_rate      scalar      lapse rate of the troposphere   K          10
                                wrt sigma units

    strat_temp      scalar      Temperature of the (isothermal) K          200   
                                stratosphere

    kdrag           scalar      Drag coeffecient representing   s-1        1/86400 s-1
                                boundary layer eddy viscosity

    krada           scalar      Radiative damping coefficient   s-1         1/(40*86400) s-1

    kradb           scalar      Radiative damping coefficient   s-1         1/(4*86400) s-1

    p_surf          scalar      Reference surface pressure      Pascal      1.e5 Pa

    eq_pole_temp    scalar      Equator-Pole temperature        K           60 K
                                difference


    * Usage
    Call instance directly to compute increments

    hs( <args> )
 
    where <args> are the following REQUIRED arguments:
    Name           Dims         Meaning                        Units     Default                Notes

    U               3           zonal winds                     ms-1

    V               3           meridional winds                ms-1

    theta           3           Virtual temperature             K

    p               3           atmospheric pressure            Pa
    
    ps              2           surface pressure                Pa


    * Outputs that are accessible as hs.<Name>
    Name           Dims         Meaning                        Units     Default                Notes

    Udot            3           zonal wind tendency             ms-2

    Vdot            3           meridional wind tendency        ms-2

    thetadot        3           Virtual temperature tendency    Ks-1

    """

    def __init__(self, **kwargs):

        if 'latitudes' not in kwargs:
            raise IndexError, '\n\n latitudes is a required argument'
        else:
            self.latitudes = np.radians(kwargs.pop('latitudes'))

        if 'lapse_rate' in kwargs:
            self.lapse_rate = kwargs.pop('lapse_rate')
        else:
            self.lapse_rate = 10.

        if 'strat_temp' in kwargs:
            self.strat_temp = kwargs.pop('strat_temp')
        else:
            self.strat_temp = 200

        if 'kdrag' in kwargs:
            self.kdrag = kwargs.pop('kdrag')
        else:
            self.kdrag = 1./86400.

        if 'krada' in kwargs:
            self.krada = kwargs.pop('krada')
        else:
            self.krada = 1./(40.*86400)

        if 'kradb' in kwargs:
            self.kradb = kwargs.pop('kradb')
        else:
            self.kradb = 1./(4.*86400)

        if 'p_surf' in kwargs:
            self.p_surf = kwargs.pop('p_surf')
        else:
            self.p_surf = 1.e5

        if 'eq_pole_temp' in kwargs:
            self.eq_pole_temp = kwargs.pop('eq_pole_temp')
        else:
            self.eq_pole_temp = 60.

        if 'sigbot' in kwargs:
            self.sigbot = kwargs.pop('sigbot')
        else:
            self.sigbot = 0.7

        #Initialise the rest of the fields
        self.Name = 'held_suarez_forcing'
        self.LevType = 'p'
        self.SteppingScheme = 'explicit'
        self.ToExtension    = ['U','V','theta','q','ps','p']
        self.FromExtension  = ['Uinc','Vinc','thetainc','qinc','psinc']
        #self.FromExtension  = ['U','V','theta','q','ps','p']
        self.Required       = ['U','V','theta','q','ps','p']
        self.Diagnostic     = ['p']
        self.Prognostic     = ['U','V','theta','q','ps']

        Component.__init__(self,**kwargs)


    def driver(self, u, v, temp, q, surf_press, press, simTime=-1):

        dt = self.UpdateFreq
        lats = self.latitudes


        blprof = press/surf_press[:,:,np.newaxis]

        rad_equil_temp = (press/self.p_surf)**(2./7.)*\
            (315.-self.eq_pole_temp*np.sin(lats[np.newaxis,:,np.newaxis])**2\
             - self.lapse_rate*np.log(press/self.p_surf)*np.cos(lats[np.newaxis,:,np.newaxis])**2)

        blprof = (blprof - self.sigbot)/(1-self.sigbot)
        blprof[blprof<0] = 0

        self.rad_equil_temp = rad_equil_temp
        self.blprof = blprof

        rad_equil_temp[rad_equil_temp < self.strat_temp] = self.strat_temp

        temp_tend = (self.krada + (self.kradb - self.krada) * \
                     blprof*np.cos(lats[np.newaxis,:,np.newaxis])**4)*(rad_equil_temp - temp) * dt

        v_tend = -blprof*self.kdrag*v * dt
        u_tend = -blprof*self.kdrag*u * dt

        self.temp_tend = temp_tend
        self.u_tend = u_tend
        self.v_tend = v_tend

        lnps_tend = np.zeros(surf_press.shape, dtype=np.double, order='F')

        q_tend = np.zeros(q.shape, dtype=np.double, order='F')

        #print "In HS: ut,vt,virtt,qt,pst: ", abs(u_tend[:,:,8]).max(), abs(v_tend).max(),\
        #    abs(temp_tend).max(), abs(q_tend).max(), abs(lnps_tend).max()


        return (u_tend,v_tend,temp_tend,q_tend,lnps_tend)

