#!/usr/bin/env python

from component  import Component
from numpy import *
import numpy as np

class dynamics(Component):
    """
    Interface to atmospheric dynamical cores. 

    * Instantiation:
    
      x=climt.dynamics( <args> )
      
      where <args> are the following OPTIONAL arguments:
      Name           Dims  Meaning                  Units     Default         Notes
      scheme            0  Dynamical core          (string) 'axisymmetric'  Choices are: 'axisymmetric'
      T               1-3  Temperature                K          283.15             
      q               1-3  Specific humidity          g/kg       1.e-5
      U               1-3  Zonal wind                 m/s        0.
      V               1-3  Meridional wind            m/s        0.
        
    * Usage:
      Call instance directly to compute dynamical tendencies.

      x( <args> )

      where <args> are as above.
        
    * Output (accessible as x.swflx etc.):
      Name       Meaning                         Units   Notes            
      Tdot       Turbulent heating rate          K s-1
      qdot       Turbulent humidification rate   g/kg s-1
      Udot       Turbulent drag                  m s-2
      Vdot       Turbulent drag                  m s-2
      SrfSenFlx  Surface sens. heat flux         W m-2
      SrfLatFlx  Surface lat. heat flux          W m-2
      taux       Surface stress                  Pa
      tauy       Surface stress                  Pa
    """
    def __init__(self, scheme = 'axisymmetric', **kwargs):
        # Initialize scheme-dependent attributes
#        if scheme not in ['axisymmetric','two_column']:
#            raise ValueError,'\n \n ++++ CliMT.dynamics: Scheme %s unknown' % scheme
#        exec('self.__%s_dynamics__init__()' % string.lower(scheme))

        if scheme == 'axisymmetric':
            self.__axisymmetric_dynamics__init__();
            # Initialize fields etc.
            Component.__init__(self, **kwargs);
        if scheme == 'two_column':
            self.__two_column_dynamics__init__();
            # Initialize fields etc.
            Component.__init__(self, **kwargs);
        elif scheme == 'gfs':
            args = self.__gfs_dynamics__init__(**kwargs)
            # Initialize fields etc.
            Component.__init__(self,**args);
        else: raise ValueError,'\n \n ++++ CliMT.dynamics: Scheme %s unknown' % scheme
        
    def __axisymmetric_dynamics__init__(self):
        # Load extension
        try: import _axisymmetric_dynamics
        except: raise ImportError, \
          '\n \n ++++ CliMT.dynamics: Could not load axisymmetric scheme'
        # Define some attributes
        self.Name           = 'axisymmetric_dynamics'
        self.LevType        = 'p'
        self.Extension      = _axisymmetric_dynamics
        self.driver         = _axisymmetric_dynamics.driver
        self.SteppingScheme = 'semi-implicit'
        self.ToExtension    = ['Rd','Cpd','r','omega','delh','delv','Newt','dt',
                               'lat','lev','T','U','V','q','Vold']
        self.FromExtension  = ['V','Tinc','Uinc','Vinc','qinc','psi','theta','Te','W',
                               'TdotDyn','UdotDyn','VdotDyn','qdotDyn']
        self.Required       = ['T','U','V','q']
        self.Prognostic     = ['T','U','V','q']
        self.Diagnostic     = ['V','psi','theta','Te','W','TdotDyn','UdotDyn','VdotDyn','qdotDyn']


    def __two_column_dynamics__init__(self):
        # Load extension
        try: import _two_column_dynamics
        except: raise ImportError, \
            '\n \n ++++ CliMT.dynamics: Could not load two-column scheme'
    # Define some attributes
        self.Name = 'two_column_dynamics'
        self.LevType = 'p'
        self.Extension = _two_column_dynamics
        self.driver = _two_column_dynamics.driver
        self.SteppingScheme = 'explicit'
        self.ToExtension = ['dt','Rd','Rv','Cpd','Cpv','g','p','z0','V','T','q']
        self.FromExtension = ['Vinc','Tinc','qinc','z0inc','W']
        self.Required = ['p','z0','V','T','q']
        self.Prognostic = ['V','T','q','z0']
        self.Diagnostic = ['W']


    def __gfs_dynamics__init__(self, **kwargs):
        # Load extension
        try:
            from _gfs_dynamics import _gfs_dynamics;

            # If you want to define your grid dimensions, it has to be done
            # somehow. For now, it reads from namelist
            if 'nlat' in kwargs:
                self.nlat = kwargs['nlat'];
            else:
                self.nlat = 94;
            if 'nlon' in kwargs:
                self.nlon = kwargs['nlon'];
            else:
                self.nlon = 192;
            if 'dt' in kwargs:
                self.timestep = kwargs['dt']
            else:
                self.timestep = 1200.0


            _gfsDycore = _gfs_dynamics(self.nlon,self.nlat, timestep=self.timestep, climt_mode=True)
            _gfsDycore.initModel()
            u,v,t,tr,ps,p,pint = _gfsDycore.getResult()

            lats = np.asarray(_gfsDycore.latitudes)[0,:]
            lons = np.asarray(_gfsDycore.longitudes)[:,0]

            args = {}

            args['lat'] = lats
            args['lon'] = lons
            args['lev'] = 100000*np.linspace(1,0.003, _gfsDycore.numLevs)

            args['p'] = p
            args['pint'] = pint
            args['ps'] = ps
            args['U'] = u
            args['V'] = v
            args['T'] = t
            args['q'] = tr
            args['CanIntegrate'] = True
            args['Integrates'] = ['U','V','T','q','ps']
            args.update(kwargs)

        except: raise ImportError, \
         '\n \n ++++ CliMT.dynamics: Could not load GFS dynamical core'
        # Define some attributes
        self.Name           = 'gfs_dynamics'
        self.LevType        = 'p'
        self.Extension      = _gfsDycore
        #self.driver         = _gfsDycore.driver
        self.integrate      = _gfsDycore.integrateFields
        self.SteppingScheme = 'explicit'
        self.ToExtension    = ['U','V','T','q','ps']
        self.FromExtension  = ['U','V','T','q','ps','p','pint']
        self.Required       = ['U','V','T','q','ps','p','pint']
        self.Diagnostic     = ['p','pint']
        self.Prognostic     = ['U','V','T','q','ps']

        return args
