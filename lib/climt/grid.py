
from numpy import *
from state   import KnownFields

class Grid:
    '''
    '''
    def __init__(self, Component, **kwargs):
        '''
        '''
        # Initialise dicts
        self.value = {}
        self.units = {}
        self.long_name = {}
        self.grid_layout = ['lon','lat','lev']

        # Get shape appropriate for component
        self.Shape3D = Component._getShape3D(**kwargs)

        #print 'in Grid init. Shape of ', Component.Name, ' is ', self.Shape3D

        # Levels
        self.value['nlev'] = self.Shape3D[self.grid_layout.index('lev')]
        self.long_name['lev'] = 'level'
        if Component.LevType == 'p':
            self.long_name['lev'] = 'pressure'
            self.units['lev'] = 'Pa'
        if Component.LevType is None: self.units['lev'] = '-'
        self._setAxis('lev', LevType=Component.LevType, **kwargs)

        # Latitude
        self.value['nlat'] = self.Shape3D[self.grid_layout.index('lat')]
        self.long_name['lat'] = 'latitude'
        self.units['lat'] = 'degrees'
        self._setAxis('lat', **kwargs)

        # Longitude
        self.value['nlon'] = self.Shape3D[self.grid_layout.index('lon')]
        self.long_name['lon'] = 'longitude'
        self.units['lon'] = 'degrees'
        self._setAxis('lon', **kwargs)

    def _setAxis(self, AxisName, LevType=None, **kwargs):
        '''
        Sets the value of a grid coordinate axis.
        '''
        i = self.grid_layout.index(AxisName)
        n = self.Shape3D[i]
        if AxisName in kwargs:
            self.value[AxisName] = array(kwargs[AxisName])
            # ensure axis is an array
            if self.value[AxisName].ndim == 0:
                self.value[AxisName] = array([kwargs[AxisName],])
        elif AxisName is 'lev' and LevType == 'p' and 'p' in kwargs:
            # this gets first column of p (do it like this because we don't know dims of p)
            self.value['lev'] = array(kwargs['p']).copy().flat[:n]
        else:
            if AxisName is 'lon':
                self.value[AxisName] = (arange(n)+0.5)*360./n
                if n == 1: self.value[AxisName] = array([0.])
            if AxisName is 'lat':
                self.value[AxisName] = (arange(n)+0.5)*180./n -90.
            if AxisName is 'lev':
                if LevType == 'p' :
                    self.value[AxisName] = (arange(n)+0.5)[::-1]*100000./n
                elif LevType is None:
                    self.value[AxisName] = arange(n)
                else: raise ValueError, \
                      '\n\n ++++ CliMT.Grid.init: LevType %s not recognized' % LevType
        assert n == len(self.value[AxisName]), \
               '\n\n ++++ CliMT.Grid.init: Length of input %s does not match Shape3D' % AxisName

    def __getitem__(self,key):
        try: return self.value[key]
        except: raise IndexError,'\n\n ++++ CliMT.Grid: %s not in Grid' % str(key)

    def __setitem__(self,key,value):
        self.value[key] = value

    def keys(self):
        return self.value.keys()

    def __iter__(self):
        return self.value.__iter__()
