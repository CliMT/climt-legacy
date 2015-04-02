cimport numpy as cnp;
import numpy as np;
from numpy import empty;

#ctypedef cnp.ndarray[complex, ndim=3] Cplx3d
#ctypedef cnp.ndarray[complex, ndim=2] Cplx2d
#ctypedef cnp.ndarray[complex, ndim=1] Cplx1d

#ctypedef cnp.ndarray[double, ndim=4] Real4d
#ctypedef cnp.ndarray[double, ndim=3] Real3d
#ctypedef cnp.ndarray[double, ndim=2] Real2d
'''
cdef extern: 
    void initialiseSpectralArrays(\
        cnp.ndarray[complex, ndim=2] *pyVrtSpec, \
        cnp.ndarray[complex, ndim=2] *pyDivSpec,\
        cnp.ndarray[complex, ndim=3] *pyTracerSpec,\
        cnp.ndarray[complex, ndim=2] *pyVirtTempSpec,\
        cnp.ndarray[complex, ndim=1] *pyTopoSpec,\
        cnp.ndarray[complex, ndim=1] *pyLnPsSpec,\
        cnp.ndarray[complex, ndim=1] *pyDissSpec,\
        cnp.ndarray[complex, ndim=1] *pyDmpProf,\
        cnp.ndarray[complex, ndim=1] *pyDiffProf)

cdef extern: 
    void initialiseGridArrays(\
        cnp.ndarray[double, ndim=3] pyUg,\
        cnp.ndarray[double, ndim=3] pyVg,\
        cnp.ndarray[double, ndim=3] pyVrtg,\
        cnp.ndarray[double, ndim=3] pyDivg,\
        cnp.ndarray[double, ndim=3] pyVirtTempg,\
        cnp.ndarray[double, ndim=4] pyTracerg,\
        cnp.ndarray[double, ndim=3] pyDlnpdtg,\
        cnp.ndarray[double, ndim=3] pyEtaDotg,\
        cnp.ndarray[double, ndim=2] pyLnPsg,\
        cnp.ndarray[double, ndim=2] pyPhis,\
        cnp.ndarray[double, ndim=2] pyDPhisdx,\
        cnp.ndarray[double, ndim=2] pyDPhisdy,\
        cnp.ndarray[double, ndim=2] pyDlnpsdt)
'''

# Typedef for function pointer returning void and taking no arguments (for now)

ctypedef void (*pyPhysicsCallback)()

# Variables for grid sizes from the library
cdef extern:
    int nlats, nlons, nlevs, ntrunc, ndimspec, ntrac

# Variables for simulation time control
cdef extern:
    int fhmax, fhout

# time step and full time from library
cdef extern:
    double dt, deltim, t

# Variables that control the behaviour of the model
# bint is Cython's boolean type
cdef extern:
    bint dry, adiabatic, heldsuarez

#Variable for initial pressure
cdef extern:
    double pdryini

# Function definitions to do our work
cdef extern:
    void gfsReadNamelist()

# Function to init dynamics
cdef extern:
    void gfsInitDynamics()

# Function to init physics
cdef extern:
    void gfsInitPhysics()

#Function to calculate pressure fields after other fields have been updated
cdef extern:
    void gfsCalcPressure()

# Function to step fields by one dt
cdef extern:
    void gfsTakeOneStep()

#Function to register an external callback
cdef extern:
    void gfsRegisterPhysicsCallback(pyPhysicsCallback callback)

#Function to deallocate arrays in physics, etc.,

cdef extern:
    void gfsFinalise()

#Function to convert the input grid arrays to spectral arrays
# The data must be available in py{Ug,Vg,VirtTempg,Tracerg,Lnpsg}
cdef extern:
    void gfsConvertToSpec()

#Function to convert the spectral arrays to grid arrays
# The result will be available in py{Ug,Vg,VirtTempg,Tracerg,Lnpsg}
cdef extern:
    void gfsConvertToGrid()

#Function to initialise the arrays for computation
# These are allocated in Python and passed in
cdef extern: 
    void initialiseSpectralArrays(\
        complex *pyVrtSpec, \
        complex *pyDivSpec,\
        complex *pyVirtTempSpec,\
        complex *pyTracerSpec,\
        complex *pyTopoSpec,\
        complex *pyLnPsSpec,\
        complex *pyDissSpec,\
        complex *pyDmpProf,\
        complex *pyDiffProf)

cdef extern: 
    void initialiseGridArrays(\
        double *pyUg,\
        double *pyVg,\
        double *pyVrtg,\
        double *pyDivg,\
        double *pyVirtTempg,\
        double *pyTracerg,\
        double *pyDlnpdtg,\
        double *pyEtaDotg,\
        double *pyLnPsg,\
        double *pyPhis,\
        double *pyDPhisdx,\
        double *pyDPhisdy,\
        double *pyDlnpsdt)

cdef extern:
    void initialisePressureArrays(\
        double *pySurfPressure,\
        double *pyPressGrid)

cdef void testFunc():
    print 'a'


cdef class _gfs_dynamics:

# Grid space arrays, cython for initialisation
    cdef public cnp.double_t[::1,:,:,:] pyTracerg
    cdef public cnp.double_t[::1,:,:] pyUg, pyVg, pyVrtg, pyDivg,\
               pyVirtTempg, pyDlnpdtg, pyEtaDotg
    cdef public cnp.double_t[::1,:] pyLnPsg, pyPhis, pyDPhisdx,\
               pyDPhisdy, pyDlnpsdt
    
# Pressure arrays, in grid space
    cdef public cnp.double_t[::1,:] pySurfPressure
    cdef public cnp.double_t[::1,:,:] pyPressGrid


# Spectral arrays, using cython for declaration
    cdef public cnp.complex_t[::1,:,:] pyTracerSpec
    cdef public cnp.complex_t[::1,:] pyVrtSpec, pyDivSpec, pyVirtTempSpec
    cdef public cnp.complex_t[:] pyTopoSpec, pyLnPsSpec, \
                pyDissSpec, pyDmpProf, pyDiffProf

# Grid size
    cdef public int numLats, numLons, numTrunc, numLevs, spectralDim, numTracers

# Model state
    cdef int modelIsInitialised

    def __init__(self, numLons=192, numLats=94, \
                    simTimeHours=24, timestep=1200.0,\
                    useHeldSuarez=True, dryPressure=1e5,\
                    numTracers=1):

        global adiabatic, dry, nlats, nlons, nlevs,\
            ntrunc, ndimspec, ntrac, fhmax, deltim,\
            heldsuarez, dt, pdryini, fhout
# Read in the namelist
#        gfsReadNamelist();


# Read in the grid sizes (mainly to improve readability)
        if(numLats):
            nlats = <int>numLats;
        
        if(numLons):
            nlons = <int>numLons;
            ntrunc = <int>numLons/3 - 2;
            ndimspec = (ntrunc+1)*(ntrunc+2)/2;

        if(timestep):
            deltim = <double>timestep;
            dt = <double>timestep;


        nlevs = 28;
        pdryini = <double> dryPressure;
        ntrac = <int>numTracers;
        heldsuarez = useHeldSuarez;
        fhmax = 9600;
        fhout = 24;
        adiabatic = False;
        dry = True;

        self.numLats = nlats;
        self.numLons = nlons;
        self.numLevs = nlevs;
        self.numTrunc = ntrunc;
        self.spectralDim = ndimspec;
        self.numTracers = ntrac;

        print 'Lats, lons, levs, trunc, dims, tracers', nlats, nlons, nlevs,\
            ntrunc, ndimspec, ntrac;


# method to reconfigure model after instantiation. HAVE to call init model
# afterwards            
    def configureModel(self, numLons=None, numLats=None):
        
        global adiabatic, dry, nlats, nlons, nlevs, ntrunc, ndimspec, ntrac

        if(numLats):
            self.numLats = <int>numLats;
            nlats = self.numLats;

        if(numLons):
            print self.numLons;
            self.numLons = <int>numLons;
            nlons = self.numLons;
            self.numTrunc = self.numLons/3;
            self.spectralDim = (self.numTrunc+1)*(self.numTrunc+2)/2;
            
            ntrunc = self.numTrunc;
            ndimspec = self.spectralDim;

        print 'Current Lats, lons, trunc, dims', nlats, nlons, ntrunc, ndimspec
        print 'Current Lats, lons, trunc, dims', self.numLats, self.numLons,\
        self.numTrunc, self.spectralDim

# Initialise arrays and dynamics and physics
    def initModel(self):

        if(self.modelIsInitialised):
            self.shutDownModel();

        self.initSpectralArrays();
        self.initGridArrays();
        self.initPressureArrays();

# Now that the arrays are initialised, call dynamics and physics init

        gfsInitDynamics();
        gfsInitPhysics();
        self.modelIsInitialised = 1;

# Create the spectral arrays (defined in spectral_data.f90)

    def initSpectralArrays(self):
        global adiabatic, dry, nlats, nlons, nlevs, ntrunc, ndimspec, ntrac

        self.pyTracerSpec = np.zeros((ndimspec, nlevs, ntrac),dtype=complex, order='F');

#self.pyVrtSpec = np.zeros((ndimspec, nlevs),dtype=complex, order='F');
        self.pyVrtSpec = np.array(np.arange(ndimspec*nlevs).reshape(ndimspec, nlevs),dtype=complex, order='F');
        self.pyVirtTempSpec = np.zeros((ndimspec, nlevs),dtype=complex, \
                order='F');
        self.pyDivSpec = np.zeros((ndimspec, nlevs),dtype=complex,order='F');

        
        self.pyTopoSpec = np.zeros(ndimspec,dtype=complex);
        self.pyLnPsSpec = np.zeros(ndimspec,dtype=complex);
        self.pyDissSpec = np.zeros(ndimspec,dtype=complex);
        
        self.pyDmpProf = np.zeros(nlevs,dtype=complex);
        self.pyDiffProf = np.zeros(nlevs,dtype=complex);

        if(ntrac > 0):
            initialiseSpectralArrays(\
                 <double complex *>&self.pyVrtSpec[0,0], \
                 <double complex *>&self.pyDivSpec[0,0],\
                 <double complex *>&self.pyVirtTempSpec[0,0],\
                 <double complex *>&self.pyTracerSpec[0,0,0],\
                 <double complex *>&self.pyTopoSpec[0],\
                 <double complex *>&self.pyLnPsSpec[0],\
                 <double complex *>&self.pyDissSpec[0],\
                 <double complex *>&self.pyDmpProf[0],\
                 <double complex *>&self.pyDiffProf[0]);
        else:
            initialiseSpectralArrays(\
                 <double complex *>&self.pyVrtSpec[0,0], \
                 <double complex *>&self.pyDivSpec[0,0],\
                 <double complex *>0,\
                 <double complex *>&self.pyVirtTempSpec[0,0],\
                 <double complex *>&self.pyTopoSpec[0],\
                 <double complex *>&self.pyLnPsSpec[0],\
                 <double complex *>&self.pyDissSpec[0],\
                 <double complex *>&self.pyDmpProf[0],\
                 <double complex *>&self.pyDiffProf[0]);


# Create the grid arrays (defined in grid_data.f90)

    def initGridArrays(self):
        global adiabatic, dry, nlats, nlons, nlevs, ntrunc, ndimspec, ntrac


        self.pyTracerg = np.zeros((nlons, nlats, nlevs, ntrac),\
                dtype=np.double, order='F');

        self.pyUg = np.zeros((nlons, nlats, nlevs), dtype=np.double, order='F');
        self.pyVg = np.zeros((nlons, nlats, nlevs), dtype=np.double, order='F');
        self.pyVrtg = np.zeros((nlons, nlats, nlevs), dtype=np.double, order='F');
        self.pyDivg = np.zeros((nlons, nlats, nlevs), dtype=np.double, order='F');
        self.pyVirtTempg = np.zeros((nlons, nlats, nlevs), dtype=np.double, order='F');
        self.pyDlnpdtg = np.zeros((nlons, nlats, nlevs), dtype=np.double, order='F');
        self.pyEtaDotg = np.zeros((nlons, nlats, nlevs+1), dtype=np.double, order='F');


        self.pyLnPsg = np.zeros((nlons, nlats), dtype=np.double, order='F');
        self.pyPhis = np.zeros((nlons, nlats), dtype=np.double, order='F');
        self.pyDPhisdx = np.zeros((nlons, nlats), dtype=np.double, order='F');
        self.pyDPhisdy = np.zeros((nlons, nlats), dtype=np.double, order='F');
        self.pyDlnpsdt = np.zeros((nlons, nlats), dtype=np.double, order='F');

        if(ntrac > 0):
            initialiseGridArrays(\
                 <double *>&self.pyUg[0,0,0],\
                 <double *>&self.pyVg[0,0,0],\
                 <double *>&self.pyVrtg[0,0,0],\
                 <double *>&self.pyDivg[0,0,0],\
                 <double *>&self.pyVirtTempg[0,0,0],\
                 <double *>&self.pyTracerg[0,0,0,0],\
                 <double *>&self.pyDlnpdtg[0,0,0],\
                 <double *>&self.pyEtaDotg[0,0,0],\
                 <double *>&self.pyLnPsg[0,0],\
                 <double *>&self.pyPhis[0,0],\
                 <double *>&self.pyDPhisdx[0,0],\
                 <double *>&self.pyDPhisdy[0,0],\
                 <double *>&self.pyDlnpsdt[0,0]);
        else:
            initialiseGridArrays(\
                 <double *>&self.pyUg[0,0,0],\
                 <double *>&self.pyVg[0,0,0],\
                 <double *>&self.pyVrtg[0,0,0],\
                 <double *>&self.pyDivg[0,0,0],\
                 <double *>&self.pyVirtTempg[0,0,0],\
                 <double *>0,\
                 <double *>&self.pyDlnpdtg[0,0,0],\
                 <double *>&self.pyEtaDotg[0,0,0],\
                 <double *>&self.pyLnPsg[0,0],\
                 <double *>&self.pyPhis[0,0],\
                 <double *>&self.pyDPhisdx[0,0],\
                 <double *>&self.pyDPhisdy[0,0],\
                 <double *>&self.pyDlnpsdt[0,0]);

#Intialise pressure arrays
    def initPressureArrays(self):
        global adiabatic, dry, nlats, nlons, nlevs, ntrunc, ndimspec, ntrac

        self.pySurfPressure = np.zeros((nlons, nlats), dtype=np.double, order='F');
        self.pyPressGrid = np.zeros((nlons, nlats, nlevs), dtype=np.double, order='F');

        initialisePressureArrays(\
                <double *>&self.pySurfPressure[0,0],\
                <double *>&self.pyPressGrid[0,0,0]);


# Take one step
    
    def oneStepForward(self):

        gfsTakeOneStep();

# Register a callback which calculates the physics (to be used in stand-alone
# mode only)

    def setPhysicsCallback(self,physicsFnPtr):

        gfsRegisterPhysicsCallback(testFunc);
    '''
    Does not work!
    def setInitialConditions(self, inputList):

        myug,myvg,myvirtempg,mytracerg,mylnpsg = inputList;

        self.pyUg[:] = myug[:];
        self.pyVg[:] = myvg[:];
        self.pyVirtTempg[:] = myvirtempg[:];
        self.pyTracerg[:] = mytracerg[:];
        self.pyLnPsg[:] = mylnpsg[:];

        gfsConvertToSpec();
    '''
    def getResult(self):

        gfsConvertToGrid();
        gfsCalcPressure();
        
        outputList = [];

        outputList.append(np.asarray(self.pyUg).copy(order='F'));
        outputList.append(np.asarray(self.pyVg).copy(order='F'));
        outputList.append(np.asarray(self.pyVirtTempg).copy(order='F'));
        outputList.append(np.asarray(self.pyTracerg[:,:,:,0]).copy(order='F'));
        outputList.append(np.asarray(self.pyLnPsg).copy(order='F'));
        outputList.append(np.asarray(self.pyPressGrid).copy(order='F'));

        return(outputList);

# method to override the parent class (Component) method (to be used in CliMT
# mode only)
    def driver(self, myug, myvg, myvirtempg, myqg, mylnpsg, mypress, double simTime=-1.):
       
        
        cdef cnp.double_t[::1,:,:] tempug,tempvg,tempvtg
        cdef cnp.double_t[::1,:,:] tempqg
        cdef cnp.double_t[::1,:] templnpsg


        if(simTime >= 0):
            global t
            t = simTime;

#myug,myvg,myvirtempg,myqg,mylnpsg,mypress = inputArgs;

        myug = np.asfortranarray(myug);
        myvg = np.asfortranarray(myvg);
        myvirtempg = np.asfortranarray(myvirtempg);
        myqg = np.asfortranarray(myqg);
        mylnpsg = np.asfortranarray(mylnpsg);

#Convert to memory view so that assignment can be made to arrays
        tempug = myug;
        tempvg = myvg;
        tempvtg = myvirtempg;

        tempqg = myqg;
        templnpsg = mylnpsg;

#Assign to model arrays
        self.pyUg[:] = tempug;
        self.pyVg[:] = tempvg;
        self.pyVirtTempg[:] = tempvtg;
        self.pyTracerg[:,:,:,0] = tempqg;
        self.pyLnPsg[:] = templnpsg;

#Convert to spectral space
        gfsConvertToSpec();

#Step forward in time
        self.oneStepForward();

#Convert back to grid space
        gfsConvertToGrid();

# only ln(Ps) is calculated in the dynamics. This calculates
# the values on the full grid
        gfsCalcPressure();

        '''
        ugInc = self.pyUg - myug;
        vgInc = self.pyVg - myvg;
        virtempgInc = self.pyVirtTempg - myvirtempg;
        tracergInc = self.pyTracerg[:,:,:,0] - myqg;
        lnpsgInc = self.pyLnPsg - mylnpsg;
        pressInc = self.pyPressGrid - mypress;

        return(ugInc,vgInc,virtempgInc,tracergInc,lnpsgInc,pressInc);
        '''
        ug = self.pyUg.copy();
        vg = self.pyVg.copy();
        virtempg = self.pyVirtTempg.copy();
        tracerg = self.pyTracerg[:,:,:,0].copy();
        lnpsg = self.pyLnPsg.copy();
        press = self.pyPressGrid.copy();

        return(ug,vg,virtempg,tracerg,lnpsg,press);

    def printTimes(self):
        global dt,t;
        print 'Timestep: ',dt, 'Total Time:', t;

    def get_nlat(self):
        return self.numLats;

    def get_nlon(self):
        return self.numLons;

    def get_nlev(self):
        return self.numLevs;


    def shutDownModel(self):

        gfsFinalise();    
