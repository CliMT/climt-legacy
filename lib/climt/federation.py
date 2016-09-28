#!/usr/bin/env python

from types      import *
from numpy import *
from component  import Component
from parameters import Parameters
from state      import State
from plot       import Monitor, Plot
from inout      import IO
from _timestep  import asselin
import _grid

class federation(Component):
    """
    Combine components to create time-dependent model.

    * Instantiation:
    
      x = climt.federation(C1, C2, ..., Cn, <args> )

      where C1, ..., Cn are instances of CliMT components 
      and <args> are any keywork arguments relevant to the
      components included.

      The parameters and state of constituent components is re-initialized.

    * Running the federation:

      x.step()     will evolve the federation 1 timestep
      x.step(100)  will evolve the federation 100 timesteps
      x.step(100.) will evolve the federation 100 seconds
    """
    def __init__(self, *components, **kwargs):
        """
        """
        # Check input components
        if len(components) < 2:
            raise \
        '\n\n +++ CliMT.federation: you must give me more than 1 component to federate!\n\n'
        else:
            for component in components:
                assert isinstance(component, Component), \
                '\n\n +++CliMT.federation: Input item %s is not an instance.\n\n' % str(component) 
                     
        # Re-order components: diagnostic, semi-implicit, explicit, implicit
        components = list(components)
        """
        for i in range(len(components)):
            if len(components[i].Prognostic) > 0:
                components.append(components.pop(i))
        for scheme in ['semi-implicit', 'explicit', 'implicit']:
            for i in range(len(components)):
                if components[i].SteppingScheme == scheme:
                    components.append(components.pop(i))
        """
        self.components = components
        
        # Federation's Required is union of all components' Required;
        # same for Prognostic and Diagnostic
        self.Required    = []
        self.Prognostic  = []
        self.Diagnostic  = []
        for component in components:
            self.Required   = list(set(self.Required).union(component.Required))
            self.Prognostic = list(set(self.Prognostic).union(component.Prognostic))
            self.Diagnostic = list(set(self.Diagnostic).union(component.Diagnostic))

        #Check if any components carry an integrator
        self.Integrators = []
        self.Integrates = []
        self.FromExtension = []
        self.CanIntegrate = False
        for component in components:
            if component.CanIntegrate:
                self.CanIntegrate = component.CanIntegrate

                common_fields = set(component.Integrates).intersection(self.Integrates)
                if common_fields:
                    # Two components are trying to integrate the same field.
                    # don't allow
                    raise IndexError, "\n\n Two components are trying to integrate the fields ",\
                                                     common_fields

                for field in component.Integrates:
                    self.Integrates.append(field)

                for field in component.FromExtension:
                    self.FromExtension.append(field)

                self.Integrators.append(component)
                print component.Name, ' can integrate ', component.Integrates

        if self.CanIntegrate:
            print 'All fields integrated by federation members: ', self.Integrates
            print 'All fields returned by federation members: ', self.FromExtension

        # Other attributes
        self.Name      = 'federation'
        self.Extension = None

        # Set LevType to None if all components are None, else p
        self.LevType  = None         
        for component in components:
            if component.LevType == 'p': self.LevType = 'p'

        # Initialize self.Fixed (subset of self.Prognostic which will NOT be time-marched)
        if 'Fixed' in kwargs: self.Fixed = kwargs.pop('Fixed')
        else: self.Fixed = []

        # Instantiate I/O
        self.Io = IO(self, **kwargs)

        # Get values from restart file, if available
        if 'RestartFile' in kwargs:
            ParamNames = Parameters().value.keys()
            FieldNames = self.Required
            kwargs = self.Io.readRestart(FieldNames, ParamNames, kwargs)

        # Initialize scalar parameters
        self.Params  = Parameters(**kwargs)

        # Initialize State
        self.State = State(self, **kwargs)
        self.Grid = self.State.Grid
        if 'grid' in kwargs:
            self.Grid = kwargs.pop('grid')

        # Set some redundant attributes (mainly for backward compatibility)
        self.nlon = self.Grid['nlon']
        self.nlat = self.Grid['nlat']
        self.nlev = self.Grid['nlev']
        try: self.o3 = self.State['o3']
        except: pass

        self.componentGrids = []
        # Check if components enforce axis dimensions, ensure consistency
        for component in self.components:
            if component.Grid not in self.componentGrids:
                self.componentGrids.append(component.Grid)
            for AxisName in ['lev','lat','lon']:
                exec('n_fed = self.n%s' % AxisName)
                try: exec('n_com = component.Extension.get_n%s()' % AxisName)
                except: n_com = n_fed
                assert n_com == n_fed, \
                '\n\n ++++ CliMT.federation.init: recompile with %i %ss to run this federation\n'\
                % (n_fed,AxisName)

        # Dictionary to hold increments on prognos fields
        # We need three increments for a third order Adams-Bashforth
        self.Inc = {}
        self.IncOld = {}
        self.IncOlder = {}

        # Adjust components' attributes
        for component in self.components:
            component.Monitoring = False
            component.Io.OutputFreq = self.Io.OutputFreq
            component.Fixed.extend(self.Fixed)
            if component.UpdateFreq == component['dt']: component.UpdateFreq = self['dt']
            component.Params     = self.Params
            component.Grid       = self.State.Grid 
            component.State      = self.State 
            component.Inc        = {}
            # insolation component gets special treatment because 
            # of need to set orb params in common block (yes, this is ugly)
            try: component.setOrbParams(**kwargs)
            except: pass
        self.compute(ForcedCompute=True)

        # Create output file
        self.Io.createOutputFile(self.State, self.Params.value)
        
        # Write out initial state
        if not self.Io.Appending: self.write()
        
        # Initialize plotting facilities
        self.Plot = Plot()

        # Initialize runtime monitor
        self.Monitor = Monitor(self,**kwargs)

        # Notify user of unused input quantities
        self._checkUnused(kwargs)
                
        # Print out report
        #self.report()

    def compute(self, ForcedCompute=False):
        """
        Update federation's diagnostics and increments.
        """
        ## New = self.State.Old.copy()
        ## for component in self.components:            
        ##     # enforce time-splitting of implicit and semi-implicit components
        ##     self.State.Old.update(New)            
        ##     # bring component's diagnostics and increments up to date
        ##     component.compute(ForcedCompute=ForcedCompute)
        ##     # accumulate increments
        ##     for key in component.Inc:
        ##         New[key] += component.Inc[key]
        ## for key in self.State.Old:
        ##     self.Inc[key] = New[key]  - self.State.Old[key]

        self.Inc = self.State.Old.copy()
        for key in self.Inc: self.Inc[key] = self.Inc[key]*0.
        for component in self.components:            
            # bring component's diagnostics and increments up to date
            component.compute(ForcedCompute=ForcedCompute)
            # accumulate increments
            for key in component.Inc:
                self.Inc[key] += component.Inc[key]

    def integrate(self, field_list, increment_list):

        Input = []
        InputTend = []
        Output = list(self.FromExtension)

        for component in self.Integrators:

            #print 'In Fed: ', component.Integrates

            for field in component.Integrates:

                index = self.Integrates.index(field)

                Input.append(field_list[index])
                InputTend.append(increment_list[index])

            OutputValues = component.integrate(Input, InputTend)

            for field in component.FromExtension:

                fed_index = self.FromExtension.index(field)
                comp_index = component.FromExtension.index(field)
                #print 'In Fed: ', field, fed_index, comp_index

                Output[fed_index] = OutputValues[comp_index]

        return Output
