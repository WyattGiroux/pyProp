############################################################
#                                                          #
#               Flowstation Class and Functions            #
#               Author: Wyatt Giroux                       #
#               Date: 1/10/25                              #
#                                                          #
############################################################
import numpy as np
from copy import deepcopy

from ...data.species import *
from ...utils.compressible import *  
from ...utils.newton import newton_base, newton_relax
from ..gas import Gas
from .port import Port

class FlowStation(Port):
    ''' FlowStation port
    
    Attributes:
        static (:ref:`gas`): Gas object storing static fluid properties.
        total (:ref:`gas`): Gas object storing total (a.k.a stagnation) fluid properties.
        ondesign (bool): True if on-design, False if off-design.
        MN (float): Fluid Mach number
        V (float): Fluid speed [m/s]
        A (float): Passage area [m^2]
        mdot (float): Passage mass flow [kg/s]
    '''
    __slots__ = ('static', 'total', 'ondesign', 'MN', 'V', 'A', 'mdot')
    
    def __init__(self, name, species=Air, Tref=298.15, Pref=101325):
        '''FlowStation constructor
        
        Args:
            name (str): FlowStation name.
            species (:ref:`oxy`): Fluid species in the passage. Defaults to air.
            Tref (float): Reference temperature (default 298.15 K)
            Pref (float): Reference pressure (default 101325 Pa)
        '''
        
        super().__init__(name)
        
        # STATIC AND STAGNATION GAS PROPERTIES
        self.static = Gas(species, Tref, Pref)
        self.total  = Gas(species, Tref, Pref)
        
        # On or off design? On design calculates Area; off-design calculates MN
        self.ondesign = True
        
        # PHYSICAL FLOW QUANTITIES
        self.MN   = None
        self.V    = None
        self.A    = None
        self.mdot = None
    
    
    def copy(self):
        '''Copy a FlowStation object
        
        Returns:
            A FlowStation object with identical values to the copied one.
        '''
        return deepcopy(self)
    
    
    def evaluate(self):
        ''' Updates the physical flow and static quantities.
        
        Sets either area (on-design) or Mach number (off-design) of the fluid passing through a FlowStation. Once complete
        updates the static quantities using the physical flow properties and the total fluid properties. Mass flow and the
        appropriate area or Mach number value must be set before calling this function.        
        '''
        if self.mdot is None:
            raise ValueError("mdot must be defined prior to setting A or MN")
        if self.ondesign and self.MN is not None:
            self.A = self.mdot * np.sqrt(self.Tt) / (self.Pt * Dm(self.MN, self.gt, self.Rt))
        elif not self.ondesign and self.A is not None:
            D = self.mdot * np.sqrt(self.Tt) / (self.Pt * self.A)
            Dmax = Dm(1, self.gt, self.Rt)
            
            if D > Dmax: 
                raise ValueError('Corr. Mass Flow per Unit Area is Greater Than Maximum Value') # FIXME: Extra logic is needed here to identify cases where supersonic flow is possible
                                                                                                #        Possibly just check if D = Dmax?
            
            self.MN = MfromD(D, supersonic=False, g=self.gt, R=self.Rt) 
            if abs(self.MN - 1) < 1e-6:
                self.MN = 1 # Flow is choked
        else:
            raise ValueError("Unable to calculate A or MN due to conflict with design vs. off-design state")
        
        self.__set_static()
        self.V = self.MN * self.static.a
        
    
    def __set_static(self):
        if not self.verifyDefined():
            raise ValueError("Static quantities cannot be computed without mdot and either MN or A being set")
        
        T = TqTt(self.MN, self.gt) * self.Tt
        P = pqPt(self.MN, self.gt) * self.Pt
        self.static.set_TP(T, P)
        
        
    def setTotal_TP(self, Tt, Pt):
        ''' Sets total fluid properties based on inpupt total temperature and pressure.
        
        Args:
            Tt (float): Total temperature [K]
            Pt (float): Total pressure [Pa]
        '''
        self.total.set_TP(Tt, Pt)
        if self.verifyDefined():
            self.evaluate()
        
        
    def setTotal_hP(self, ht, Pt):
        ''' Sets total fluid properties based on inpupt total enthalpy and pressure.
        
        Args:
            ht (float): Total enthalpy [J/kg]
            Pt (float): Total pressure [Pa]
        '''
        self.total.set_hp(ht, Pt)
        if self.verifyDefined():
            self.evaluate()
        
        
    def setTotal_hs(self, ht, st):
        ''' Sets total fluid properties based on inpupt total enthalpy and entropy.
        
        Args:
            ht (float): Total temperature [J/kg]
            st (float): Total entropy [J/kgK]
        '''
        self.total.set_hs(ht, st)
        if self.verifyDefined():
            self.evaluate()
    
    
    def setTotal_sP(self, st, Pt):
        ''' Sets total fluid properties based on inpupt total entropy and pressure.
        
        Args:
            st (float): Total entropy [J/kgK]
            Pt (float): Total pressure [Pa]
        '''
        self.total.set_sp(st, Pt)
        if self.verifyDefined():
            self.evaluate()
        
        
    def verifyDefined(self):
        '''Verifies whether the physical flow quantities are defined enough to be evaluated.
        
        Returns:
            True if ready for evaluation, False if not.
        '''
        if self.mdot is None:
            return 0
        if self.MN is None and self.A is None:
            return 0
        else:
            return 1
            
        
    def __getattr__(self, name):
        var = name[:-1]
        if name[-1] == 't':
            return self.total.getvar(var)
        elif name[-1] == 's':
            return self.static.getvar(var)
        else:
            raise ValueError(f'Unrecognized variable: {name}')
        
    
    def __deepcopy__(self, memo):
        new_obj = self.__class__.__new__(self.__class__)
        memo[id(self)] = new_obj
        
        for slot in self.__slots__:
            setattr(new_obj, slot, deepcopy(getattr(self, slot), memo))

        return new_obj