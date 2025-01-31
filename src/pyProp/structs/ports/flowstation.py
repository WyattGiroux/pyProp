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
    __slots__ = ('static', 'total', 'ondesign', 'MN', 'V', 'A', 'mdot')
    
    def __init__(self, name, species=Air, Tref=298.15, Pref=101325):
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
    
    
    def evaluate(self):
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
        self.total.set_TP(Tt, Pt)
        if self.verifyDefined():
            self.evaluate()
        
        
    def setTotal_hP(self, ht, Pt):
        self.total.set_hp(ht, Pt)
        if self.verifyDefined():
            self.evaluate()
        
        
    def setTotal_hs(self, ht, s):
        self.total.set_hs(ht, s)
        if self.verifyDefined():
            self.evaluate()
    
    
    def setTotal_sP(self, s, Pt):
        self.total.set_sp(s, Pt)
        if self.verifyDefined():
            self.evaluate()
        
        
    def verifyDefined(self):
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