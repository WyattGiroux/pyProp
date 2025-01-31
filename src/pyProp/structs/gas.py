############################################################
#                                                          #
#     Gas Class and Functions                              #
#     Python recreation of IdealGases.jl Gas               #
#     https://mit-lae.github.io/IdealGases.jl/stable/      #
#     Author: Wyatt Giroux                                 #
#     Date: 1/10/25                                        #
#                                                          #
############################################################
import numpy as np
import bisect as bi
from copy import deepcopy

from ..data.species import *
from ..data.constants import *  
from ..utils.newton import *

class Gas:
    def __init__(self, species:Species, T, P, Tref=298.15, Pref=101325, ref0=False):
        '''
        Class representing a gaseous fluid at a specified pressure and temperature using the NASA-9
        polynomials. Currently only supports gases consisting of one entry in the `thermo.py` file.
        All units are in SI. Gas constants, enthalpy, etc. are in terms of J/kg or J/kgK.
        
        NASA-9 Report (McBride et al., 2002): https://ntrs.nasa.gov/citations/20020085330
        
        species  : gas species as defined in `src/data/thermo.py`
        T        : gas temperature [K]
        P        : gas pressure [Pa]
        Tref     : reference temperature [K]; default 298.15 K
        Pref     : reference pressure [Pa]; default 101325 Pa
        ref0     : Reference enthalpy against 0 K (True) or 298.15 K (False); default False
        '''
        self.__dict__['species'] = species
        self.__dict__['ref0'] = ref0
        
        # Identify the number of discrete temperature ranges on which the NASA-9 coefficients are
        # discretizeda and store in a sorted set
        rangeTuples = list(species.ranges.keys())
        rangeList = np.zeros(len(rangeTuples)*2)
        for i, r in enumerate(rangeTuples):
            rangeList[2*i:2*i+2] = [r[0], r[1]]
        self.__dict__['ranges'] = sorted(set(rangeList))
        
        # Store reference quantities and calculate the species-specific gas constant
        self.__dict__['Tref'] = Tref
        self.__dict__['Pref'] = Pref
        self.__dict__['R'] = Rmol / self.species.molecweight * 1000
        
        # Store the initial pressure and temperature
        self.__dict__['T'] = T
        self.__dict__['P'] = P
        
        # Update the fluid state (calculate all relevant gas properties)
        self.__updateState()
    
        
    def __updateState(self):
        ''' 
        Updates Tarray, cp, cv, g, h, phi (entropy complement), and s based on temp and pressure. 
        
        Updated Variables:
        Tarray : list of T**e values used in NASA-9 cp calculation. ln(T) is appended to the end for convenience
        cp     : Specific heat at constant pressure [J/kgK]
        cv     : Specific heat at constant volume [J/kgK]
        g      : Specific heat ratio [-]
        h      : Enthalpy [J/kg] referenced against h(298.15 K)
        phi    : Entropy complement [J/kgK]
        s      : Entropy [J/kgK]
        a      : Speed of sound [m/s]
        '''
        self.__dict__['Tarray'], self.__dict__['cp'], self.__dict__['cv'],\
        self.__dict__['g'], self.__dict__['h'], self.__dict__['phi'],\
        self.__dict__['s'], self.__dict__['a'] = self.__calcState(self.T, self.P)
        
        
    def __calcState(self, T, P):
        ''' 
        Calculates Tarray, cp, cv, g, h, phi (entropy complement), and s based on temp and pressure. 
        
        Inputs:
        T      : temperature [K]
        P      : pressure [Pa]
        
        Returns:
        Tarray : list of T**e values used in NASA-9 cp calculation. ln(T) is appended to the end for convenience
        cp     : Specific heat at constant pressure [J/kgK]
        cv     : Specific heat at constant volume [J/kgK]
        g      : Specific heat ratio [-]
        h      : Enthalpy [J/kg] referenced against h(298.15 K)
        phi    : Entropy complement [J/kgK]
        s      : Entropy [J/kgK]
        '''
        rangeKey = self.__findRange(T)
        Tarray = self.__setTarray(T)
        coef = self.species.ranges[rangeKey]['coeffs']

        dh0 = self.species.ranges[rangeKey]['dh0']
        
        cp = self.__calc_cp(Tarray, coef)
        cv = self.__calc_cv(cp)
        g = self.__calc_g(cp, cv)
        h = self.__calc_h(Tarray, coef, dh0, T)

        phi = self.__calc_phi(Tarray, coef)
        s = self.__calc_s(phi, P)
        
        a = self.__calc_a(g, T)
        return Tarray, cp, cv, g, h, phi, s, a
    
    
    def __findRange(self, T):
        """
        Determines the NASA-9 temperature range within which T [K] lies.
        If out of bounds, the NASA-9 models don't apply and an error is 
        returned.
        """
        i = bi.bisect(self.ranges, T)
        if i == 0 or i == len(self.ranges):
            raise ValueError('Temperature is out of valid NASA-9 range')
        else:
            return (self.ranges[i-1], self.ranges[i])
        
        
    def __setTarray(self, T):
        """
        Calculates list of T**e values used in NASA-9 cp calculation. ln(T) 
        is appended to the end for convenience.
        """
        rangeKey = self.__findRange(T)
        exp = self.species.ranges[rangeKey]['exp']
        Tarray_temp = [T**e for e in exp]
        Tarray_temp.append(np.log(T))
        return Tarray_temp
        
        
    def __calc_cp(self, Tarray, coef):
        """ Calculates cp given Tarray and the NASA-9 coefficients """
        cp = (coef[0]*Tarray[0] + coef[1]*Tarray[1] + coef[2] + coef[3]*Tarray[3] + \
              coef[4]*Tarray[4] + coef[5]*Tarray[5] + coef[6]*Tarray[6]) * Rmol / self.species.molecweight * 1000
        return cp
        
        
    def __calc_cv(self, cp):
        """ Calculates cv given cp and gas constant """
        return cp - self.R
        
        
    def __calc_g(self, cp, cv):
        """ Calculates specific heat ratio (cp/cv) given cp and cv """
        return cp / cv
          
                  
    def __calc_h(self, Tarray, coef, dh0, T):
        """ 
        Calculates enthalpy based on Tarray, NASA-9 coefficients, 
        h(298.15 K) - h(0 K), and temperature [K] 
        """
        shift = 0
        if self.ref0:
            shift = dh0
        h = (-coef[0]*Tarray[0]   + coef[1]*Tarray[1]*Tarray[7] + coef[2]                   + coef[3]*Tarray[3]/2 + \
                               coef[4]*Tarray[4]/3 + coef[5]*Tarray[5]/4               + coef[6]*Tarray[6]/5 + coef[7]/T) * Rmol * T /\
                               self.species.molecweight * 1000 + shift
        return h
        
        
    def __calc_phi(self, Tarray, coef):
        """ Calculates entropy complement given Tarray and NASA-9 coefficients """
        phi = (-coef[0]*Tarray[0]/2 - coef[1]*Tarray[1]   + coef[2]*Tarray[7]   + coef[3]*Tarray[3] + \
                                 coef[4]*Tarray[4]/2 + coef[5]*Tarray[5]/3 + coef[6]*Tarray[6]/4 + coef[8]) * Rmol / self.species.molecweight * 1000
        return phi
        
        
    def __calc_s(self, phi, P):
        """ Calculates entropy given entropy complement and pressure """
        s = phi - np.log(P/self.Pref) * Rmol / self.species.molecweight * 1000
        return s
    
    
    def __calc_a(self, g, T):
        """ Calculate's speed of sound [m/s] """
        return np.sqrt(g*self.R*T)
        
        
    def getvar(self, name):
        return deepcopy(self.__getattribute__(name))
        
        
    def set_TP(self, T, P):
        """ Sets gas state based on temperature and pressure """
        self.__dict__['T'] = T
        self.__dict__['P'] = P
        self.__updateState()
        return 1
        
    
    def set_h(self, hset, method=newton_relax, rtol=1e-6, atol=1e-6):
        """ Uses a non-linear newton solver to find T for a specified h """
        # Create 1x1 array of initial conditions (for compatibility w/ n-dim solver)
        x0 = np.array([self.Tref])
        
        # Resudiual function
        def f(Ttest):
            rangeKey = self.__findRange(*Ttest)
            Tarray = self.__setTarray(*Ttest)
            coef = self.species.ranges[rangeKey]['coeffs']
            dh0 = self.species.ranges[rangeKey]['dh0']
            return np.array([self.__calc_h(Tarray, coef, dh0, *Ttest) - hset])
        
        # Jacobian
        def J(Ttest):
            rangeKey = self.__findRange(*Ttest)
            Tarray = self.__setTarray(*Ttest)
            coef = self.species.ranges[rangeKey]['coeffs']
            cp = [self.__calc_cp(Tarray, coef)]
            return np.array([cp])
        
        # Choose newton solver method (baseline, relaxed) and run
        if method == newton_base:
            sol, _, _, _, _ = method(x0, f, J, reltol=rtol, abstol=atol)
        elif method == newton_relax:
            # Trust region for relaxed solver is the range covered by NASA-9
            xlow = np.array([self.ranges[0]])
            xhigh = np.array([self.ranges[-1]])
            sol, _, _, _, _ = method(x0, f, J, xlow, xhigh, reltol=rtol, abstol=atol)
        else:
            raise ValueError('Invalid solve method passed to set_h')
        
        # If the sovler converged, update the temperature and other quantities.
        # Otherwise, retain the current T
        if sol is not None:
            self.__dict__['T'] = sol[0]
            self.__updateState()
            return 1
        else:
            print('Solution not found. Keeping current T')
            return 0
        
        
    def set_s_constP(self, sset, method=newton_relax, rtol=1e-6, atol=1e-6):
        """ Uses a non-linear newton solver to T for a given s and constant pressure """
        # Create 1x1 array of initial conditions (for compatibility w/ n-dim solver)
        x0 = np.array([self.Tref])
        
        # Resudiual function
        def f(Ttest):
            rangeKey = self.__findRange(*Ttest)
            Tarray = self.__setTarray(*Ttest)
            coef = self.species.ranges[rangeKey]['coeffs']
            return np.array([self.__calc_s(self.__calc_phi(Tarray, coef), self.P) - sset])
        
        # Jacobian
        def J(Ttest):
            rangeKey = self.__findRange(*Ttest)
            Tarray = self.__setTarray(*Ttest)
            coef = self.species.ranges[rangeKey]['coeffs']
            cp_T = [self.__calc_cp(Tarray, coef) / Ttest[0]]
            return np.array([cp_T])
        
        # Choose newton solver method (baseline, relaxed) and run
        if method == newton_base:
            sol, _, _, _, _ = method(x0, f, J, reltol=rtol, abstol=atol)
        elif method == newton_relax:
            # Trust region for relaxed solver is the range covered by NASA-9
            xlow = np.array([self.ranges[0]])
            xhigh = np.array([self.ranges[-1]])
            sol, _, _, _, _ = method(x0, f, J, xlow, xhigh, reltol=rtol, abstol=atol)
        else:
            raise ValueError('Invalid solve method passed to set_h')
        
        # If the sovler converged, update the temperature and other quantities.
        # Otherwise, retain the current T
        if sol is not None:
            self.__dict__['T'] = sol[0]
            self.__updateState()
            return 1
        else:
            print('Solution not found. Keeping current T')
            return 0
            
            
    def set_hs(self, hset, sset, method=newton_relax, rtol=1e-6, atol=1e-6):
        """ Sets gas state using enthalpy and entropy """
        # Set the temperature based on specified enthalpy first
        conv = self.set_h(hset, method=method, rtol=rtol, atol=atol)
        
        # If the set_h failed, do not continue with the set
        if not conv:
            print('Cancelling s_set as well')
            return 0
        
        # Create 1x1 array of initial conditions (for compatibility w/ n-dim solver)
        x0 = np.array([self.Pref])
        
        # Residual function
        def f(Ptest):
            return np.array([self.__calc_s(self.phi, Ptest[0]) - sset])
        
        # Jacobian
        def J(Ptest):
            return np.array([-self.R / Ptest[0]])
        
        # Choose newton solver method (baseline, relaxed) and run
        if method == newton_base:
            sol, _, _, _, _ = method(x0, f, J, reltol=rtol, abstol=atol)
        elif method == newton_relax:
            # Trust region ensures no negative pressure
            xlow = np.array([0])
            xhigh = np.array([np.inf])
            sol, _, _, _, _ = method(x0, f, J, xlow, xhigh, reltol=rtol, abstol=atol)
        else:
            raise ValueError('Invalid solve method passed to set_h')
        
        # If the sovler converged, update the pressure and other quantities.
        # Otherwise, retain the current P
        if sol is not None:
            self.__dict__['P'] = sol[0]
            self.__updateState()
            return 1
        else:
            print('Solution not found. Keeping current P')
            return 0
        

    def set_hp(self, hset, Pset, method=newton_relax, rtol=1e-6, atol=1e-6):
        """ Sets gas state based on enthalpy and pressure """
        conv = self.set_h(hset, method=method, rtol=rtol, atol=atol)
        if not conv:
            print('Cancelling P_set as well')
            return 0
        
        self.__dict__['P'] = Pset
        self.__updateState()        
        
        
    def set_sp(self, sset, Pset, method=newton_relax, rtol=1e-6, atol=1e-6):
        """ Sets gas state using enthalpy and entropy """
        # Set pressure
        self.__dict__['P'] = Pset
        
        # Set the temperature based on specified entropy
        conv = self.set_s_constP(sset, method=method, rtol=rtol, atol=atol)
        if not conv:
            print('Cancelling update')
            return 0
        
        self.__updateState()
        
        
    def __setattr__(self, name, val):
        """ Allows the user to set T and P directly with Gas.`name` = `val` """
        if name not in ['T', 'P']:
            raise TypeError('Can only set T and P')
        self.__dict__[name] = val
        self.__updateState()