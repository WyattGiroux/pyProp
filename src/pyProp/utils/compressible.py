############################################################
#                                                          #
#               Compressible Flow Utility Functions        #
#               Author: Wyatt Giroux                       #
#               Date: 9/28/23                              #
#                                                          #
############################################################

import numpy as np
from .general import bisection

############################
# Computation of Constants #
############################
def cp(g, R):
    '''Specific heat at constant pressure given sp. heat ratio and gas constant'''
    return R*g/(g-1)


def cv(g, R):
    '''Specific heat at constant volume given sp. heat ratio and gas constant'''
    return R/(g-1)


########################
# Isentropic Relations #
########################
def pqPt(M, g=1.4):
    '''Isentropic relation for p/Pt as a function of mach number'''
    if M < 0: raise ValueError('Mach number must be greater than 0')
    return (1+((g-1)/2)*M**2)**(-g/(g-1))


def TqTt(M, g=1.4):
    '''Isentropic relation for T/Tt as a function of mach number'''
    if M < 0: raise ValueError('Mach number must be greater than 0')
    return (1+((g-1)/2)*M**2)**(-1)


def MfromP(PtqP, g=1.4):
    '''Isentropic relation for mach number as a function of PtqP'''
    if PtqP < 1: raise ValueError('Stagnation pressure must be greater than static pressure')
    return np.sqrt((2/(g-1))*(PtqP**((g-1)/g)-1))


def MfromT(TtqT, g=1.4):
    '''Isentropic relation for mach number as a function of TtqT'''
    if TtqT < 1: raise ValueError('Stagnation temperature must be greater than static temperature')
    return np.sqrt((2/(g-1))*(TtqT-1))


############################
# Corr. Flow Per Unit Area #
############################
def cmf(mdot, pt, tt, a):
    '''Corrected Mass Flow per Unit Area'''
    return mdot*np.sqrt(tt)/(pt*a)


def Dm(M, g=1.4, R=287):
    '''Calculates corrected mass flow per unit area given a mach number'''
    return float(np.sqrt(g/R)*M*(1+(g-1)/2*M**2)**(-1*(g+1)/(2*(g-1))))


def MfromD(D, supersonic=True, g=1.4, R=287, Mhigh=10, verbose=False):
    '''Gets mach number from corrected mass flow using bisection'''
    if D > Dm(1, g, R): raise ValueError('Given corrected mass flow exceeds maximum for given air properties')
    if supersonic: ML, MR = 1, Mhigh
    else: ML, MR = 0, 1
    return bisection(ML, MR, Dm, D, direction=supersonic, verbose=verbose)


#######################
# Basic Relationships #
#######################
def geth(T, g=1.4, R=287):
    '''Gets enthalpy given temperature and specific heat ratio'''
    return T * cp(g, R)


def ds(P,T, gam = 1.4, P0 = 101325, T0=288.15):
    '''Get entropy change from sea level conditions using Gibbs equation'''
    gmi = gam - 1
    if P/P0 < 0 or T/T0 < 0: raise ValueError('Fluid state must be real (P/P0 > 0 and T/T0 > 0)')
    return np.log(T/T0) - gmi/gam * np.log(P/P0)


def impulse(p,A,mdot,u):
    '''Gets impulse of a channel flow'''
    return p*A + mdot*u


def htrArea(D, htr):
    '''Compressor area given diameter and hub-to-tip ratio'''
    return np.pi*(D/2)**2 * (1-htr**2)


def htrAreaMean(rmean, htr):
    '''Compressor area given diameter and hub-to-tip ratio'''
    return 4 * np.pi * rmean**2 * (1-htr) / (1+htr)


#######################
# Expanded Flow Power # - OBSOLETE, IMPLEMENT WITH NASA-9 LOGIC
####################### 
# def getU(ht,h):
#     '''Calculate flow velocity based on stagnation enthalpy definition'''
#     return np.sqrt(2*(ht-h))


# def sig(mdot, V):
#     '''Ambient station power, sigma, given mass flow and velocity'''
#     return 0.5*mdot*V**2


# def sig2(mdot, ht, h):
#     '''Ambient station power, sigma, given total and static enthalpy (at ambient p)'''
#     return mdot*(ht - h)


# def gethAmb(P0, Pt, Tt, g=1.4, R=287):
#     cpi = cp(g, R)
    
#     MN = MfromP(Pt/P0, g=g)
#     Ts = TqTt(MN, g=g) * Tt
#     hs = cpi * Ts
    
#     return hs


# def getSigs(eng, mKey, htKey, hpKey): 
#     sig = {}
#     for i in eng[mKey].keys():
#         sig[i] = sig2(eng[mKey][i], eng[htKey][i], eng[hpKey][i])
    
#     return sig
                                                    
