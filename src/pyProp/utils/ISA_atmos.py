##################################################################
#                                                                #
#          File: ISA Calculations                                #
#          Original Author: Prakash Prashanth                    #
#          Notes: Intended for the Intro Prop Package            #
#          https://github.mit.edu/prash/IntroPropulsion          #
#                                                                #
##################################################################

import numpy as np

#Basic constants
R = 287.05287 #Gas constant for air
gee = 9.80665  #Acc. due to gravitation
βT  = -0.0065  #Adiabatic temp lapse rate
gam = 1.4      
gmi = gam - 1

# MSL conditions
TSL = 288.15
PSL = 101325
ρSL = 1.225
aSL = 340.294
Hptrop = 11000

def atmos(Hp, form="ISA"):
    """Std. ISA properties per BADA. Uses standard ISA above 24994 m"""
    if Hp>24995:
        raise ValueError("`atmos` is only valid for altitudes below 24.995 km")
    if Hp<=Hptrop:
        T = TSL + βT*Hp  
        P = PSL*(T/TSL)**(-gee/βT/R)
    elif Hp<24994:
        T = TSL + βT*Hptrop
        Ptrop = PSL*(T/TSL)**(-gee/βT/R)
        P = Ptrop*np.exp(-gee/R/T*(Hp - Hptrop))
#     else:
#         T = (TSL + βT*Hptrop) + 0.0029892*(Hp - 24994)
#         P = 249.02*(216.65/T)**11.8

    ρ = P/R/T
    a = np.sqrt(gam*R*T)    

    return T,P,ρ,a