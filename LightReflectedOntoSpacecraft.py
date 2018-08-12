

from sympy import *
import numpy as np 



Rp = 6371*1000. #m
h = 629*1000. #m

eta = np.arccos(Rp/(Rp+h))

RG = Rp*np.cos(eta)
RK = Rp*np.sin(eta)


r_s_pHAT = np.asarray([1.,0.,0.])
#r_sc_p = np.asarray()
#r_sc_pHAT = np.asarray([np.sqrt(2)/2,-np.sqrt(2)/2,0])
r_sc_pHAT = np.asarray([0.,-1.,0.])

x = Symbol('x')
y = Symbol('y')
z = Symbol('z')
eqn1 = r_s_pHAT[0]*x + r_s_pHAT[1]*y + r_s_pHAT[2]*z
eqn2 = r_sc_pHAT[0]*x + r_sc_pHAT[1]*y + r_sc_pHAT[2]*z - RG

eqn3 = eqn2 - eqn1#line equation of two circles

# r_barHAT
# r_eHAT
# r_fHAT

# #Spacecraft Viewing Circle
# SVC = RG*r_barHAT + RK*(r_eHAT*np.cos(theta3) + r_fHAT*np.sin(theta3))

# #Sun Illumination Circle
# SIC = Rp*(s_eHAT*np.cos(theta4) + s_fHAT*np.sin(theta4))






