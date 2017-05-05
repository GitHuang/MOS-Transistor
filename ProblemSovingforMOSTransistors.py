# -*- coding: utf-8 -*-
"""
Created on Sat Mar 18 10:48:50 2017

@author: e0012979
"""
#get_ipython().magic('reset -sf') # clear variable explorer

# Problem set 2 for the MOS transistor course
# Problem1
get_ipython().magic('reset -sf')
import numpy as np

N_A = 5E17
n_i = 1.5E10
phi_t = 0.0259 # T=300 K, the thermal voltage is 0.0259V
phi_poly = -0.56 # degenerated n-type poly gate
Qo = 2E-9   # Dieletric trapped charge C/cm^2 
kox = 3.9 # dieletric for SiO2
epsilon0 = 8.854E-14  # F/cm
tox = 2 # nm
Cox = kox*epsilon0/tox*1E7 # the unit is F/cm^2, ie. F per cm^2
gamma0 = 0.53*(tox/10)*np.sqrt(N_A*1E-17) # Body effect gammma0 = sqrt(2q*epsilon*N_A)/Cox


#N_A = 1E18
#n_i = 1.5E10
#phi_t = 0.0259 # T=300 K, the thermal voltage is 0.0259V
#phi_poly = -0.56 # degenerated n-type poly gate
#Qo = 1E-8   # Dieletric trapped charge C/cm^2 
#kox = 3.9 # dieletric for SiO2
#epsilon0 = 8.854E-14  # F/cm
#tox = 2 # nm
#Cox = kox*epsilon0/tox*1E7 # the unit is F/cm^2, ie. F per cm^2
#gamma0 = 0.53*(tox/10)*np.sqrt(N_A*1E-17) # Body effect gammma0 = sqrt(2q*epsilon*N_A)/Cox



phi_F = phi_t*np.log(N_A/n_i)  # Fermi potential of the p-body
phi_MS = phi_poly - phi_F
V_FB = phi_MS - Qo/Cox

# Problem 6
phi_0 = 2*phi_F+6*phi_t  # maximum surface potential
VT0 = V_FB + (phi_0) + gamma0*np.sqrt(phi_0)  # external bias VGB for thresh hold voltage


