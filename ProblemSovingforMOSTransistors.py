# -*- coding: utf-8 -*-
"""
Created on Sat Mar 18 10:48:50 2017

@author: e0012979
"""

# Problem set 2 for the MOS transistor course
# Problem1
import numpy as np

N_A = 5E17
n_i = 1.5E10
phi_t = 0.0259 # T=300 K, the thermal voltage is 0.0259V
phi_poly = -0.56 # degenerated n-type poly gate
Qo = 2E-9   # Dieletric trapped charge C/cm^2 
kox = 3.9 # dieletric for SiO2
epsilon0 = 8.854E-3  # fF/um
tox = 2 # nm
Cox = kox*epsilon0/tox*10E-4 # the unit is F/cm^2, ie. F per cm^2


phi_F = phi_t*np.log(N_A/n_i)  # Fermi potential of the p-body
phi_MS = phi_poly - phi_F
V_FB = phi_MS - Qo/Cox

New_ni = N_A*np.exp(-0.477/0.0259)