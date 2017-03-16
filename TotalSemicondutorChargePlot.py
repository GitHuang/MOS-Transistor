# -*- coding: utf-8 -*-
"""
Created on Wed Mar 15 18:12:02 2017

@author: e0012979
"""
# Problem of Operation and Modeling of MOS Transistors Page 112, Problem 2.3
# This function is to plot the Total semiconductor charge vs. surface potential and constrint imposed by the rest of the system.


import numpy as np
import matplotlib.pyplot as plt


phi_t = 0.0259 # T=300 K, the thermal voltage is 0.0259eV
N_A = 5E17 # the dopant density of the p-type substrate per cm^-3
n_i = 1.5E10 # the carrier density of intirnsic Si per cm^-3
tox = 2 # the thickness of the SiO2  unit: nm
V_FB = -0.8 # the flatband voltage is -0.8V
kox = 3.9 # dieletric for SiO2
epsilon0 = 8.854E-3  # fF/um


Cox = kox*epsilon0/tox*10E-4 # the unit is F/cm^2, ie. F per cm^2
gamma0 = 0.53*(tox/10)*np.sqrt(N_A*10E-17) # Body effect gammma0 = sqrt(2q*epsilon*N_A)/Cox
phi_F = phi_t*np.log(N_A/n_i)

def Qc_psi(phi_F,phi_t,gamma0,Cox,psi):
    HDE = phi_t*np.exp(-psi/phi_t) + psi - phi_t + np.exp(-2*phi_F/phi_t)*(phi_t*np.exp(psi/phi_t)-psi-phi_t)
    HDE = np.sqrt(HDE)
    Qc = -np.sign(psi)*gamma0*Cox*HDE
    return Qc   # The unit is F per cm^2
    

psi = np.arange(-0.6*phi_F, 2.7*phi_F, 0.01*phi_F)


fig = plt.figure()
fig.suptitle('Bulk Charge vs Surface Potential', fontsize=14, fontweight='bold')
plt.plot(psi,Qc_psi(phi_F,phi_t,gamma0,Cox, psi))
plt.grid()
plt.axhline(y=0, color='k')
plt.axvline(x=0, color='k')
plt.xlabel(r'$\psi_s (V)$', fontsize=18, fontweight='bold')
plt.ylabel(r"$Q_C (F/cm^2)$", fontsize=18, fontweight='bold')
plt.text(phi_F, -0.0001, r'$\phi_F$', fontsize=15)
plt.text(2*phi_F, -0.0001, r'$2\phi_F$', fontsize=15)
plt.gcf().subplots_adjust(left=0.18)
plt.savefig('Bulk Charge vs Surface Potential.pdf')
