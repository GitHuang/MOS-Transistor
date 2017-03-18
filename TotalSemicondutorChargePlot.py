# -*- coding: utf-8 -*-
"""
Created on Wed Mar 15 18:12:02 2017

@author: e0012979
"""
# Problem of Operation and Modeling of MOS Transistors Page 112, Problem 2.3


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
    return Qc # The unit is C per cm^2
    

psi = np.arange(-5*phi_t, 2*phi_F+6*phi_t, 0.01*phi_F)
Qc= Qc_psi(phi_F,phi_t,gamma0,Cox, psi)
# plot the the Total semiconductor charge vs. surface potential and constrint imposed by the rest of the system.
fig = plt.figure()
#fig.suptitle('Bulk Charge vs Surface Potential', fontsize=14, fontweight='bold')

plt.plot(psi,Qc)
# Here we link the external bias VGB to the surface potential
VGB = psi -Qc/Cox + V_FB # External Bias
Vext_flat = VGB - V_FB
Qc_line = Cox*(psi-Vext_flat) 
plt.plot(psi,Qc_line,'k--')
plt.grid()
plt.axhline(y=0, color='k')
plt.axvline(x=0, color='k')
plt.xlabel(r'$\psi_s (V)$', fontsize=18, fontweight='bold')
plt.ylabel(r"$Q_C (C/cm^2)$", fontsize=18, fontweight='bold')
plt.text(phi_F, 0.000001, r'$\phi_F$', fontsize=15)
plt.text(2*phi_F, 0.000001, r'$2\phi_F$', fontsize=15)
plt.text(2.2*phi_F, 0.000003, r'$V_{GB}-V_{FB}$', fontsize=15)
plt.gcf().subplots_adjust(left=0.2)
plt.savefig('Bulk Charge vs Surface Potential.pdf')


# Here we plot VGB vs psi
Vext = psi -Qc/Cox
fig = plt.figure()
plt.plot(Vext, psi)
#fig.suptitle('VGB vs VGB-VFB', fontsize=14, fontweight='bold')
plt.grid()
plt.axhline(y=0, color='k')
plt.axvline(x=0, color='k')
plt.ylabel(r'$\psi_s (V)$', fontsize=18, fontweight='bold')
plt.xlabel(r"$V_{GB}-V_{FB}(V)$", fontsize=18, fontweight='bold')
plt.text(0.01,phi_F,  r'$\phi_F$', fontsize=15)
plt.text( 0.01,2*phi_F, r'$2\phi_F$', fontsize=15)
plt.gcf().subplots_adjust(left=0.13)
plt.savefig('External Bias vs Surface Potential.pdf')


# Here we plot Qc vs psi
fig = plt.figure()
plt.plot(Vext_flat, Qc)
#fig.suptitle('Qc vs VGB-VFB', fontsize=14, fontweight='bold')
plt.grid()
plt.axhline(y=0, color='k')
plt.axvline(x=0, color='k')
plt.ylabel(r'$Q_c (C/cm^2)$', fontsize=18, fontweight='bold')
plt.xlabel(r"$V_{GB}-V_{FB}(V)$", fontsize=18, fontweight='bold')
plt.gcf().subplots_adjust(left=0.2)
plt.savefig('Bulk Charge vs External Bias.pdf')
























