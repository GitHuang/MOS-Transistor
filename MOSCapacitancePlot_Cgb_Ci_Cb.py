# -*- coding: utf-8 -*-
"""
Created on Fri Mar 17 15:25:08 2017

@author: e0012979
"""

# THis function is to plot the capacitance of the MOS capacitance.

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

def C_c_psi(phi_F,phi_t,gamma0,Cox,psi):  # Calculation of total bulk capacitance C_c
    den = phi_t*np.exp(-psi/phi_t) + psi - phi_t + np.exp(-2*phi_F/phi_t)*(phi_t*np.exp(psi/phi_t)-psi-phi_t)
    den = 2*np.sqrt(den)    
    num = 1-np.exp(-psi/phi_t) + np.exp(-2*phi_F/phi_t)*(np.exp(psi/phi_t)-1)
    C_c = np.sign(psi)*gamma0*Cox*num/den
    return C_c # The unit is C per cm^2
               
    
def C_b_psi(phi_F,phi_t,gamma0,Cox,psi):# Calculation of deplition capacitance C_c
    den =  psi + phi_t*np.exp((psi-2*phi_F)/phi_t)
    den = 2*np.sqrt(den)
    C_b = gamma0*Cox/den
    return C_b
    
def C_i_psi(phi_F,phi_t,gamma0,Cox,psi): # Calculation of inversion layer capacitance C_c
    den =  psi + phi_t*np.exp((psi-2*phi_F)/phi_t)
    den = 2*np.sqrt(den)
    num = np.exp((psi-2*phi_F)/phi_t)
    C_i = gamma0*Cox*num/den
    return C_i
    
#def Jacobian_Psi_VGB(Cox,C_b,C_i): # This is the function to calculte the derivative of Psi against VGB
#    J_Psi_VGB = Cox/(Cox + C_b +C_i)
#    return J_Psi_VGB

def Qc_psi(phi_F,phi_t,gamma0,Cox,psi):
    HDE = phi_t*np.exp(-psi/phi_t) + psi - phi_t + np.exp(-2*phi_F/phi_t)*(phi_t*np.exp(psi/phi_t)-psi-phi_t)
    HDE = np.sqrt(HDE)
    Qc = -np.sign(psi)*gamma0*Cox*HDE
    return Qc # The unit is C per cm^2
    
# Plot C_c over all the region
psi = np.arange(0.3*phi_F, 2*phi_F+2.4*phi_t, 0.01*phi_F)
C_c = C_c_psi(phi_F,phi_t,gamma0,Cox,psi)
C_b = C_b_psi(phi_F,phi_t,gamma0,Cox,psi)
C_i = C_i_psi(phi_F,phi_t,gamma0,Cox,psi)
C_gb = C_c*Cox/(C_c+Cox) # C_gb vs VGB
Qc= Qc_psi(phi_F,phi_t,gamma0,Cox, psi)
VGB = psi -Qc/Cox + V_FB # External Bias

fig = plt.figure()
plt.plot(psi, C_c,linewidth=2.0)
plt.plot(psi, C_b,linewidth=2.0,ls='--')
plt.plot(psi, C_i,linewidth=2.0,ls = '-.')
plt.grid()
plt.axhline(y=0, color='k')
plt.axvline(x=0, color='k')
plt.ylabel(r'$C_c (F/cm^2)$', fontsize=18, fontweight='bold')
plt.xlabel(r"$\psi_s (V)$", fontsize=18, fontweight='bold')
plt.text(phi_F, 0.0000003, r'$\phi_F$', fontsize=15)
plt.text(2*phi_F, 0.0000003, r'$2\phi_F$', fontsize=15)
plt.gcf().subplots_adjust(left=0.2)
plt.legend(['$C_c$','$C_b$','$C_i$'],loc = 'upper center',fontsize=15)
plt.savefig('bulk capacitance vs Surface Potential.pdf')



# Plot C_gb vs VGB
psi = np.arange(-4*phi_t, 2*phi_F+6*phi_t, 0.01*phi_F)
C_c = C_c_psi(phi_F,phi_t,gamma0,Cox,psi)
C_b = C_b_psi(phi_F,phi_t,gamma0,Cox,psi)
C_i = C_i_psi(phi_F,phi_t,gamma0,Cox,psi)
C_gb = C_c*Cox/(C_c+Cox) # C_gb vs VGB
Qc= Qc_psi(phi_F,phi_t,gamma0,Cox, psi)
VGB = psi -Qc/Cox + V_FB # External Bias

fig = plt.figure()
plt.plot(VGB, C_gb,linewidth=2.0)
plt.grid()
plt.axhline(y=0, color='k')
plt.axvline(x=0, color='k')
plt.ylabel(r'$C_{gb} (F/cm^2)$', fontsize=18, fontweight='bold')
plt.xlabel(r"$V_{GB} (V)$", fontsize=18, fontweight='bold')
plt.text(phi_F, 0.0000007, r'$\phi_F$', fontsize=15)
plt.text(2*phi_F, 0.0000007, r'$2\phi_F$', fontsize=15)
plt.gcf().subplots_adjust(left=0.2)
plt.savefig('total capacitance vs external bias VGB.pdf')






