# -*- coding: utf-8 -*-
"""
Created on Fri Mar 17 12:32:02 2017

@author: e0012979
"""

# -*- coding: utf-8 -*-
"""
Created on Wed Mar 15 18:12:02 2017

@author: e0012979
"""
# Problem of Operation and Modeling of MOS Transistors Page 112, Problem 2.3
# plot the total charge, bulk charge, and inversion of the MOS capacitor

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
    
psi = np.arange(0.1*phi_F, 2*phi_F+7*phi_t, 0.01*phi_F)
Qc= Qc_psi(phi_F,phi_t,gamma0,Cox, psi)
VGB = psi -Qc/Cox + V_FB # External Bias

# Calculate Bulk Charge against the surface potential, the charge sheet approximation is used
Q_B = -gamma0*Cox*np.sqrt(psi)
# Calculate inversion Charge against the surface potential, the charge sheet approximation is used
Q_I = Qc - Q_B

# Comparison of the QI QB QC against surface potential
fig = plt.figure()
plt.plot(psi, -Q_B,linewidth=2.0,linestyle = '--')
plt.plot(psi, -Q_I,linewidth=2.0)
plt.plot(psi, -Qc,linewidth=2.0)
plt.grid()
plt.axhline(y=0, color='k')
plt.axvline(x=0, color='k')
plt.ylabel(r'$Q(C/cm^2)$', fontsize=18, fontweight='bold')
plt.xlabel(r"$\psi_s (V)$", fontsize=18, fontweight='bold')
plt.text(phi_F, -0.000005, r'$\phi_F$', fontsize=15)
plt.text(2*phi_F, -0.000005, r'$2\phi_F$', fontsize=15)
plt.text(2*phi_F+6*phi_t, -0.000005, r'$\phi_0$', fontsize=15)
plt.gcf().subplots_adjust(left=0.2)
plt.legend(['$|Q_B|$', '$|Q_I|$','$|Q_C|$'],loc = 'upper center',fontsize=15)
plt.savefig('Qc QI QB vs Surface Potential.pdf')



# plot the inversion layer charge against the external bias VGB
VL0 = V_FB + phi_F + gamma0*np.sqrt(phi_F)  # external bias VGB for psi = phi_F
VM0 = V_FB + 2*phi_F + gamma0*np.sqrt(2*phi_F)  # external bias VGB for psi = 2phi_F
phi_0 = 2*phi_F+6*phi_t  # maximum surface potential
VT0 = V_FB + (phi_0) + gamma0*np.sqrt(phi_0)  # external bias VGB for thresh hold voltage
VH0 = V_FB + (phi_0) + gamma0*np.sqrt(phi_0 + phi_t* np.exp((phi_0-2*phi_F)/phi_t) )  # external bias VGB for psi = phi_0
QI_VH = -Cox*(VGB-VT0)

fig = plt.figure()
plt.plot(VGB, -Q_I,linewidth=2.0)
plt.plot([n for n in VGB if n>VT0], -QI_VH[[ n for n,i in enumerate(VGB) if i>VT0 ]] ,linewidth=2.0,ls = '--',color='r')
plt.grid()
plt.axhline(y=0, color='k')
plt.ylabel(r'$|Q_I|\ (C/cm^2)$', fontsize=18, fontweight='bold')
plt.xlabel(r"$V_{GB} (V)$", fontsize=18, fontweight='bold')
plt.text(VL0 , -0.000005, r'$V_L$', fontsize=15)
plt.text(VM0*0.8, -0.000005, r'$V_M$', fontsize=15)
plt.text(VT0*1.1, -0.000005, r'$V_{T0}$', fontsize=15)
plt.text(VH0 , -0.000005, r'$V_H$', fontsize=15)
plt.legend(['$|Q_I|(exact)$','$|Q_I|(approx)$'],loc = 'left center',fontsize=15)
plt.gcf().subplots_adjust(left=0.2)
plt.savefig('QI vs External Bias VGB.pdf')


