# -*- coding: utf-8 -*-
"""
Created on Mon Mar 20 21:12:52 2017

@author: e0012979
"""
get_ipython().magic('reset -sf')
# THis script is to plot the I-V curve of a 4-port nMOSFET, the Ids expression is valid for all regions, and numerical methods is ultilized 
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize

Folder = r'C:\Users\e0012979\Google Drive\Projects\Git Project\MOS-Transistor\Figures/'

phi_t = 0.0259 # T=300 K, the thermal voltage is 0.0259eV
N_A = 5E17 # the dopant density of the p-type substrate per cm^-3
n_i = 1.5E10 # the carrier density of intirnsic Si per cm^-3
tox = 2.5 # the thickness of the SiO2  unit: nm
V_FB = -0.8 # the flatband voltage is -0.8V
kox = 3.9 # dieletric for SiO2
epsilon0 = 8.854E-3  # fF/um
Cox = kox*epsilon0/tox*1E-4 # the unit is F/cm^2, ie. F per cm^2
gamma0 = 0.53*(tox/10)*np.sqrt(N_A*1E-17) # Body effect gammma0 = sqrt(2q*epsilon*N_A)/Cox
phi_F = phi_t*np.log(N_A/n_i)
W = 10 # Gate Width of MOSFET : unit micro-meter
L = 10 # Length of MOSFET: micro-meter Gate
VSB = 0 # gate body bias (V)
mu = 400 # electron mobility cm^2/(V-s)
VDB = 0.7# Drain to Body bais (V)
VGB = 0 #  Gate-body voltage
markers = ['8', 's', '^', 'o', '*','v','h','8', 's', '^', 'o', '*','v','h','8', 's', '^', 'o', '*','v','h','8', 's', '^', 'o', '*','v','h']

# Define the equation to solve the surface potential at the source side
def Surface_Potential_Source_VS_ExtBias_Numerical(psi):
    temp = psi+ phi_t*np.exp((psi-2*phi_F-VSB)/phi_t)
    return VGB-V_FB-gamma0*np.sqrt(temp) -psi
    
# Define the equation to solve the surface potential at the drain side
def Surface_Potential_Drain_VS_ExtBias_Numerical(psi):
    temp = psi+ phi_t*np.exp((psi-2*phi_F-VDB)/phi_t)
    return VGB-V_FB-gamma0*np.sqrt(temp) -psi    

def Drift_Current_VS_Surface_Potential_Complete(psi_s0,psi_sL):
    Idrift = W/L*mu*Cox*(VGB-V_FB-0.5*(psi_s0 + psi_sL) -2/3*gamma0*(psi_s0+psi_sL+np.sqrt(psi_s0*psi_sL))/(np.sqrt(psi_sL) + np.sqrt(psi_s0)))*(psi_sL-psi_s0)
    Idiff = W/L*mu*Cox*phi_t*(1 + gamma0/(np.sqrt(psi_sL) + np.sqrt(psi_s0)))*(psi_sL-psi_s0)
    IDS = Idrift + Idiff
    return Idrift, Idiff, IDS


VGB1 = np.arange(0.2, 0.9, 0.05) # range of Gate-Body bias
Idrift = []
Idiff = []
IDS = []
psi_s0_array = []
psi_sL_array = []

for index_VGB, x in enumerate(VGB1):
    # Calculate the surface potential at the source 
    VGB = x
#    psi_s0 = scipy.optimize.broyden1(Surface_Potential_Source_VS_ExtBias_Numerical,1,f_tol = 1e-14)
#    psi_sL = scipy.optimize.broyden1(Surface_Potential_Drain_VS_ExtBias_Numerical,1,f_tol = 1e-14)
    sol_psi_s0 = scipy.optimize.root(Surface_Potential_Source_VS_ExtBias_Numerical, 1, method = 'lm')
    psi_s0 = sol_psi_s0.x     
    sol_psi_sL = scipy.optimize.root(Surface_Potential_Drain_VS_ExtBias_Numerical, 1, method = 'lm')    
    psi_sL = sol_psi_sL.x
    psi_s0_array.append(psi_s0)
    psi_sL_array.append(psi_sL)
    x0, x1, x2 = Drift_Current_VS_Surface_Potential_Complete(psi_s0,psi_sL)
    Idrift.append(x0)
    Idiff.append(x1)
    IDS.append(x2)

# Plot IDS, Idrift, Idiff vs VGB

fig = plt.figure()
plt.plot(VGB1,np.log(Idrift),linewidth=1.5,marker = 'o',label = '$I_{drift}$')
plt.plot(VGB1,np.log(Idiff),linewidth=1.5,marker = 's',label = '$I_{diff}$')
plt.plot(VGB1,np.log(IDS),linewidth=1.5,label = '$I_{DS}$')
plt.ylim((-20,-10))
plt.grid()
plt.ylabel('Log axis (A)', fontsize=18)
plt.xlabel(r"$V_{GB} (V)$", fontsize=18, fontweight='bold')
#plt.legend(['$Idrift$','$Idiff$','$IDS$'],loc = 'lower right',fontsize=15)
plt.legend(loc = 'lower right',fontsize=15)
plt.gcf().subplots_adjust(left=0.2)
#plt.savefig(Folder +'Comparision of drift and diffusion current of MOSFET for Various VGB.pdf')


# Plot the whole nMOSFET I-V curve
VGB2 = np.arange(0.6+VSB,2.7+VSB,0.3)
VDB2 = np.arange(0+VSB,2.0+VSB,0.1)
Id = np.zeros((len(VGB2),len(VDB2)))


for VGB2_indx,x in enumerate(VGB2):
    for VDB2_indx , y in enumerate(VDB2):
        VGB = x
        VDB = y
        sol_psi_s0 = scipy.optimize.root(Surface_Potential_Source_VS_ExtBias_Numerical, 1, method = 'lm')
        psi_s0 = sol_psi_s0.x     
        sol_psi_sL = scipy.optimize.root(Surface_Potential_Drain_VS_ExtBias_Numerical, 1, method = 'lm')    
        psi_sL = sol_psi_sL.x
        x0, x1, x2 = Drift_Current_VS_Surface_Potential_Complete(psi_s0,psi_sL)
        Id[VGB2_indx][VDB2_indx] = x2



fig = plt.figure()


for m in range(len(VGB2)): 
    plt.plot(VDB2-VSB,Id[m][:],linewidth=2.0,marker= markers[m],label = "$V_{GS}=$ "+ str(VGB2[m])+" V")
#plt.plot(VGB1,np.log(IDS),linewidth=2.0)
plt.grid()
plt.axhline(y=0, color='k')
plt.axvline(x=0, color='k')
plt.ylabel('$I_{DS} (A)$', fontsize=18)
plt.xlabel(r"$V_{DS} (V)$", fontsize=18, fontweight='bold')
plt.legend(bbox_to_anchor=(0, 1), loc='upper left', ncol=1)
plt.gcf().subplots_adjust(left=0.2)
#plt.savefig(Folder +'IDS for MOSFET.pdf')














