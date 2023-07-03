#data from nf_b50
#parameters for the simulation were
#gamma=0.3
#Tmax=2450 K
#width=1100 um
#maxz = 970 su

# -*- coding: utf-8 -*-
#Parameters files. contains constants and other parameters that will 
# NOT change during the model code.

#======================================================
# This code takes in the simulation parameters and problem constants, 
#and provides the viscosities and temperature profile.
#no model calculations are made!
#======================================================

import numpy as np
from matplotlib import pyplot as plt 

print('loading model parameters...')

MSC = 60 #Marginal Stability Criterion constant

#======= Fiber parameters: ============
a = 2e-6 #core radius in m

#A and B are viscosity parameters for Silicon
A = 5.7e-8 #in Pa*s
B = 61812 #in K. 
sigB = 52 #uncertainty in B

surface_tension = 1.2 #units: J/m^2


#======== Temperature profile parameters: depends on materials ============
Troom = 293.15 #room temperature, in K

Tsi = 1687.15 #in K. Melting temperature of Si
Tsoft = 1973.15 #in K. softening temperature of silica


min_maxT=Tsoft #in K. softening T of silica
max_maxT=2503.15 #in K. boiling T of silica

max_x = 10000 #in microns. distance from Si melting point

x = np.linspace(0,max_x+2,max_x+2) 

max_vf = 150 #create vf_interp = continuous variable from  1 to maximum vf.
vf_interp = np.linspace(0,max_vf,max_vf)   #interval between speeds is 1


#======================================================
#T profile AS DEFINED IN THE SIMULATION:

Tmax = 2440 #in K. Maximal temperature of proposed T profile
width = 1100 # in microns

#======================================================


T = Tmax - (Tmax-Tsi)*np.exp(-x/width) #proposed T profile
derivT = np.diff(T) #calculate dT/dx

visco_silica = A*np.exp(B/T) #silica viscosity

visco_si = np.zeros(len(T))

for i, taux in enumerate(T):
	t = taux
	if taux>1900:
		t = 1900
	
	visco_si[i] = 0.003632339 - (2.86*1E-6)*t + (6.35*1E-10)*(t**2)#viscosity of silicon
	

visco_ratio = visco_si/visco_silica
#plt.figure()
#plt.rcParams.update({'font.size': 30})
#plt.title('Viscosity Ratio as function of the temperature')
#plt.semilogy(T, visco_ratio, label = 'viscosity ratio')
#plt.xlabel('Temperature [K]')
#plt.ylabel('$\eta_{Si}/\eta_{SiO_2}$')


soft = [Tsoft] * len(vf_interp)


#plt.figure()
#plt.rcParams.update({'font.size': 30})
#plt.semilogy(T, visco_silica, label = 'SiO2 viscosity')
#plt.semilogy(T, visco_si, label = 'Si viscosity')
#plt.xlabel('Temperature [K]')
#plt.ylabel('Viscosity [Pa*s]')

#plt.figure()
#plt.rcParams.update({'font.size': 30})
#plt.title('Temperature Profile')
#plt.plot(x, T, label = 'T')
#plt.xlabel('position [$\mu$m]')
#plt.ylabel('Temperature [K]')

print('...done')
