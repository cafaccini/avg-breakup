#data from nf-b50
#parameters for the simulation were gamma=0.3
#Tmax=2450 K
#width=1100 um
#======================================================
# This code takes in the simulation results, 
#and plots them.
#
#======================================================

import numpy as np
from matplotlib import pyplot as plt 
from statistics import mean
from scipy import optimize as opt
import parameters as para

print('loading experimental data...')

beta = para.MSC #constant. unitless.
g = para.surface_tension #  Units: J/m^2
C = beta*para.a*para.A/g; #constant for simplifying maths. units: s

epsilon = 0.1  #parameter for T and Tbu comparison.

#defining functions for fitting:

def expo_fit(x,a,b,c):
    y = a + b*np.log(x/c)
    return y

def expo_dec_fit(x,a,b,c):
    y = a*np.exp(-x/b) + c
    return y

def lin_fit(x,a,b):
    y = a + b*x
    return y

def ssd_func(x,y):
    diff2 = (x - y)**2
    return np.sum(diff2)

#=================== data input: Data from simulation ===================

feed = [10, 30, 50, 90, 150] #feed velocity in microns/s
period = [115, 251, 346, 502, 615] #breakup period, in microns
error_l = [18, 10, 18, 30, 15] #error in period

# zbu = np.array([806, 540.7, 275.6, 131.8])
# xbu_sim = np.zeros(len(zbu))
# xbu_sim[:] = (970 - zbu[:] + 510)*(2/0.97)

xbu_sim = [1235, 1495, 1689, 1983, 2285]
err_xbu = [13, 20, 30, 36, 32]
#
#travel_time = [63.7, 45, 20.5, 14.5]
Tbu_sim = [2202, 2254, 2286, 2324, 2354] #Temperature at breakup position, in K.
#
#dT_sim = [0.2069, 0.1576, 0.1155, 0.0962] #derivative of T at xbu

#=======================================================================

#call x, T from parameters:
x = para.x
T = para.T
derivT = para.derivT


#data plots!


#plot breakup period:

#plt.figure()
#plt.rcParams.update({'font.size': 30})
#plt.title('Breakup period as function of feed speed')
#plt.plot(feed, period, 'go')
#plt.xlabel('$v_f$ [$\mu$m/s]')
#plt.ylabel('$\lambda$ [[$\mu$m]')

#plot T(x) and SIMULATION Tbu

fig, ax = plt.subplots()
plt.rcParams.update({'font.size': 30})
plt.title('Temperature Profile, nf_b50')
plt.plot(x, T, label = 'T', lw=3)
plt.plot(xbu_sim, Tbu_sim, 'ro', label = 'data from simulation')
ax.axvspan(xbu_sim[0], xbu_sim[-1], alpha=0.175, color = 'orange', label = 'breakup window')
plt.axhline(y=para.Tsi, lw=3, color='black', linestyle = 'dashed')
plt.axhline(y=para.max_maxT, lw=3, color='black', linestyle = 'dashed')
plt.xlabel('Distance from melting position of Si [$\mu$m]')
plt.ylabel('Temperature [K]')

#plt.figure()
#plt.rcParams.update({'font.size': 30})
#plt.title('Tbu as function of feed speed')
#plt.plot(feed[0:4], Tbu, 'go', label = 'Tbu', lw=3)


#======================================================================

#Now we define continuous variables for the feed speed:

max_vf = para.max_vf #create vf_interp = continuous variable from  1 to maximum vf.
vf_interp = para.vf_interp   #interval between speeds is 1

Tbu_data = np.zeros(len(feed)) 
dT_data = np.zeros(len(feed)) 

dT_data2 = np.zeros(len(feed)) 
dT_data3 = np.zeros(len(feed))
dT_data4 = np.zeros(len(feed))  

Tbu_data2 = np.zeros(len(feed)) 
Tbu_data3 = np.zeros(len(feed))
Tbu_data4 = np.zeros(len(feed))

#l1 = np.zeros(len(Tbu_sim))
#l2 = np.zeros(len(Tbu_sim))
#l3 = np.zeros(len(Tbu_sim))
#l4 = np.zeros(len(Tbu_sim))
#l5 = np.zeros(len(Tbu_sim))
#
#l_t = np.zeros(len(Tbu_sim))



#err_Tbu_data = np.zeros(len(feed)) 
#err_dTdata = np.zeros(len(feed))

#======================================================================
#calculate Tbu from data: Tbu = 

#for i in range(len(feed)): 
#    l = period[i]
#    v = feed[i]
#    Tbu = para.B/(np.log(l/(C*v)))
#    Tbu_data[i] = Tbu #Tbu from data
#    dT_data[i] = (Tbu*Tbu)/(para.B*l) #dT/dx at xbu from data
#    
#for i in range(len(feed)): 
#    l = period[i]
#    v = feed[i]
#    Tbu = para.B/(np.log(l/(C2*v)))
#    Tbu_data2[i] = Tbu #Tbu from data
#    dT_data2[i] = (Tbu*Tbu)/(para.B*l) #dT/dx at xbu from data    
#for i in range(len(feed)): 
#    l = period[i]
#    v = feed[i]
#    Tbu = para.B/(np.log(l/(C3*v)))
#    Tbu_data3[i] = Tbu #Tbu from data
#    dT_data3[i] = (Tbu*Tbu)/(para.B*l) #dT/dx at xbu from data   
#for i in range(len(feed)): 
#    l = period[i]
#    v = feed[i]
#    Tbu = para.B/(np.log(l/(C4*v)))
#    Tbu_data4[i] = Tbu #Tbu from data
#    dT_data4[i] = (Tbu*Tbu)/(para.B*l) #dT/dx at xbu from data      
    
#plot T(x) from simulation, Tbu_sim, and Tbu predicted
#fig, ax = plt.subplots()
#plt.rcParams.update({'font.size': 30})
#plt.title('Temperature Profile, sim.pg69, $\gamma=0.3$')
#plt.xlim(0,5000)
#plt.plot(x, T, label = 'T', lw=3)
#plt.plot(xbu_sim, Tbu_sim, 'ro', label = 'data from simulation', ms=12)
#plt.plot(xbu_sim, Tbu_data[0:4], 'g^', label = 'beta=20', ms=12)
#plt.plot(xbu_sim, Tbu_data2[0:4], 'bs', label = 'beta=40', ms=12)
#plt.plot(xbu_sim, Tbu_data3[0:4], 'mp', label = 'beta=55', ms=12)
#plt.plot(xbu_sim, Tbu_data4[0:4], 'yD', label = 'beta=70', ms=12)
#ax.axvspan(xbu_sim[0], xbu_sim[-1], alpha=0.175, color = 'orange', label = 'breakup window')
#plt.axhline(y=para.Tsi, lw=3, color='black', linestyle = 'dashed')
#plt.axhline(y=para.max_maxT, lw=3, color='black', linestyle = 'dashed')
#plt.xlabel('Distance from melting position of Si [$\mu$m]')
#plt.ylabel('Temperature [K]')



#plot derivative of T from simulation, dT_sim, and model predictions for different beta
#
#fig, ax = plt.subplots()
#plt.rcParams.update({'font.size': 30})
#plt.title('Derivative of Temperature Profile, sim.pg69, $\gamma=0.3$')
#plt.xlim(0,5000)
#plt.plot(x[:-1], derivT, label = 'dT/dx', lw=3)
#plt.plot(xbu_sim, dT_sim, 'ro', label = 'data from simulation', ms=12)
#plt.plot(xbu_sim, dT_data[0:4], 'g^', label = 'beta=20', ms=12)
#plt.plot(xbu_sim, dT_data2[0:4], 'bs', label = 'beta=40', ms=12)
#plt.plot(xbu_sim, dT_data3[0:4], 'mp', label = 'beta=55', ms=12)
#plt.plot(xbu_sim, dT_data4[0:4], 'yD', label = 'beta=70', ms=12)
#plt.xlabel('Distance from melting position of Si [$\mu$m]')
#plt.ylabel('Temperature Gradient [K/m]')



#======================================================================
#Let's interpolate Tbu (from the simulation data!)

#Interpolate Tbu_data and dT_data over a continuous vf:

#initial_par = [100, 50, 1e-8]   #initial parameter guesses for the fit (from origin)  
#popt, pcov = opt.curve_fit(expo_fit, feed, Tbu_data, p0 = initial_par) #, sigma = err_Tbu_data)
#Tbu_data_interp = expo_fit(vf_interp, *popt)
#sig_Tbu_data = np.sqrt(np.diag(pcov)) #stdev. for Tbu_data_interp
#
##plt.figure()
##plt.plot(feed, Tbu_data, 'go', label = 'Simulation data')
##plt.plot(vf_interp, Tbu_data_interp, lw=3, label = 'Exponential fit')
##plt.title('Interpolation of Tbu')
##plt.xlabel('vf [$\mu$m/s]')
##plt.ylabel('Tbu [K]')
#
#initial_par2 = [1, 20, 1e-2]   #initial parameter guesses for the fit (from origin)  
#popt2, pcov2 = opt.curve_fit(expo_fit, feed, dT_data, p0 = initial_par2)#, sigma = err_dTdata)
#dT_interp = expo_fit(vf_interp, *popt2) #fitting from discete values.
#sig_dT = np.sqrt(np.diag(pcov2))
#
#initial_par3 = [1, 20, 1e-2]   #initial parameter guesses for the fit (from origin)  
#popt3, pcov3 = opt.curve_fit(expo_fit, feed, period, p0 = initial_par3)#, sigma = err_dTdata)
#l_interp = expo_fit(vf_interp, *popt3) #fitting from discete values.
#sig_l = np.sqrt(np.diag(pcov3))

#======================================================================
#For each beta (20, 40, 55, 60), calculate l(vf) using tbu_sim, period
#
#for i in range(len(Tbu_sim)):
#	Tbu = Tbu_sim[i]
#	v = feed[i]
#	l1[i] = v*C*np.exp(para.B/Tbu)
#
#for i in range(len(Tbu_sim)):
#	Tbu = Tbu_sim[i]
#	v = feed[i]
#	l2[i] = v*C2*np.exp(para.B/Tbu)
#	
#for i in range(len(Tbu_sim)):
#	Tbu = Tbu_sim[i]
#	v = feed[i]
#	l3[i] = v*C3*np.exp(para.B/Tbu)
#	
#for i in range(len(Tbu_sim)):
#	Tbu = Tbu_sim[i]
#	v = feed[i]
#	l4[i] = v*C4*np.exp(para.B/Tbu)	
#	
## Also calculate l = tbu^2/(B*dt/dx)
#
#
#for i in range(len(Tbu_sim)):
#	Tbu = Tbu_sim[i]
#	v = feed[i]
#	dT = dT_sim[i]
#	l_t[i] = (Tbu*Tbu)/(para.B*dT)

#plot derivative of T from simulation, dT_sim, and model predictions for different beta
#
#fig, ax = plt.subplots()
#plt.rcParams.update({'font.size': 30})
#plt.title('Breakup Period, sim.pg69, $\gamma=0.3$')
#plt.xlim(0,5000)
##plt.plot(x[:-1], derivT, label = 'dT/dX', lw=3)
#plt.plot(xbu_sim, l_t, 'c*', label = 'modeled', ms=12 )
#plt.plot(xbu_sim, period[0:4], 'ro', label = 'data from simulation', ms=12)
#plt.plot(xbu_sim, l1, 'g^', label = 'beta=20', ms=12)
#plt.plot(xbu_sim, l2, 'bs', label = 'beta=40', ms=12)
#plt.plot(xbu_sim, l3, 'mp', label = 'beta=55', ms=12)
#plt.plot(xbu_sim, l4, 'yD', label = 'beta=70', ms=12)
#plt.xlabel('Distance from melting position of Si [$\mu$m]')
#plt.ylabel('$\lambda$ [$\mu$m]')
#
#
#plt.figure()
#plt.rcParams.update({'font.size': 30})
#plt.title('Breakup period, sim.pg69, $\gamma=0.3$')
#plt.plot(feed, period, 'ro', label = 'data from simulation', ms=12 )
#plt.plot(feed[0:4], l1, 'g^', label = 'beta=20', ms=12 )
#plt.plot(feed[0:4], l2, 'bs', label = 'beta=40', ms=12 )
#plt.plot(feed[0:4], l3, 'mp', label = 'beta=55', ms=12)
#plt.plot(feed[0:4], l4, 'yD', label = 'beta=70', ms=12)
#plt.xlabel('$v_f$ [$\mu$m/s]')
#plt.ylabel('$\lambda$ [$\mu$m]')


print('...done')





