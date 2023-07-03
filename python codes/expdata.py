#============================================================================
#Dataset from Alexander's experiment
#Fiber: Si-core, 4microns
#
#============================================================================
import numpy as np
from matplotlib import pyplot as plt 
from statistics import mean
from scipy import optimize as opt
import parameters as para
from scipy.interpolate import interp1d

print('loading experimental data...')

x = para.x
beta = para.MSC #constant. unitless.
g = para.surface_tension #  Units: J/m^2
C = beta*para.a*para.A/g; #constant for simplifying maths. units: s
epsilon = 0.1  #parameter for T and Tbu comparison.


#defining functions for fitting:

def expo_fit(x,a,b,c):
    y = a + b*np.log(x/c)
    return y

def expo_fit2(x,a,b,c):
    y = a - b*np.exp(-x/c)
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

#=================== data input: Experimental data ===================

feed = [2,10,30,50,90,150] #feed velocity in microns/s
period = [280,430,530,630,715,795] #breakup period, in microns
travel= [1003, 1206, 1426, 1638, 1802, 1996.5 ]
traveltime = [25.37, 18.18, 10.8, 7.35, 3.65, 3.24] #derivative of travel distance. in s.

sigL = [6, 9, 17, 7, 13, 7] #standard deviation in the period, measured. 
sig_trav = [95, 95, 95, 95, 95, 95]
#======================================================================
#Now we get all the information we can from here:
#note: *_data indicates variable is directly calculated from experimental data
    
Tbu_data = np.zeros(len(feed)) 
dT_data = np.zeros(len(feed)) 
err_Tbu_data = np.zeros(len(feed)) 
err_dTdata = np.zeros(len(feed))

max_vf = feed[-1] #create vf_interp = continuous variable from  1 to maximum vf.
vf_interp = np.linspace(1,max_vf,max_vf)   #interval between speeds is 1


#calculate Tbu from data: Tbu = B/ln(l/Cvf), C = baA/g

for i in range(len(feed)): 
    l = period[i]
    v = feed[i]
    Tbu = para.B/(np.log(l/(C*v)))
    Tbu_data[i] = Tbu #Tbu from data
    aux = (C*v*sigL[i])/(l*np.log(l/(C*v)))
    err_Tbu_data[i] = np.abs(Tbu)*np.sqrt(((para.sigB/para.B)*(para.sigB/para.B))+(aux*aux))
    dT_data[i] = (Tbu*Tbu)/(para.B*l) #dT/dx at xbu from data
    f = Tbu*Tbu
    sigf = f*np.sqrt(2*(err_Tbu_data[i]/Tbu)*(err_Tbu_data[i]/Tbu))
    g = para.B*l
    sigg = g*np.sqrt(((para.sigB/para.B)*(para.sigB/para.B))+((sigL[i]/l)*(sigL[i]/l)))
    err_dTdata[i] = (f/g)*np.sqrt((sigf/f)*(sigf/f)+(sigg/g)*(sigg/g))
    
#======================================================================
#Fit Tbu_data and dT_data over a continuous vf:

initial_par = [100, 50, 1e-8]   #initial parameter guesses for the fit (from origin)  
popt, pcov = opt.curve_fit(expo_fit, feed, Tbu_data, p0 = initial_par, sigma = err_Tbu_data)
Tbu_data_interp = expo_fit(vf_interp, *popt)
sig_Tbu_data = np.sqrt(np.diag(pcov)) #stdev. for Tbu_data_interp

if __name__ == '__main__':
    plt.figure()
    plt.rcParams.update({'font.size': 20})
    plt.title('Calculated Breakup Temperature , for beta=%i' %beta)
    #plt.ylim(1950, 2300)
    plt.plot(feed, Tbu_data, 'go', ls = 'none', label = 'Calculated $T_{bu}$ for each experiemtal data point', ms = 12)
    plt.errorbar(feed, Tbu_data, yerr = err_Tbu_data, xerr = None, ecolor = 'g', ls = 'none', capsize =5, elinewidth = 5)
    plt.plot(vf_interp, Tbu_data_interp, 'b--', lw = 5, label = 'fit of data points')
    #plt.fill_between(vf_interp, Tbu_data_interp - sig_Tbu_data[1], Tbu_data_interp + sig_Tbu_data[1], alpha=0.2)
    plt.axhline(y=para.Tsoft, xmin=0, xmax = max(feed), color = 'orange', lw = 3, label = 'Melting temperature of Si' )
    plt.axhline(y=para.max_maxT, xmin=0, xmax = max(feed), color = 'r', lw = 3, label = 'Boiling temperature of Silica' )
    plt.xlabel('Feed velocity ($\mu m/s$)')
    plt.ylabel('Breakup Temperature (K)')

#======================================================================

l_interp = vf_interp*C*np.exp(para.B/Tbu_data_interp)

if __name__ == '__main__':
    plt.figure()
    plt.rcParams.update({'font.size': 20})
    plt.title('Breakup period from Experimental data')
    plt.plot(feed, period, 'go', label = 'Experimental data')
    plt.errorbar(feed, period, yerr = sigL, xerr = None, ecolor = 'g', ls = 'none', capsize =5, elinewidth = 3)
    plt.plot(vf_interp, l_interp, 'b--', label = 'from fitted Tbu', lw = 3)
    plt.xlabel('Feed velocity ($\mu m/s$)')
    plt.ylabel('Breakup periood ($\mu m$)')


#======================================================================


initial_par3 = [1, 20, 1e-2]   #initial parameter guesses for the fit (from origin)  
popt3, pcov3 = opt.curve_fit(expo_fit, feed, dT_data, p0 = initial_par3, sigma = err_dTdata)
dT_interp = expo_fit(vf_interp, *popt3) #fitting from discete values.
sig_dT = np.sqrt(np.diag(pcov3))


if __name__ == '__main__':
    dT_interp2 = np.zeros(len(Tbu_data_interp))
    for i in range(len(Tbu_data_interp)):
     	dT_interp2[i] = (Tbu_data_interp[i]**2)/(para.B*l_interp[i])
    
    plt.figure()
    plt.rcParams.update({'font.size': 30})
    plt.title('dT/dx at xbu from Experimental data')
    plt.plot(feed, dT_data, 'go', label = 'Experimental data')
    plt.errorbar(feed, dT_data, yerr = err_dTdata, xerr = None, ecolor = 'g', ls = 'none', capsize =5, elinewidth = 3)
    plt.plot(vf_interp, dT_interp, 'b--', label = 'Interpolated', lw = 3)
    plt.fill_between(vf_interp, dT_interp - sig_dT[1], dT_interp + sig_dT[1], alpha=0.2)
    plt.plot(vf_interp, dT_interp2, 'r--', label = 'Interpolated2', lw = 3)
    plt.xlabel('Feed velocity ($\mu m/s$)')
    plt.ylabel('dT/dx at xbu ($K/\mu m$)')


print('...done')





