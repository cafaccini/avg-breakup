# MAIN CODE for model validation.
#The caldulations here DO NOT depend on the simulation data, which are plotted for reference ONLY.

import numpy as np
from matplotlib import pyplot as plt 
from scipy import optimize as opt
import parameters as para
import data as data
import expdata as expdata

def lin_fit(x,a,b):
    y = a + b*x
    return y

# The experimental data is loaded from "data"
# the model parameters are loaded from "parameters"
beta = para.MSC
g=para.surface_tension
x = para.x
vf = data.vf_interp
C = beta*para.a*para.A/g; #data.C #change if different beta value is desired!!!
epsilon = data.epsilon

#load experimental values:
Tbu_dat = data.Tbu_data
#sig_Tbu_dat = data.sig_Tbu_data
dT_dat = data.dT_data
#sig_dT_dat = data.sig_dT

#============================================================================
#define desired T(x):

T=para.T
derivT = np.diff(T) #calculate dT/dx

visco_silica = para.A*np.exp(para.B/T) #silica viscosity
soft = [para.Tsoft] * len(vf)

#============================================================================
#Given the T(x) -> find xbu!
#use the capillary number to determine xbu! at xbu, Ca = 1, dCa/dx = 0

Tb = T[:-1] #to match size of dT array
g = para.B*derivT/(Tb**2) #g=g(x)

#For quick plots:
#f2 = np.exp(-para.B/Tb)/(vf[2]*C) 
#f10 = np.exp(-para.B/Tb)/(vf[10]*C) 
#f50 = np.exp(-para.B/Tb)/(vf[50]*C)
#f100 = np.exp(-para.B/Tb)/(vf[100]*C)  
#
#plt.figure()
#plt.rcParams.update({'font.size': 30})
#plt.plot(x[:-1],g, label = 'g(x)')
#plt.plot(x[:-1],f2, label = 'f(x) for vf =2')
#plt.plot(x[:-1],f10, label = 'f(x) for vf =10')
#plt.plot(x[:-1],f50, label = 'f(x) for vf =50')
#plt.plot(x[:-1],f100, label = 'f(x) for vf =100')
#plt.xlabel('x [$\mu m$]')
#============================================================

f = np.zeros((len(Tb),len(vf)))
fgdiff = np.zeros((len(Tb),len(vf)))

for i in range(len(vf)):
	f[:,i] = np.exp(-para.B/Tb)/(vf[i]*C)  #f = f(x,vf)
	fgdiff[:,i] = np.abs(f[:,i]-g[:])

fgzeros = fgdiff < 1e-5   
idxs_fgzeros = np.where(fgzeros)
idxs_fgzeros = list(zip(idxs_fgzeros[0],idxs_fgzeros[1]))
#check to see if all columns of fgzeros have AT LEAST one true value
all_have_0 = np.zeros(len(vf), dtype=bool)

for i in range(len(vf)):
	all_have_0[i] = np.any(fgzeros[:,i])
print('all columns have at least one zero:', np.all(all_have_0))

#now we find xbu for each vf:
#declare variables:
mean_xbu = np.zeros(len(vf))
std_xbu = np.zeros(len(vf))
mean_Tbu = np.zeros(len(vf))
std_Tbu = np.zeros(len(vf))
mean_dT = np.zeros(len(vf))
std_dT = np.zeros(len(vf))

range_of_xbu = np.zeros(len(vf))
range_of_xbu[:] = False
range_of_xbu = list(range_of_xbu)

range_of_T = [False for i in range(len(vf))]
range_of_dT = [False for i in range(len(vf))]

for coor in idxs_fgzeros:
	if not range_of_xbu[coor[1]]:
		range_of_xbu[coor[1]] =[x[coor[0]]]
		range_of_T[coor[1]] =[T[coor[0]]]
		if coor[0] != 10001:
			range_of_dT[coor[1]] =[derivT[coor[0]]]
	else:
		range_of_xbu[coor[1]].append(x[coor[0]])
		range_of_T[coor[1]].append(T[coor[0]])
		if coor[0] != 10001:
			range_of_dT[coor[1]].append(derivT[coor[0]])

for i, (lx, lT, ldT) in enumerate(zip(range_of_xbu,range_of_T, range_of_dT)):
	mean_xbu[i] = np.mean(lx)
	std_xbu[i] = np.std(lx)
	
	mean_Tbu[i] = np.mean(lT)
	std_Tbu[i] = np.std(lT)
	mean_dT[i] = np.mean(ldT)
	std_dT[i] = np.std(ldT)


#plt.figure()
#plt.rcParams.update({'font.size': 30})
#plt.title('Breakup position')
#plt.plot(vf, mean_xbu, 'r-', label = 'Predicted values', lw=3)
#plt.vlines(2, min(mean_xbu), max(mean_xbu))
#plt.vlines(10, min(mean_xbu), max(mean_xbu))
#plt.vlines(30, min(mean_xbu), max(mean_xbu))
#plt.vlines(50, min(mean_xbu), max(mean_xbu))
#plt.vlines(90, min(mean_xbu), max(mean_xbu))
#plt.vlines(150, min(mean_xbu), max(mean_xbu))
#plt.xlabel('Feed velocity [$\mu m/s$]')
#plt.ylabel('$x_{bu}$ $[\mu m]$') 

#============================================================================
#Get Ca, dCa/dx and plot  
Ca = np.zeros(len(x))
dCa = np.zeros(len(x))

#============================================================================
#Predict breakup period: l_model

l_model = np.zeros(len(vf))
l_2 = np.zeros(len(vf))
std_l = np.zeros(len(vf))
std_l2 = np.zeros(len(vf))

for i in range(len(vf)):
	l_model[i] = C*vf[i]*np.exp(para.B/mean_Tbu[i]) 
	f = para.B/mean_Tbu[i]
	sigf = np.abs(f)*np.sqrt(((std_Tbu[1]/mean_Tbu[i])**2)+((para.sigB/para.B)**2))
	std_l[i] = l_model[i]*sigf
	l_2[i] = (mean_Tbu[i]*mean_Tbu[i])/(para.B*mean_dT[i]) 
	g = mean_Tbu[i]**2
	sigg = g*np.sqrt(2*((std_Tbu[1]/mean_Tbu[i])**2))
	h = para.B * mean_dT[i]
	sigh = np.abs(h)*np.sqrt(((para.sigB/para.B)**2)+((std_dT[i]/mean_dT[i])**2))
	std_l2[i] = np.abs(g/h)*np.sqrt(((sigg/g)**2)+((sigh/h)**2))
	
#note: l_model and l_2 MUST be the same. They are just two different ways of obtaining l
# (from 1st and 2nd equation of the model, respectively.)

#============================================================================
#Calculate travel time = dxbu/dvf

travel_time = np.diff(mean_xbu)  


# logt = np.log(travel_time)
# logv = np.log(vf[1:])
# fit_par = [1, 10]   #initial parameter guesses for the fit (from origin)  
# popt3, pcov3 = opt.curve_fit(lin_fit, logv, logt, p0 = fit_par)
# tt_fit = lin_fit(logv, *popt3)
# tt_fit = np.exp(tt_fit)

# err_tt = np.zeros(len(travel_time))  #travel time error: have to improve this!
# for i in range(len(travel_time)):
#  	err_tt[i] = np.abs(tt_fit[i]-travel_time[i])

# err_tt = np.mean(err_tt)


#============================================================================
# @@@@@@@@@@@@@@@@@@@@@@@
# COLOR HERE

color = "#e6e8e8"
#============================================================================

plt.figure(figsize=(6,5))
plt.rcParams.update({'font.size': 18})
# plt.title('Simulation vs model, $\gamma$=0.3 N/m, $\beta$= %i' %beta)
plt.xlim(0,max(vf)+10)
# plt.plot(expdata.feed, expdata.period, 'go', label = 'experiment')
# plt.errorbar(expdata.feed, expdata.period, yerr = expdata.sigL, xerr = None, ecolor = 'g', ls = 'none', capsize =5, elinewidth = 4)
plt.plot(data.feed, data.period, 'ro', label = 'simulation')
plt.errorbar(data.feed, data.period, yerr = data.error_l, xerr = None, ecolor = 'r', ls = 'none', capsize =5, elinewidth = 4)
plt.plot(vf, l_model, 'b-', label = 'model using beta=%i'%beta, lw = 4)
# plt.fill_between(vf, l_model - std_l, l_model + std_l, alpha=0.175, color = 'blue', label = '1 $\sigma$')
# plt.fill_between(vf, l_model - 2*std_l, l_model + 2*std_l, alpha=0.1, color = 'blue', label = '2 $\sigma$')
# plt.plot(data.feed, data.lsim_calc, 'c*')
plt.xlim(0,160)
plt.ylim(30,900)
#plt.xlabel(' $v_f$ [$\mu m/s$]')
#plt.ylabel('$\lambda_{bu}$ $[\mu m]$')
#plt.text(10, 700, 'width=%i um'%width) 
# plt.legend()

ax = plt.gca()
ax.set_facecolor(color)

plt.tight_layout()
plt.show()   

plt.figure(figsize=(6,5))
plt.rcParams.update({'font.size': 18})
# plt.title('Breakup position: nf-b20 vs model, gamma=0.3 N/m, beta= %i' %beta)
# plt.plot(expdata.feed, expdata.travel, 'go', label = 'experiment')
# plt.errorbar(expdata.feed, expdata.travel, yerr = expdata.sig_trav, xerr = None, ecolor = 'g', ls = 'none', capsize =5, elinewidth = 4)
plt.plot(data.feed, data.xbu_sim, 'ro', label = 'simulation')
plt.errorbar(data.feed, data.xbu_sim, yerr = data.err_xbu, xerr = None, ecolor = 'r', ls = 'none', capsize =5, elinewidth = 4)
plt.plot(vf, mean_xbu, 'b-', label = 'model using beta=%i'%beta, lw = 4)
#plt.text(10, 750, 'Tmax=%i K'%Tmax)
#plt.text(10, 700, 'width=%i um'%width)
plt.xlim(0,160)
plt.ylim(650,2700)
#plt.xlabel('$v_f$ [$\mu m/s$]')
#plt.ylabel('$x_{bu}$ $[\mu m]$') 
# plt.legend()

ax = plt.gca()
ax.set_facecolor(color)

plt.tight_layout()
plt.show()
   

# plt.figure()
# plt.rcParams.update({'font.size': 30})
# plt.title('Breakup temperature: nf-b20 vs model, gamma=0.3 N/m, beta= %i' %beta)
# plt.xlim(0,max(vf)+10)
# plt.ylim(para.Tsi, max(mean_Tbu)+100)
# plt.plot(data.feed, data.Tbu_sim, 'ro', label = 'simulation values')
# #plt.errorbar(data.feed, data.period, yerr = None, xerr = None, ecolor = 'g', ls = 'none', capsize =5, elinewidth = 3)
# plt.plot(vf, mean_Tbu, 'b-', label = 'predicted values', lw = 2)
# #plt.fill_between(vf, l_model - std_l, l_model + std_l, alpha=0.175, color = 'blue', label = '1 $\sigma$')
# #plt.fill_between(vf, l_model - 2*std_l, l_model + 2*std_l, alpha=0.1, color = 'blue', label = '2 $\sigma$')
# plt.xlabel(' vf [$\mu m/s$]')
# plt.ylabel('$T_bu$ [K]')

# plt.figure()
# plt.rcParams.update({'font.size': 20})
# plt.title('Travel Time: nf-b20 vs model, gamma=0.3 N/m, beta= %i' %beta)
# plt.xlim(0,max(vf)+10)
# plt.ylim(0,50)
# plt.plot(data.feed, data.traveltime, 'ro', label = 'simulation values')
# plt.plot(expdata.feed, expdata.traveltime, 'go', label = 'experimetal values')
# plt.plot(vf[1:], travel_time, 'b-')
# # plt.plot(vf[1:], tt_fit, 'b-', label = 'predicted values', lw=2)
# # plt.fill_between(vf[1:], tt_fit - err_tt, tt_fit + err_tt, alpha=0.175, color = 'blue', label= '1 $\sigma$')
# # plt.fill_between(vf[1:], tt_fit - 2*err_tt, tt_fit + 2*err_tt, alpha=0.1,  color = 'blue', label= '2 $\sigma$')
# plt.xlabel('Feed velocity [$\mu m/s$]')
# plt.ylabel('travel time [s]')  


