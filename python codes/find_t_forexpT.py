# AUX CODE for temperature profile determination.

import numpy as np 
from matplotlib import pyplot as plt 
from scipy import optimize as opt
import parameters as para
import data as data

beta = para.MSC #constant. might be changed later.
g = para.surface_tension # migth be changed later. Units: J/m^2
C = beta*para.a*para.A/g; #constant for simplifying maths. units: s
epsilon = 1 #parameter for T and Tbu comparison.

#defining exponential function for fitting:
def expo_fit(x,a,b,c):
    y = a + b*np.log(x/c)
    return y

#function: sum of squared difference
def ssd_func(x,y):
    diff2 = (x - y)**2
    return np.sum(diff2)

x = para.x

#load experimental data: values interpolated for vf from 0 to 150 microns/s
Tbu_data = data.Tbu_data_interp
vf = data.vf_interp
l_data = data.l_interp

dT_data_interp = data.dT_interp
sig_dT = data.sig_dT


xbu = np.zeros(len(Tbu_data))
b = np.zeros(len(Tbu_data))
#============================================================================
#Up to here nothing depends on the Temperature profile.
Tmax0 = 2550 #in K. Maximal temperature of proposed T profile
width0 = 2000 #guess value! units: microns

initial = (Tmax0, width0)

fixed_para =[Tbu_data,dT_data_interp]

def profile(tw,*fixed_para):  #tw = (Tmax, width)
    Tbu = fixed_para[0]
         
    #delta_T = (tw[0] - para.Troom)/2    
    #shift = np.arctanh((2*para.Tsi-tw[0]-para.Troom)/(tw[0]-para.Troom))
    
    T = tw[0] - (tw[0]-para.Tsi)*np.exp(-x/tw[1]) #Generate T profile.   
	    
    derivT = np.diff(T) #calculate dT/dx
        
    for i in range(len(Tbu)):
	    var = 10
	    dist = 1
	    while var > epsilon and dist < para.max_x:
		    var = np.abs(T[dist]-Tbu[i])
		    dist = dist + 1
	    xbu[i] = dist
	    b[i] = derivT[dist]   # obtain dT/dx at the xbu values:  b(i)
	            
    return T, b

def profile_optimizer(tw,fixed_para):
	print(tw)
	T, b = profile(tw,*fixed_para)
	
	return ssd_func(b,fixed_para[1])
    
bnds = ((max(Tbu_data),3000), 
	(1000,para.max_x/2))    

min_results = opt.minimize(profile_optimizer,
			   initial,
			   args = fixed_para,
			   method = 'L-BFGS-B',
			   bounds = bnds)


Tmax, width = min_results.x 


#delta_T = (Tmax - para.Troom)/2
#shift = np.arctanh((2*para.Tsi-Tmax-para.Troom)/(Tmax-para.Troom))

T = Tmax - (Tmax-para.Tsi)*np.exp(-x/width) #Generate T profile. 
derivT = np.diff(T)

for i in range(len(Tbu_data)):
	var = 10
	dist = 1
	while var > epsilon and dist < para.max_x:
	    var = np.abs(T[dist]-Tbu_data[i])
	    dist = dist + 1
	xbu[i] = dist
	b[i] = T[dist]   # obtain T at the xbu values:  b(i)


plt.figure()
plt.rcParams.update({'font.size': 20})
plt.plot(vf, Tbu_data, 'b-', label = 'from experimental data', lw = 3)
plt.plot(vf, b, 'r-', label = 'Tbu at $x_{bu}$', lw = 3)
plt.xlabel('Feed velocity ($\mu m$)')
plt.ylabel('T (K)')
plt.title('proposed T profile at the breakup position $x_{bu}$, as function of the feed speed')




























