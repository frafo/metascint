# -*- coding: utf-8 -*-
"""
Created on Tue Jan 18 10:50:49 2022

@author: yocon
"""
#Energy spectrum correction

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import scipy
from scipy.stats import norm
from scipy import optimize
import codecs
import pandas as pd
import os


def pp_fit(ev_photopeak,bins, caption):
    CTRhisto,bins = np.histogram(Etotal, bins)
    n = np.size(CTRhisto)                        #the number of data
    mean = np.argmax(CTRhisto)             #note this correction
    sigma = 120      #note this correction
    y = np.array(CTRhisto/np.max(CTRhisto))
    x = np.arange(0,np.size(y)*(bins[2]-bins[1]),(bins[2]-bins[1]))
    popt,pcov = curve_fit(gaus,x,y,p0=[1,mean,sigma])
    mx= np.argmax(gaus(x,*popt))*(bins[2]-bins[1])
    residuals = y- gaus(x, *popt)
    ss_res = np.sum(residuals**2)
    #  You can get the total sum of squares (ss_tot) with
    ss_tot = np.sum((y-np.mean(y))**2)
    # And finally, the r_squared-value with,
    r_squared = 1 - (ss_res / ss_tot)
    plt.plot(x-mx,gaus(x,*popt),label="FWHM: " + str(round(abs(2.355*popt[2]),1))+ " ps, R^2: "+str(round(r_squared,3)))
    return round(abs(2.355*popt[2]),1)


def valley_value(ff):
    Etotal= event_registry_max_real[:,4]-ff*event_registry_max_real[:,2]
    bin_num=80
    Ehisto,bins = np.histogram(Etotal, bins=bin_num)
    bin_size= (bins[1]-bins[0])
    Ehisto_norm =Ehisto/np.max(Ehisto)
    Loc_max = np.argmax(Ehisto_norm)
    Loc_min = int(np.floor(Loc_max*0.7))
    Loc_val = Loc_min+np.argmin(Ehisto_norm[Loc_min:Loc_max])
    Loc_min_real= int(np.ceil(Loc_min*bin_size))
    Loc_val_real= int(np.ceil(Loc_val*bin_size))
    Loc_max_real= int(np.ceil(Loc_max*bin_size))
    return np.min(Ehisto_norm[Loc_min:Loc_max])

def FWHM(ff):
    Ehisto_norm = event_registry_max_real[:,4]-ff*event_registry_max_real[:,2]
    Loc_max = np.argmax(Ehisto_norm) # corresponds to 511 keV    
    return np.argmin(abs(Ehisto_norm[Loc_max:int(np.ceil(Loc_max*1.3))]-0.5))/Loc_max

def FWHM_func(ff):
    Etotal= event_registry_max_real[:,4]-ff*event_registry_max_real[:,2]
    bin_num=80
    Ehisto,bins = np.histogram(Etotal, bins=bin_num)
    bin_size= (bins[1]-bins[0])
    Ehisto_norm =Ehisto/np.max(Ehisto)
    Loc_max = np.argmax(Ehisto_norm)
    Loc_min = int(np.floor(Loc_max*0.7))
    Loc_val = Loc_min+np.argmin(Ehisto_norm[Loc_min:Loc_max])
    Loc_min_real= int(np.ceil(Loc_min*bin_size))
    Loc_val_real= int(np.ceil(Loc_val*bin_size))
    Loc_max_real= int(np.ceil(Loc_max*bin_size))
    Loc_FWHM_real = int((np.argmin(abs(Ehisto_norm[Loc_max:int(np.ceil(Loc_max*1.3))]-0.5)))*bin_size)
    return Loc_FWHM_real*2/Loc_max_real

def gaus(x,a,x0,sigma):
    return a*np.exp(-(x-x0)**2/(2*sigma**2))


def FWHM_fit(ff):
    Etotal= event_registry_max_real[:,4]-ff*event_registry_max_real[:,2]
    bin_num=80
    Ehisto,bins = np.histogram(Etotal, bins=bin_num)
    bin_size= (bins[1]-bins[0])
    sigma= int(bin_num/20)

    Ehisto_norm =Ehisto/np.max(Ehisto)
    Loc_max = np.argmax(Ehisto_norm)
    
    popt, pcov = scipy.optimize.curve_fit(gaus, np.linspace(Loc_max-4*sigma,Loc_max+4*sigma,1), Ehisto_norm[Loc_max-2*sigma: Loc_max+2*sigma], p0=[1.,Loc_max,sigma])
    if popt[2]>0:
        output= 2.355*popt[2]*bin_size/Loc_max
    else:
        output =150*bin_size/Loc_max
    return popt

# step=1000
# ff_vec=np.zeros(step)
# for i in range(step):
#     ff= -5+0.01*i
#     ff_vec[i] = FWHM_fit(ff)

# plt.plot(np.linspace(-5,-5+0.01*step,step), ff_vec, label= str(-5+0.01*np.argmin(ff_vec)))

Optimization = optimize.minimize(FWHM, x0= np.array([0.25]), bounds=[(-2,0.75)], method='Nelder-Mead')

ff=0.25 # Optimization.x
Etotal= event_registry_max_real[:,4]-ff*event_registry_max_real[:,2]
bin_num=60
Ehisto,bins = np.histogram(Etotal, bins=bin_num)
bin_size= (bins[1]-bins[0])

Ehisto_norm =Ehisto/np.max(Ehisto)
Loc_max = np.argmax(Ehisto_norm)
Loc_min = int(np.floor(Loc_max*0.7))
Loc_val = Loc_min+np.argmin(Ehisto_norm[Loc_min:Loc_max])
Loc_min_real= int(np.ceil(Loc_min*bin_size))
Loc_val_real= int(np.ceil(Loc_val*bin_size))
Loc_max_real= int(np.ceil(Loc_max*bin_size))
Loc_FWHM_real = int((Loc_max + np.argmin(abs(Ehisto_norm[Loc_max:int(np.ceil(Loc_max*1.3))]-0.5)))*bin_size)

plt.vlines(Loc_min_real,0,1,'r','dashed')
plt.vlines(Loc_val_real,0,1,'r','dashed')
plt.vlines(Loc_max_real,0,1,'r','dashed')
plt.vlines(Loc_FWHM_real,0,1,'r','dashed')
plt.hlines(0.5,0,bin_size*bin_num,'b','dashed')
valley = np.min(Ehisto_norm[Loc_min:Loc_max])
plt.hlines(valley,0,bin_size*bin_num,'r','dashed')
FWHM= np.argmin(Ehisto_norm[Loc_max:int(np.ceil(Loc_max*1.3))]-0.5) 
# plt.vlines(FWHM,0,1,'dashed')
popt =FWHM_fit(ff)
x_vec=np.arange(Loc_max_real-4*popt[2],Loc_max_real+4*popt[2],1)/bin_size
plt.plot(bins[0:bin_num],Ehisto_norm,label=str(ff) + ' ' +str(Loc_max_real) +' '+str(Loc_val_real)+" "+ str(Loc_min_real)+' '+str(round((Loc_FWHM_real-Loc_max_real)*2/Loc_max_real,2)))  
plt.plot(x_vec, gaus(x_vec,popt[0],popt[1],popt[2]))

# CTRhisto,bins = np.histogram(EtotalPH, bins=30)
# plt.plot(bins[0:np.size(CTRhisto)],CTRhisto/np.max(CTRhisto),label="photopeak")  
# plt.xlabel("Etotal")
# plt.ylabel('#')
plt.legend()
plt.show()    
# CTRhisto,bins = np.histogram(Esharing, bins=50)
# plt.plot(bins[0:np.size(CTRhisto)],CTRhisto/np.max(CTRhisto),label="all events")  
# plt.xlabel("Esharing")
# plt.ylabel('#')
# plt.legend()
# plt.show()    
# ev_photopeak = event_registry_max_real[event_registry_max_real[:,4] >5700]

# LYratio=2.9/8.4
# t_decay=400
# t_gate1=50
# t_gate2=450

# EtotalPH= ev_photopeak[:,4]-ff*ev_photopeak[:,2]
# f1=np.exp(-t_gate1/t_decay)
# f2=np.exp(-t_gate2/t_decay)
# Eslow=(event_registry_max_real[:,4]-event_registry_max_real[:,3])/(f1-f2)
# Efast=(event_registry_max_real[:,4]-(1-f2)*Eslow)/LYratio
# Etotal=Eslow+Efast
# Esharing = Efast/Eslow
# CTRhisto,bins = np.histogram(Eslow, bins=80)
# plt.plot(bins[0:np.size(CTRhisto)], CTRhisto/np.max(CTRhisto),label="all events")      
# plt.xlabel("Eslow")
# plt.ylabel('#')
# plt.legend()
# plt.show()
# CTRhisto,bins = np.histogram(Efast, bins=50)
# plt.plot(bins[0:np.size(CTRhisto)], CTRhisto/np.max(CTRhisto),label="all events")  
# plt.xlabel("Efast")
# plt.ylabel('#')
# plt.legend()
# plt.show()   
