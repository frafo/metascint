# -*- coding: utf-8 -*-
"""
Created on Wed Jul 28 16:01:55 2021

@author: yocon
"""


import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm
from scipy.optimize import curve_fit


def namestr(obj, namespace):
    return [name for name in namespace if namespace[name] is obj]

def expon(x,weight,conv,size,offset):
    return weight*np.exp(-(x-offset)*size)+conv

def sort_fit(ev_photopeak, factor, capt):
    ev_ph_sort=ev_photopeak[ev_photopeak[:,factor].argsort()]
    plt.scatter(ev_ph_sort[:,factor],ev_ph_sort[:,5],dot_size)
    fit1 = np.polyfit(ev_ph_sort[:,factor],ev_ph_sort[:,5],1)
    trendline1 = np.poly1d(fit1)
    plt.plot(ev_ph_sort[:,factor],trendline1(ev_ph_sort[:,factor]),'--',label= capt + ' fit')
    return fit1

def sort_fit_exp(ev_photopeak, factor, capt):
    ev_ph_sort=ev_photopeak[ev_photopeak[:,factor].argsort()]
    plt.plot(ev_ph_sort[:,factor],ev_ph_sort[:,5])#, label= namestr(ev_photopeak, globals()))
    x = np.arange(np.min(ev_photopeak[:,factor]),np.max(ev_photopeak[:,factor]),np.size(ev_ph_sort))
    y =  np.array(ev_ph_sort[:,5])
    conv=100
    weight=100
    size=0.01
    offset=np.min(ev_photopeak[:,factor])
    popt,pcov = curve_fit(expon,x,y,p0=[weight,conv,size,offset])
    mx= np.argmax(expon(x,*popt))
    plt.plot(x-mx,expon(x-mx,*popt),label=capt)

fac=2
weight=1
Ph= sort_fit(ev_photopeak, fac,'ph')
# Fa=sort_fit(ev_fast,fac,'Fast')
Uf=sort_fit(ev_ultrafast,fac,'Ultrafast')
# Sh=sort_fit(ev_shared,fac,'Shared')
# BG=sort_fit(ev_BGO,fac,'BGO')
expo=sort_fit_exp(ev_photopeak, fac,'ph exponential')
# plt.xlim(1700,2700)

# x = np.arange(750,2000,10)
# y = 100*np.exp(-0.01*(x-750))+100
# plt.plot(x,y)
plt.legend()
plt.show()


