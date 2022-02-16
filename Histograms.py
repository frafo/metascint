# -*- coding: utf-8 -*-
"""
Created on Mon Apr 19 11:12:07 2021

@author: yocon
"""


import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm
from scipy.optimize import curve_fit
import codecs
import pandas as pd
import os


def histoplot(vector,bins,label):
    CTRhisto,bins = np.histogram(vector, bins=bins)
    normaliz=np.max(CTRhisto)        
    plt.plot(bins[0:np.size(CTRhisto)], CTRhisto/normaliz,label=label)
#"""

root = 'D:\Downloads\multiwave\Metapixel_datasets\BB15a_VUV'
file = 'signal_ref2.csv.npy'
fname = os.path.join(root,file)
contentfile2 = np.load(fname)

file = 'signal_ref3.csv.npy'
fname = os.path.join(root,file)
contentfile3 = np.load(fname)

file = 'signal_ref4.csv.npy'
fname = os.path.join(root,file)
contentfile4 = np.load(fname)

file = 'signal_ref5.csv.npy' #remember, this one must have been analysed with different mega_analysis. Versioning and logging!
fname = os.path.join(root,file)
contentfile1 = np.load(fname)

event_registry_max_real= np.concatenate((contentfile1,contentfile2, contentfile3, contentfile4), axis=0)
# event_registry_max_real =contentfile2
#setting filtering conditions
event_registry_max_real= event_registry_max_real[event_registry_max_real[:,8] != 0]
event_registry_max_real= event_registry_max_real[(event_registry_max_real[:,5]) <1000] 
event_registry_max_real= event_registry_max_real[(event_registry_max_real[:,5]) >-500]

#choose subsets
ev_photopeak = event_registry_max_real[event_registry_max_real[:,4] >9000] #photopeak

for i in range(9):
    # if i>1 and i<5:
    histoplot(event_registry_max_real[:,i], 50, "all events")
    histoplot(ev_photopeak[:,i],30, "photopeak")
    plt.xlabel(str(i))
    plt.ylabel('#')
    plt.legend()
    plt.show()
        
