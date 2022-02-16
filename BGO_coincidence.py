# -*- coding: utf-8 -*-
"""
Created on Wed Mar 31 11:11:05 2021
CTR Timing from metapixels
@author: yocon
"""

import os
import numpy as np
from scipy.stats import norm
import matplotlib.pyplot as plt
import csv
import pylab as plb
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy import asarray as ar,exp
import codecs
import pandas as pd

#load data from file. Start with lite for ease.
#"""

Snap = 100000    #this is the frame of each event, corresponds to 500 ns with 10 ps resolution
reso = 10       #oscilloscope resolution in ps
file = 'test_100.csv'
root = 'D:\Downloads\multiwave\Data_acquisitions_2022\BGO15mm'
fname = os.path.join(root,file)
df = pd.read_csv(fname, encoding ='utf-8', sep=' ')
 # = pad.concat(chunk)
content=df.values

t1 = np.array(content[:,0])
e1 = np.array(content[:,1])
t2 = np.array(content[:,2])
e2 = np.array(content[:,1])

#"""

# """
#Find the framing. 
t1_start=0
t2_start=0
length_chain=int(np.size(t1)/Snap) #this is the number of events
event_registry = np.zeros((length_chain, 9)) #initialize pulse features vector
#event framework

#Get features. Max values of both pulses, energy content etc
    
for i in range(length_chain):
    startROI = (Snap*i)
    startROI_sh=startROI #add here if frame is not well adjusted
    endROI = (Snap*(i+1)-1)
    endROI_sh= startROI+15000 #1/10 of the frame, just the rising slope
    eventID=i    
    # Take noise levels at the beginning of each event, for callibration and baselining
    noise_bl_t1 = np.max(t1[startROI:startROI+2000])
    noise_bl_t2 = np.max(t2[startROI:startROI+2000])
    noise_bl_e1 = np.max(e1[startROI:startROI+2000])
    noise_bl_e2 = np.max(e2[startROI:startROI+2000])
    t1_max = np.max(t1[startROI:endROI])
    e1_max = np.max(e1[startROI:endROI])
    t2_max = np.max(t2[startROI:endROI])
    e2_max = np.max(e2[startROI:endROI])

    t1_area = sum(y for y in t1[startROI:endROI] if y >noise_bl_t1)
    e1_area = sum(y for y in e1[startROI:endROI] if y >noise_bl_e1)
    t2_area = sum(y for y in t2[startROI:endROI] if y >noise_bl_t2)
    e2_area = sum(y for y in e2[startROI:endROI] if y >noise_bl_e2)   
      
    
#Find an algorithm to go down to the timing pulse onset to establish a CTR value
    t1_stop1 = next((x[0] for x in enumerate(t1[startROI:endROI]) if x[1] > 0.017),0) #High threshold
    t1_star = next((x[0] for x in enumerate(reversed(t1[startROI:startROI_sh+t1_stop1])) if x[1] < noise_bl_t1),0) #stepping back till noise
    t1_start = t1_stop1 -t1_star #integrating the results

    t2_stop1 = next((x[0] for x in enumerate(t2[startROI:endROI]) if x[1] > 0.017),0)
    t2_star = next((x[0] for x in enumerate(reversed(t2[startROI:startROI_sh+t2_stop])) if x[1] < noise_bl_t2),0)
    t2_start= t2_stop1 -t2_star
      
    CTR = (reso*(t1_start-t2_start)) #put here the resolution in ps

    e1_start = next((x[0] for x in enumerate(e1[startROI_sh:endROI] ) if x[1] > noise_bl_e1+0.005),0) #beginning of energy pulse for integration windows
    e1_area1 = sum(y for y in e1[(startROI+e1_start):(startROI+e1_start+5000)] if y >noise_bl_e1) #first 50 ns of the pulse
    e1_area2 = sum(y for y in e1[(startROI+e1_start):(startROI+e1_start+9000)] if y >noise_bl_e1) #first 90 ns of the pulse
    
    e2_start = next((x[0] for x in enumerate(e2[startROI_sh:endROI] ) if x[1] > noise_bl_e2+0.005),0) #beginning of energy pulse for integration windows
    e2_area1 = sum(y for y in e2[(startROI+e2_start):(startROI+e1_start+5000)] if y >noise_bl_e2) #first 50 ns of the pulse
    e2_area2 = sum(y for y in e2[(startROI+e2_start):(startROI+e1_start+9000)] if y >noise_bl_e2) #first 90 ns of the pulse
    
#Load all information in a database per event.
#    if (abs(CTR)<500):# and (e1_area <20000): #and (t1_slope <1000): # in case there are known features to remove from registry
#                        0         1       2        3        4        5         6        7         8
    event_registry[i]=[eventID, e1_max, e1_area2, e1_area1, e1_area, CTR, t1_star , e2_area, t1_area]
event_registry_max_real= event_registry[np.all(event_registry[:,6:8] != 0, axis=1)]


np.save(fname+'.npy', event_registry_max_real) #save db in anpy array to avoid repeating the process

CTRhisto,bins = np.histogram(event_registry_max_real[:,4], bins=40) # print a histo to check how things look 
normaliz=np.max(CTRhisto)        
plt.plot(bins[0:np.size(CTRhisto)], CTRhisto/np.max(CTRhisto))
plt.xlabel("energy (A.U.)")
plt.show()
#"""
                                            #plot event i
event_id=287                # print an event to check how things look 
startROI = (Snap*event_id)
endROI = startROI +50000
t1a = np.array(content[startROI:endROI,0])
e1a = np.array(content[startROI:endROI,1])
t2a = np.array(content[startROI:endROI,2])
e2a = np.array(content[startROI:endROI,1])

x = np.arange(0,np.size(t1a),1) #4 picosecond of resolution
plt.plot(x,t1a,label="t1")
plt.plot(x,e1a,label="e1")
plt.plot(x,t2a,label="t2")
plt.plot(x,e2a,label="e2")
plt.xlabel('time (ps)')
plt.ylabel('Voltage (V)')
#plt.ylim(0.170,0.19)
plt.legend()
plt.show()
#"""
