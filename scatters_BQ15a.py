# -*- coding: utf-8 -*-
"""
Created on Mon Apr 19 11:30:55 2021

@author: yocon
"""

import numpy as np
import matplotlib.pyplot as plt
#Now there should be a bit of second level analysis. Take the registry and normalize it, populate it, calculate energy content etc etc
#This config is for dataset BC15a full pulses
photopeak_cut=9700
# BGO_cut = 140
# BaF2_cut = 50
# cut1 = 85
ev_photopeak = event_registry_max_real[event_registry_max_real[:,4]>photopeak_cut]
ev_photopeak = ev_photopeak[ev_photopeak[:,5] >-400]
ev_photopeak = ev_photopeak[ev_photopeak[:,5] < 500]
# ev_photopeak = ev_photopeak[ev_photopeak[:,6] < 500]

# ev_compton=event_registry_max_real[event_registry_max_real[:,4] < photopeak_cut]
# ev_shared= ev_photopeak[ev_photopeak[:,6] <73]
# ev_BGO= ev_photopeak[ev_photopeak[:,6] >BGO_cut]
# ev_shared = ev_photopeak[ev_photopeak[:,6] <BGO_cut+1]
# ev_shared = ev_shared[ev_shared[:,6] > cut1]
# ev_shared1 = ev_photopeak[ev_photopeak[:,6] < cut1+1]
# ev_shared1 = ev_shared1[ev_shared1[:,6] >BaF2_cut]
# ev_BaF2 = ev_photopeak[ev_photopeak[:,6] < BaF2_cut+1]

# #Choose values using the total energy
# BGO_cut =10600
# Share_cut = 11500
# BaF2_cut =12400
# ev_BGO= ev_photopeak[ev_photopeak[:,4] <BGO_cut]

# ev_shared_less = ev_photopeak[ev_photopeak[:,4] >BGO_cut]
# ev_shared_less = ev_shared_less[ev_shared_less[:,4] <Share_cut+1]
# ev_shared_more = ev_photopeak[ev_photopeak[:,4] >Share_cut]
# ev_shared_more = ev_shared_more[ev_shared_more[:,4] <BaF2_cut+1]

# ev_BaF2 = ev_photopeak[ev_photopeak[:,4] > BaF2_cut+1]

#Choose values using the max energy
fac1=3
BGO_cut = 410
fast_cut = 580
ev_photopeak = event_registry_max_real[event_registry_max_real[:,4] >5500]
ev_photopeak = ev_photopeak[ev_photopeak[:,4]<7500]
NLim1=0
Lim1=250
Lim2= 500


# ev_fast1 = ev_photopeak[ev_photopeak[:,3]<0.115*ev_photopeak[:,4]-Lim3]
# ev_fast1 = ev_fast1[ev_fast1[:,3]>0.115*ev_fast1[:,4]-Lim1]
# ev_fast2 = ev_photopeak[ev_photopeak[:,3]<0.115*ev_photopeak[:,4]-Lim1]
# ev_fast2 = ev_fast2[ev_fast2[:,3]>1050]
# ev_fast = np.concatenate((ev_fast1, ev_fast2), axis=0)

ev_ultrafast = ev_photopeak[ev_photopeak[:,3]>0.115*ev_photopeak[:,4]-NLim1]
ev_fast = ev_photopeak[ev_photopeak[:,3]<0.115*ev_photopeak[:,4]]
ev_fast = ev_fast[ev_fast[:,3]>0.115*ev_fast[:,4]-Lim1]
ev_shared = ev_photopeak[ev_photopeak[:,fac1] < 0.115*ev_photopeak[:,4]-Lim1]
ev_shared = ev_shared[ev_shared[:,fac1] > 0.105*ev_shared[:,4]-Lim2]
ev_BGO = ev_photopeak[ev_photopeak[:,fac1] <0.105*ev_photopeak[:,4]-Lim2]



dot_size=3
# Here I am choosing data based on their slope and total energy
ref_var=4

#rise time scatters
# # plt.scatter(ev_compton[:,ref_var],ev_compton[:,6]/100, dot_size, label="Compton shared")
# plt.scatter(ev_photopeak[:,ref_var],ev_photopeak[:,6]/100, dot_size, label="Photoelectric")
plt.scatter(ev_ultrafast[:,ref_var],ev_ultrafast[:,2], dot_size, label="ultrafast")
plt.scatter(ev_fast[:,ref_var],ev_fast[:,2], dot_size, label="fast")

plt.scatter(ev_shared[:,ref_var],ev_shared[:,2], dot_size, label="Shared mid")
plt.scatter(ev_BGO[:,ref_var],ev_BGO[:,2], dot_size, label="BGO")

plt.ylabel('Gate 90 ns (A.U.)')
plt.xlabel('Gate 450 ns (A.U.)')
# plt.xlim(photopeak_cut,14000)
plt.legend()
plt.show()

#CTR scatters
# plt.scatter(ev_compton[:,ref_var],ev_compton[:,5], dot_size, label="Compton shared")
# plt.scatter(ev_photopeak[:,ref_var],ev_photopeak[:,5], dot_size, label="Photoelectric")
plt.scatter(ev_ultrafast[:,ref_var],ev_ultrafast[:,5], dot_size, label="ultrafast")
plt.scatter(ev_fast[:,ref_var],ev_fast[:,5], dot_size, label="fast")

plt.scatter(ev_shared[:,ref_var],ev_shared[:,5], dot_size, label="Shared")
plt.scatter(ev_BGO[:,ref_var],ev_BGO[:,5], dot_size, label="BGO")

plt.ylabel('CTR (ps)')
plt.xlabel('Gate 450 ns (A.U.)')
# plt.xlim(photopeak_cut,14000)
# plt.ylim(-500,500)
plt.legend()
plt.show()

# #max energy scatters
# # # plt.scatter(ev_compton[:,ref_var],ev_compton[:,6]/100, dot_size, label="Compton shared")
# # plt.scatter(ev_photopeak[:,ref_var],ev_photopeak[:,6]/100, dot_size, label="Photoelectric")
# plt.scatter(ev_ultrafast[:,ref_var],ev_ultrafast[:,3], dot_size, label="ultrafast")
# plt.scatter(ev_fast[:,ref_var],ev_fast[:,3], dot_size, label="fast")

# plt.scatter(ev_shared[:,ref_var],ev_shared[:,3], dot_size, label="Shared")
# plt.scatter(ev_BGO[:,ref_var],ev_BGO[:,3], dot_size, label="BGO")


# # # Slope scatters 6 ,8
Fig1=5
plt.scatter(ev_ultrafast[:,ref_var],ev_ultrafast[:,Fig1]*10, dot_size, label="ultrafast")
plt.scatter(ev_fast[:,ref_var],ev_fast[:,Fig1]*10, dot_size, label="fast")

plt.scatter(ev_shared[:,ref_var],ev_shared[:,Fig1]*10, dot_size, label="Shared mid")
plt.scatter(ev_BGO[:,ref_var],ev_BGO[:,Fig1]*10, dot_size, label="BGO")

plt.ylabel('rise time (ps)')
plt.xlabel('Gate 450 ns (A.U.)')
plt.legend()
plt.show()

# Fig2=8
# plt.scatter(ev_ultrafast[:,ref_var],ev_ultrafast[:,Fig2]-ev_ultrafast[:,Fig1], dot_size, label="ultrafast")
# plt.scatter(ev_fast[:,ref_var],ev_fast[:,Fig2]-ev_fast[:,Fig1], dot_size, label="fast")

# plt.scatter(ev_shared[:,ref_var],ev_shared[:,Fig2]-ev_shared[:,Fig1], dot_size, label="Shared mid")
# plt.scatter(ev_BGO[:,ref_var],ev_BGO[:,Fig2]-ev_BGO[:,Fig1], dot_size, label="BGO")

# plt.ylabel('Gate 90 ns (A.U.)')
# plt.xlabel('Gate 450 ns (A.U.)')
# plt.legend()
# plt.show()

# # plt.ylim(0.2,0.7)
# plt.ylabel('Gate 50 ns (A.U.)')
# # plt.xlim(photopeak_cut,14000)
# plt.xlabel('Gate 450 ns (A.U.)')
# plt.legend()
# plt.show()
# for i in range(9):
#     if i>4 and i<7:
#         plt.scatter(event_registry_max_real[:,ref_var],event_registry_max_real[:,i]/100, dot_size, label="BGO")
#         plt.scatter(ev_shared[:,ref_var],ev_shared[:,i]/100, dot_size, label="Compton shared")
#         plt.scatter(ev_photopeak[:,ref_var],ev_photopeak[:,i]/100, dot_size, label="Photoelectric shared")
#         plt.scatter(ev_fast3[:,ref_var],ev_fast3[:,i]/100, dot_size, label="BaF2")


# #        plt.scatter(ev_fast1[:,ref_var],ev_fast1[:,i], dot_size, label="fast1")
# #        plt.scatter(ev_fast2[:,ref_var],ev_fast2[:,i], dot_size, label="fast2")
# #        plt.scatter(ev_fast3[:,ref_var],ev_fast3[:,i], dot_size, label="fast3")
#         plt.xlabel('Energy (A.U.)')
#         plt.ylabel(str(i))
# #        plt.xlim(0,0.001)
#         if i==5:
#             plt.ylim(0,3)
#             plt.ylabel('CTR (ns)')
#         if i==3:
#             plt.ylim(0,0.00015)
# #            plt.xlim(10000,19000)
#         if i==6:
#             plt.ylim(0,2)
#             plt.ylabel('rise time to 0.18 mV (ns)')
# #        if i==1:
# #            x=np.arange(10000,19000,100)
# #            y=0.000107*x-0.39
# #            plt.plot(x,y,'r--',linewidth=3)
# #            x=np.arange(9000,14500,100)
# #            y=0.000107*x+0.18
# #            plt.plot(x,y,'r--',linewidth=3)
# #            x=np.arange(7000,11000,100)            
# #            y=0.000107*x+0.55
# #            plt.plot(x,y,'r--',linewidth=3)
#         plt.legend()
#         plt.show()
'''
#fit1 = np.polyfit(ev_shared[:,4],0.001*np.log(ev_shared[:,1])+4,1)
#trendline1 = np.poly1d(fit1)
x = np.arange(50,350,1)
y=0.015*np.log(x-40)-0.013
plt.plot(x,y,'r--',linewidth=3)
x = np.arange(75,200,1)
y=0.015*np.log(x+100)-0.008
plt.plot(x,y,'r--',linewidth=3)
#x = np.arange(25,400,1)
#y = -0.000032*x+0.082
#plt.plot(x,y,'r--',linewidth=3)
x = np.arange(25,230,1)
y = -0.000012*x+0.069
plt.plot(x,y,'r--',linewidth=3)
#x = np.arange(100,295,1)
#y = -0.000032*x+0.0755
#plt.plot(x,y,'b--',linewidth=3)

#'''