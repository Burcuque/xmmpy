#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 23 17:17:07 2019

Check NGC1566 two observations

@author: ivaltchanov
"""
import os

import numpy as np

from astropy.table import Table
import matplotlib
import matplotlib.pylab as plt
import seaborn as sns
sns.set(style="white")

matplotlib.rcParams['axes.formatter.useoffset'] = False
plt.rc('text', usetex=False)
plt.rc('font', family='serif')
#
#%%
#
wdir = '/home/ivaltchanov/IVAN/PN_SW/sources/ngc1566'
os.chdir(wdir)

t1 = Table.read('0800840201_xspec_dump.csv')
t2 = Table.read('0820530401_xspec_dump.csv')
#
redshift = 0.005036
fek = 6.399/(1.0+redshift)
#
ax.axvline(fek,color='black',linestyle='--')
#%%
#
fig = plt.figure(figsize=(12,8))
ax = fig.add_subplot(111)
xx1 = t1['energy']
yy1 = t1['norm']
ix = np.where((xx1 <= 6.0) + (xx1 >= 6.7))[0]
#
z = np.polyfit(xx1[ix], yy1[ix], 1)
p = np.poly1d(z)
zz1 = p(xx1)
#
ax.errorbar(xx1,yy1,xerr=(t1['energyErr'],t1['energyErr']),\
            yerr=(t1['normErr'],t1['normErr']),fmt='+',label='0800840201',color='blue')
ax.plot(xx1,t1['model'],color='blue',label='')
#ax.errorbar(xx1,yy1-zz1,xerr=(t1['energyErr'],t1['energyErr']),\
#            yerr=(t1['normErr'],t1['normErr']),fmt='+',label='0800840201',color='blue')
#ax.plot(xx1,t1['model']-zz1,color='blue')
#

xx2 = t2['energy']
yy2 = t2['norm']
ix = np.where((xx2 <= 6.0) + (xx2 >= 6.7))[0]
#
z = np.polyfit(xx2[ix], yy2[ix], 1)
p = np.poly1d(z)
zz2 = p(xx1)
ax.errorbar(xx2,yy2,xerr=(t2['energyErr'],t2['energyErr']),\
            yerr=(t2['normErr'],t2['normErr']),fmt='+',label='0820530401',color='red')
ax.plot(xx2,t2['model'],color='red',label='')
#ax.errorbar(xx2,yy2-zz2,xerr=(t2['energyErr'],t2['energyErr']),\
#            yerr=(t2['normErr'],t2['normErr']),fmt='+',label='0820530401',color='red')
#ax.plot(xx2,t2['model']-zz2,color='red')
#
#
ax.axvline(fek,color='black',linestyle='--')
#
t = Table.read('ngc1566_output_xspec_cti49_t28.csv',data_start=0,\
               names=("obsid","rev","delta_time","submode",\
                                  "full","xfilt","inst","ontime",\
                                  "lineE","lineE_low","lineE_up","cstat","chi2r","dof"))

ax.axvspan(t['lineE'][0] - t['lineE_low'][0],t['lineE'][0] + t['lineE_up'][0],alpha=0.2,color='red')
ax.axvspan(t['lineE'][1] - t['lineE_low'][1],t['lineE'][1] + t['lineE_up'][1],alpha=0.2,color='blue')

ax.grid(True)
ax.legend()
ax.set_xlim((5.5,7.5))
ax.set_xlabel("Energy (keV)")
ax.set_ylabel("Normalised (counts/s/keV)")
plt.title("NGC 1566")
plt.savefig('ngc15644_FeKa_comparison.png',dpi=100)
plt.show()
plt.close()

