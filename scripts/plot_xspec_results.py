#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 20 12:27:08 2018

@author: ivaltchanov
"""
import os

import numpy as np

from astropy.table import Table

import matplotlib.pylab as plt

plt.rc('text', usetex=False)
plt.rc('font', family='serif')
#
#target = "mrk1048"
target = "WR140"
wdir = "/home/ivaltchanov/XMM/CAL/PN_SW/{}".format(target)
#wdir = "/home/ivaltchanov/XMM/CAL/PN_SW/WR140"
os.chdir(wdir)
#%%
t = Table.read('output_xspec.csv')
nt = len(t)
#
rev = t['revol'].data
ff = t['ff'].data
inst = t['instrument'].data
line = t['line'].data
lineErrUp = t['lineErr_up'].data
lineErrLo = t['lineErr_lo'].data
#
# now set up the windowed modes
##
#
m1 = (inst == 'mos1')
i1 = np.where(m1)[0]
i1x = np.where(m1*ff)[0]
m2 = (inst == 'mos2')
i2 = np.where(m2)[0]
i2x = np.where(m2*ff)[0]
pn = (inst == 'pn')
i3 = np.where(pn)[0]
i3x = np.where(pn*ff)[0]
#
# now the ratio of all lines to MOS1
#

#
msize=10
fig = plt.figure(figsize=(10,8))
ax = fig.subplots()
#
#ax.errorbar(np.arange(len(rev[i1])),line[i1],yerr=(lineErrLo[i1],lineErrUp[i1]),fmt='-o',markersize=msize,\
#            color='blue',fillstyle='none',label='MOS1')
#ax.errorbar(np.arange(len(rev[i2]))+0.1,line[i2],yerr=(lineErrLo[i2],lineErrUp[i2]),fmt='-o',markersize=msize,\
#            color='green',fillstyle='none',label='MOS2')
#ax.errorbar(np.arange(len(rev[i3]))+0.2,line[i3],yerr=(lineErrLo[i3],lineErrUp[i3]),fmt='-o',markersize=msize,\
#            color='red',fillstyle='none',label='PN')
#
#ax.errorbar(np.arange(len(rev[i1x])),line[i1x],yerr=(lineErrLo[i1x],lineErrUp[i1x]),fmt='-o',markersize=msize,\
#            color='blue',fillstyle='full')
#ax.errorbar(np.arange(len(rev[i2x]))+0.1,line[i2x],yerr=(lineErrLo[i2x],lineErrUp[i2x]),fmt='-o',markersize=msize,\
#            color='green',fillstyle='full')
#ax.errorbar(np.arange(len(rev[i3x]))+0.2,line[i3x],yerr=(lineErrLo[i3x],lineErrUp[i3x]),fmt='-o',markersize=msize,\
#            color='red',fillstyle='full')
#
#for j in np.arange(len(rev[i1])):
#    ax.text(j,6.05,rev[i1[j]])
#
ax.errorbar(rev[i1],line[i1],yerr=(lineErrLo[i1],lineErrUp[i1]),fmt='-o',markersize=msize,\
            color='blue',fillstyle='none',label='MOS1')
ax.errorbar(rev[i2]+1,line[i2],yerr=(lineErrLo[i2],lineErrUp[i2]),fmt='-o',markersize=msize,\
            color='green',fillstyle='none',label='MOS2')
ax.errorbar(rev[i3]+2,line[i3],yerr=(lineErrLo[i3],lineErrUp[i3]),fmt='-o',markersize=msize,\
            color='red',fillstyle='none',label='PN')

ax.errorbar(rev[i1x],line[i1x],yerr=(lineErrLo[i1x],lineErrUp[i1x]),fmt='-o',markersize=msize,\
            color='blue',fillstyle='full')
ax.errorbar(rev[i2x]+1,line[i2x],yerr=(lineErrLo[i2x],lineErrUp[i2x]),fmt='-o',markersize=msize,\
            color='green',fillstyle='full')
ax.errorbar(rev[i3x]+2,line[i3x],yerr=(lineErrLo[i3x],lineErrUp[i3x]),fmt='-o',markersize=msize,\
            color='red',fillstyle='full')
#
ax.set_ylabel(r"Fe Line Energy (keV)")
ax.set_xlabel("Revolution index")
ax.set_title("{} Analysis".format(target.capitalize()))
#ax.set_title("WR140 Analysis")
#
ay = fig.subplots()
ay.plot(rev[i1],line[i1]/line[i1],fmt='-o',markersize=msize,\
            color='blue',fillstyle='none',label='MOS1')
ay.plot(rev[i2]+1,line[i2]/line[i1],fmt='-o',markersize=msize,\
            color='green',fillstyle='none',label='MOS2')
ay.plot(rev[i3]+2,line[i3]/line[i1],fmt='-o',markersize=msize,\
            color='red',fillstyle='none',label='PN')

ay.errorbar(rev[i1x],line[i1x]/line[i1x],yerr=(lineErrLo[i1x],lineErrUp[i1x]),fmt='-o',markersize=msize,\
            color='blue',fillstyle='full')
ay.errorbar(rev[i2x]+1,line[i2x],yerr=(lineErrLo[i2x],lineErrUp[i2x]),fmt='-o',markersize=msize,\
            color='green',fillstyle='full')
ay.errorbar(rev[i3x]+2,line[i3x],yerr=(lineErrLo[i3x],lineErrUp[i3x]),fmt='-o',markersize=msize,\
            color='red',fillstyle='full')
#
plt.grid(True)
plt.legend(loc=2)
#plt.savefig('{}/{}_fe_fit_results.png'.format(wdir,target),dpi=100)
plt.show()
#plt.close()
#



