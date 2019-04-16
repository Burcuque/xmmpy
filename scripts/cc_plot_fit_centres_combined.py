#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 17 15:55:36 2018

@author: ivaltchanov
"""

import os
import numpy as np

from astropy.table import Table

import matplotlib.pyplot as plt

plt.rc('text', usetex=False)
plt.rc('font', family='serif')

home = os.path.expanduser('~')
#wdir = home + '/XMM/CAL'
wdir = '/xdata/xcaldata/XMM/IVAN/PN_calclosed'

#
ff_results = "Jan2018/PN_FF_fit_results.csv"
t = Table.read('%s/%s'%(wdir,ff_results))
nt = len(t)
#
sw_results = "Jan2018/PN_SW_fit_results.csv"
tx = Table.read('%s/%s'%(wdir,sw_results))
nx = len(tx)
#
fig = plt.figure(figsize=(10,10))
# Al line
lab = 1.486 #eV
ax1 = fig.add_subplot(311)
ax1.errorbar(t['revolution'],(t['al_cen']/200.0 - lab)*1000.0,yerr=1000*t['al_cen_err']/200.0,fmt='o',label='PN FF, n={}'.format(nt))
ax1.errorbar(tx['revolution'],(tx['al_cen']/200.0 - lab)*1000,yerr=1000*tx['al_cen_err']/200.0,fmt='o',color='red',label='PN SW, n={}'.format(nx))
#med = np.median(t['al_cen'])
ax1.axhline(y=0.0,color='k',ls='dashed')
#ax1.axhline(y=med,color='r',ls='dashed',label='median')
#ax1.axhline(y=lab,color='r',ls='dashed',label='theoretical')
ax1.set_ylim([-60.0,60.0])
ax1.set_xlabel("revolution")
ax1.set_ylabel(r"Al K$_\alpha$ centre (eV)")
ax1.grid(True)
ax1.legend(loc=2)
#plt.title("PN CalClosed in Small Window (n=52)")
plt.title("PN CalClosed")
# Mn1 line
lab = (0.162*5.888 + 0.6*5.899)/0.762 #eV probabilities from Wikipedia on Iron-55
ax2 = fig.add_subplot(312)
ax2.errorbar(t['revolution'],(t['mn1_cen']/200.0 - lab)*1000.0,yerr=1000*t['mn1_cen_err']/200.0,fmt='o',label='PN FF, n={}'.format(nt))
ax2.errorbar(tx['revolution'],(tx['mn1_cen']/200.0 - lab)*1000,yerr=1000*tx['mn1_cen_err']/200.0,fmt='o',color='red',label='PN SW, n={}'.format(nx))
#med = np.median(t['al_cen'])
#ax1.axhline(y=med,color='r',ls='dashed',label='median')
#ax1.axhline(y=lab,color='r',ls='dashed',label='theoretical')
ax2.axhline(y=0.0,color='k',ls='dashed')
ax2.set_ylim([-60.0,60.0])
ax2.set_xlabel("revolution")
ax2.set_ylabel(r"Mn K$_{\alpha}$ centre (eV)")
ax2.grid(True)
ax2.legend(loc=2)
# Mn2
lab = 6.490 #eV
ax3 = fig.add_subplot(313)
ax3.errorbar(t['revolution'],(t['mn2_cen']/200.0 - lab)*1000.0,yerr=1000*t['mn2_cen_err']/200.0,fmt='o',label='PN FF, n={}'.format(nt))
ax3.errorbar(tx['revolution'],(tx['mn2_cen']/200.0 - lab)*1000,yerr=1000*tx['mn2_cen_err']/200.0,fmt='o',color='red',label='PN SW, n={}'.format(nx))
#med = np.median(t['al_cen'])
#ax1.axhline(y=med,color='r',ls='dashed',label='median')
#ax1.axhline(y=lab,color='r',ls='dashed',label='theoretical')
ax3.axhline(y=0.0,color='k',ls='dashed')
ax3.set_ylim([-60.0,60.0])
ax3.set_xlabel("revolution")
ax3.set_ylabel(r"Mn K$_{\beta}$ centre (eV)")
ax3.grid(True)
ax3.legend(loc=2)
#
plt.savefig("{}/PN_FF_SW_CalClosed_bin5_lineCentres.png".format(wdir),dpi=100)
plt.show()
#plt.savefig("{}/PN_SW_CalClosed_bin5_lineCentres.png".format(wdir),dpi=100)
plt.close()


