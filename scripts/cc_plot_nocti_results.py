#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 20 12:27:08 2018

Plot the line fit results, the difference is with respect to the 
redshifted Fe Kalpha line

@author: ivaltchanov
"""
import os

import numpy as np

from astropy.table import Table

import matplotlib.pylab as plt
#from matplotlib.patches import Patch
#from matplotlib.lines import Line2D

plt.rc('text', usetex=False)
plt.rc('font', family='serif')
#

wdir = "/xdata/xcaldata/XMM/IVAN/PN_calclosed"
os.chdir(wdir)

#%%
t = Table.read("pn_sw_calclosed_fit_results_nocti.csv")
nt = len(t)
#
rev = t['rev'].data
rtime = t['time'].data
#
line1_lab = 1.486 # keV
line1 = t['al_cen'].data/1000.0
line1Err = t['al_cen_err'].data/1000.0
#
line2_lab = (0.162*5.888 + 0.6*5.899)/0.762 #eV probabilities from Wikipedia on Iron-55
line2 = t['mn1_cen'].data/1000.0
line2Err = t['mn1_cen_err'].data/1000.0
#
line3_lab = 6.490 # keV
line3 = t['mn2_cen'].data/1000.0
line3Err = t['mn2_cen_err'].data/1000.0
#
#%%
msize=10
fig = plt.figure(figsize=(12,8))
ax = fig.subplots(3,sharex=True)

#
ax[0].errorbar(rtime,line1,yerr=([line1Err,line1Err]),fmt='o',\
           markersize=msize)
ax[0].axhline(line1_lab,color='k',ls='dashed')
ax[0].set_ylabel(r"Al Line Energy (keV)")
ax[0].grid(True)
#
ax[1].errorbar(rtime,line2,yerr=([line2Err,line2Err]),fmt='o',\
           markersize=msize)
ax[1].axhline(line2_lab,color='k',ls='dashed')
ax[1].set_ylabel(r"Mn K$_\alpha$ Line Energy (keV)")
ax[1].grid(True)
#
ax[2].errorbar(rtime,line3,yerr=([line3Err,line3Err]),fmt='o',\
           markersize=msize)
ax[2].axhline(line3_lab,color='k',ls='dashed')
ax[2].set_ylabel(r"Mn K$_\beta$ Line Energy (keV)")

ax[2].set_xlabel(r"Years since 2000-01-01T00:00:00")
ax[2].grid(True)
#fig.subplots_adjust(hspace=0)
#plt.savefig('{}/{}_fe_fit_results.png'.format(wdir,target),dpi=100)
plt.show()
#plt.close()
#
#fout.close()
