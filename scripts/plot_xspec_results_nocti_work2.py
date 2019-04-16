#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 20 12:27:08 2018

Plot the line fit results, the difference is with respect to the 
redshifted Fe Kalpha line

Will read the results form selected OBS_ID in a table.

@author: ivaltchanov
"""
import os

import numpy as np

from astropy.table import Table

import matplotlib.pylab as plt
#from matplotlib.patches import Patch
from matplotlib.lines import Line2D

plt.rc('text', usetex=False)
plt.rc('font', family='serif')
#
#target = "ngc4593"
#target = "ngc4151"
#target = "mrk1048"
#target = "ngc3227"
#target = "ngc3783"
#target = "ngc5506"
#target = "mcg-5-23-16"
#target = "ngc5548"
#target = "ngc3516"
#target = "WR140"
#targets = ["ngc5506","ngc5548", "ngc4151","ngc3516","ngc2992","ngc3227",'ngc3783','ngc1566']
#
# redshifts
#
redshift = {'ngc4151': 0.003262, 'ngc3227': 0.00386, 'mrk1048': 0.0427, 'ngc3783': 0.009755,\
            'ngc4593': 0.008344, 'ngc5506': 0.00589, 'mcg-5-23-16': 0.008226, 'ngc3516': 0.008816,\
            'ngc5548': 0.01627, 'ngc2992':  0.007296, 'ngc1566': 0.005036}
feK = 6.399 # the weitghed mean of Kalpha_2 at 6.3908 (intensity 50) and Kalpha_1 at 6.40308 (intensity 100)
#if (target == 'WR140'):
#    lineX = 6.697 # see Sugawara et al. (2015)
#
wdir = "/xdata/xcaldata/XMM/IVAN/PN_SW"
#wdir = "/home/ivaltchanov/XMM/CAL/PN_SW/[]".format(target)
#wdir = "/home/ivaltchanov/XMM/CAL/PN_SW/WR140"
os.chdir(wdir)
#
#%%
#
# read the table
#
tx = Table.read(f"{wdir}/selected_targets_nocti.csv")
targets = tx['target']
nt = len(targets)
#
#%%
msize=10
fig = plt.figure(figsize=(12,8))
ax = fig.subplots()
#
rev = []
rtime = []
line = []
lineErrUp = []
lineErrLo = []
ratio = []
ratioErrUp = []
ratioErrLo = []
#
for i in np.arange(nt):
    #
    target = targets[i]
    lineX =  feK/(1.0 + redshift[target]) # redshifted line position
    rev.append(tx['rev'].data[i])
    rtime.append(tx['delta_time'].data[i])
    line.append(tx['lineE'].data[i])
    lineErrLo.append(tx['lineE_low'].data[i])
    lineErrUp.append(tx['lineE_up'].data[i])
    ratio.append(tx['lineE'].data[i]/lineX)
    ratioErrLo.append(tx['lineE_low'].data[i]/lineX)
    ratioErrUp.append(tx['lineE_up'].data[i]/lineX)
#
ax.errorbar(rtime,ratio,yerr=([ratioErrLo,ratioErrUp]),fmt='o',\
               markersize=msize)
#ax.plot(rtime,ratio,'o')
#
ax.axhline(1.0,color='k',ls='dashed')
#ax.set_ylabel(r"Fe Line Energy (keV)")
ax.set_ylabel("Ratio")
ax.set_xlabel(r"Years since 2000-01-01T00:00:00")
#ax.set_title("[] Analysis".format(target.capitalize()))
ax.grid(True)
#
# Now add the calclosed data, only Mn Kalpha double results
#
t = Table.read("/xdata/xcaldata/XMM/IVAN/PN_calclosed/pn_sw_calclosed_fit_results_nocti.csv")
nt = len(t)
#
rev = t['rev'].data
rtime = t['time'].data
#
line2_lab = (0.162*5.888 + 0.6*5.899)/0.762 #eV probabilities from Wikipedia on Iron-55
line2 = t['mn1_cen'].data/1000.0
line2Err = t['mn1_cen_err'].data/1000.0
ax.errorbar(rtime,line2/line2_lab,yerr=([line2Err/line2_lab,line2Err/line2_lab]),fmt='^',\
           markersize=msize,color='black', label='PN SW calclosed')

#fig.subplots_adjust(hspace=0)
plt.legend()
#plt.savefig('[]/pn_sw_nocti_results.png'.format(wdir,target),dpi=100)
ax.set_title("PN SmallWindow Energy Scale Analysis\n Long-term CTI off")
plt.show()
#plt.close()
#
#fout.close()
