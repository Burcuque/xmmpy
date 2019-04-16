#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 27 16:46:34 2018

Using the NGC4151 results in SW mode for PN
comparing CTI v48 and v49 results 

plot results

@author: ivaltchanov
"""
import numpy as np

from astropy.table import Table
import matplotlib.pylab as plt
from matplotlib.lines import Line2D

#wdir = "/home/ivaltchanov/XMM/CAL/PN_SW"

wdir = "/xdata/xcaldata/XMM/IVAN/PN_SW/ngc4151"
#
feK = 6.399 # keV
#
# NGC 4151
redshift =  0.003262
lineX = feK/(1.0 + redshift)
#
v48all = Table.read('{}/output_xspec.csv'.format(wdir))
# need to pick up pn and SmallWindow
#
m1 = v48all['instrument'] == "pn"
m2 = v48all['ff'] == 0
ix = np.where(m1*m2)[0]
v48=v48all[ix]
#
v49 = Table.read('{}/output_xspec_cti49.csv'.format(wdir))
v49y = Table.read('{}/output_xspec_cti49y.csv'.format(wdir))

#%%
msize=10
fig = plt.figure(figsize=(12,8))
ax = fig.subplots()
#
# Current CTI v48 results
#
xrev48 = v48['revol']
#xrev48 = v48['delta_time']
line48 = v48['line']
line48_errL = v48['lineErr_lo']*1000
line48_errU = v48['lineErr_up']*1000
delta_v48 = (line48 - lineX)*1000.0
#
#xrev49 = v49['rev']
xrev49 = v49['delta_time']
line49 = v49['lineE']
line49_errL = v49['lineE_low']*1000
line49_errU = v49['lineE_up']*1000
delta_v49 = (line49 - lineX)*1000.0
#
xrev49y = v49y['delta_time']
line49y = v49y['lineE']
line49y_errL = v49y['lineE_low']*1000
line49y_errU = v49y['lineE_up']*1000
delta_v49y = (line49y - lineX)*1000.0
#
#ax.errorbar(xrev49,delta_v48,yerr=(line48_errL,line48_errU),marker='o',markersize=15,color='blue',
#            fillstyle='none',linestyle='none',label='CTI CCF v48')
ax.errorbar(xrev49,delta_v49,yerr=(line49_errL,line49_errU),marker='o',markersize=15,color='red',
            fillstyle='full',linestyle='none',label='CTI CCF v49')
ax.errorbar(xrev49y,delta_v49y,yerr=(line49y_errL,line49y_errU),marker='o',markersize=15,color='blue',
            fillstyle='none',linestyle='none',label='CTI CCF v49')
#ax.plot(xrev49,delta_v48,marker='o',markersize=15,color='blue',
#            fillstyle='none',linestyle='none')
#ax.plot(xrev49,delta_v49,marker='o',markersize=15,color='red',
#            fillstyle='full',linestyle='none')
ax.axhline(0.0,color='k',ls='dashed')
#
ax.grid(True)
ax.set_xlabel(r"Years since 2000-01-01T00:00:00")
ax.set_ylabel(r"$\Delta$ Fe K$\alpha$ (obs-lab) (eV)")
ax.set_title("NGC4151 in PN SW mode")
ax.legend(loc=2)

plt.savefig(f'{wdir}/ngc4151_pnsw_cti49y_results.png',dpi=100)
plt.show()
plt.close()
