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

home = os.path.expanduser('~')
wdir = home + '/XMM/CAL/calclosed'
#
fitfile = "PN_FF_fit_results.csv"
t = Table.read('%s/%s'%(wdir,fitfile))
nx = len(t)
#t = Table.read('%s/PN_SW_fit_results.csv'%wdir)
#
fig = plt.figure(figsize=(10,10))
# Al line
ax1 = fig.add_subplot(311)
ax1.errorbar(t['revolution'],t['al_fwhm'],yerr=t['al_fwhm_err'],fmt='o')
med = np.median(t['al_fwhm'])
ax1.axhline(y=med,color='r',ls='dashed',label='median')
ax1.set_ylim([med-15,med+15])
ax1.set_xlabel("revolution")
ax1.set_ylabel("Al line FWHM (channel)")
ax1.grid(True)
ax1.legend(loc=2)
plt.title("{} (n={})".format(fitfile,nx))
# Mn1 line
ax2 = fig.add_subplot(312)
ax2.errorbar(t['revolution'],t['mn1_fwhm'],yerr=t['mn1_fwhm_err'],fmt='o')
med = np.median(t['mn1_fwhm'])
ax2.axhline(y=med,color='r',ls='dashed',label='median')
ax2.set_ylim([med-15,med+15])
ax2.set_xlabel("revolution")
ax2.set_ylabel("Mn1 line FWHM (channel)")
ax2.legend(loc=2)
ax2.grid(True)
# Mn2
ax3 = fig.add_subplot(313)
ax3.errorbar(t['revolution'],t['mn2_fwhm'],yerr=t['mn2_fwhm_err'],fmt='o')
med = np.median(t['mn2_fwhm'])
ax3.axhline(y=med,color='r',ls='dashed',label='median')
ax3.set_ylim([med-15,med+15])
ax3.set_xlabel("revolution")
ax3.set_ylabel("Mn2 line FWHM (channel)")
#ax3.set_ylim([5,12])
ax3.legend(loc=2)
ax3.grid(True)

plt.savefig("{}/PN_FF_CalClosed_bin5_lineFwhm.png".format(wdir),dpi=100)
#plt.savefig("{}/PN_SW_CalClosed_bin5_lineFwhm.png".format(wdir),dpi=100)
plt.close()
#plt.show()

