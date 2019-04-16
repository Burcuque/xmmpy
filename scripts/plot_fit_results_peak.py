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
t = Table.read('%s/PN_SW_fit_results.csv'%wdir)
#
fig = plt.figure(figsize=(10,10))
# Al line
ax1 = fig.add_subplot(311)
ax1.errorbar(t['revolution'],t['al_peak'],yerr=t['al_peak_err'],fmt='o')
#med = np.median(t['al_peak'])
#ax1.axhline(y=med,color='r',ls='dashed',label='median')
#ax1.set_ylim([58,60])
ax1.set_xlabel("revolution")
ax1.set_ylabel("Al line peak (cts/s)")
ax1.set_yscale('log')
ax1.grid(True)
ax1.legend(loc=1)
plt.title("PN CalClosed in Small Window (n=52)")
# Mn1 line
ax2 = fig.add_subplot(312)
ax2.errorbar(t['revolution'],t['mn1_peak'],yerr=t['mn1_peak_err'],fmt='o')
#med = np.median(t['mn1_peak'])
#ax2.axhline(y=med,color='r',ls='dashed',label='median')
ax2.set_xlabel("revolution")
ax2.set_ylabel("Mn1 line peak (cts/s)")
ax2.set_yscale('log')
ax2.legend(loc=1)
ax2.grid(True)
# Mn2
ax3 = fig.add_subplot(313)
ax3.errorbar(t['revolution'],t['mn2_peak'],yerr=t['mn2_peak_err'],fmt='o')
#med = np.median(t['mn2_peak'])
#ax3.axhline(y=med,color='r',ls='dashed',label='median')
ax3.set_xlabel("revolution")
ax3.set_ylabel("Mn2 line peak (cts/s)")
#ax3.set_ylim([0,0.4])
ax3.set_yscale('log')
ax3.legend(loc=1)
ax3.grid(True)

plt.savefig("{}/PN_SW_CalClosed_bin5_linePeak.png".format(wdir),dpi=100)
plt.close()
#plt.show()

