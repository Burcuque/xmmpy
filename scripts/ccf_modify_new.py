#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov  9 10:36:31 2018

Modify CCF file "by hand "

@author: ivaltchanov
"""
import os
from astropy.io import fits
import numpy as np
import matplotlib.pylab as plt
#from matplotlib.patches import Patch
#from matplotlib.lines import Line2D

import seaborn as sns
sns.set(style="white")

plt.rc('text', usetex=False)
plt.rc('font', family='serif')
#%%
wdir = "/home/ivaltchanov/IVAN/PN_SW/ccfdev"

os.chdir(wdir)
ccf49_file = 'EPN_CTI_0049.CCF_test25'
#
ccf48_file = '/xdata/ccf/valid/EPN_CTI_0048.CCF'

hdu48 = fits.open(f"{ccf48_file}")
ltc48 = hdu48['LONG_TERM_CTI']
trel = hdu48['LTC_TIMES'].data["TIME"][0]
ix48 = np.where(ltc48.data['MODE_ID'] == 3)[0]
tcoeff = ltc48.data['T_COEFF']
gt48 = tcoeff[ix48[1]]
print (gt48)

hdu49 = fits.open(f"{ccf49_file}")
ltc49 = hdu49['LONG_TERM_CTI']
ix49 = np.where(ltc49.data['MODE_ID'] == 3)[0]

tcoeff = ltc49.data['T_COEFF']
gt49 = tcoeff[ix49[1]]
print (gt49)
iq = np.where(trel <= 10.0)[0]
gt49x = gt49.copy()
gt49x[iq] = gt49x[iq]/

fig = plt.figure(figsize=(12,8))
ax = fig.subplots()
ax.plot(trel,gt48,'ro',label='EPN_CTI_v0048')
ax.plot(trel,gt49,'b^',label='EPN_CTI_v0049')
ax.set_ylabel("g(t)")
ax.set_xlabel(r"Years since 2000-01-01T00:00:00")
ax.grid(True)
ax.legend()

#%%
#
# now update

#hdu48['LONG_TERM_CTI'].data['T_COEFF'][ix48[1]] = gt49
#hdu48.writeto("/xdata/xcaldata/XMM/IVAN/PN_SW/ccfdev/EPN_CTI_0049.CCF_test25",overwrite=True)


