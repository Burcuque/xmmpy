#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov  9 10:36:31 2018

Modify CCF file "by hand "

@author: ivaltchanov
"""

from astropy.io import fits
import numpy as np

import matplotlib.pylab as plt

wdir = "/home/ivaltchanov/IVAN/PN_SW/ccfdev"

ccf49_file = 'EPN_CTI_0049.CCF_orig'
ccf49x_file = 'EPN_CTI_0049.CCF'
ccf49xx_file = 'EPN_CTI_0049.CCF_test2'
#
ccf48_file = '/xdata/ccf/valid/EPN_CTI_0048.CCF'

hdu48 = fits.open(f"{ccf48_file}")
hdu49 = fits.open(f"{wdir}/{ccf49_file}")
hdu49x = fits.open(f"{wdir}/{ccf49x_file}")
hdu49xx = fits.open(f"{wdir}/{ccf49xx_file}")
#
ltc48 = hdu48['LONG_TERM_CTI']
ix48 = np.where(ltc48.data['MODE_ID'] == 3)[0]
ltc49 = hdu49['LONG_TERM_CTI']
ix49 = np.where(ltc49.data['MODE_ID'] == 3)[0]
ltc49x = hdu49x['LONG_TERM_CTI']
ix49x = np.where(ltc49x.data['MODE_ID'] == 3)[0]
ltc49xx = hdu49xx['LONG_TERM_CTI']
ix49xx = np.where(ltc49xx.data['MODE_ID'] == 3)[0]
#
xtime = hdu49['LTC_TIMES'].data['TIME'][0]
y48 = ltc48.data['T_COEFF'][ix48[1]]
y49 = ltc49.data['T_COEFF'][ix49[1]]
y49x = ltc49x.data['T_COEFF'][ix49x[1]]
y49a = ltc49.data['T_COEFF'][ix49[2]]
#y49b = ltc49x.data['T_COEFF'][ix49x[2]]
y49c = ltc49xx.data['T_COEFF'][ix49xx[2]]
#
fig = plt.figure()
ax = fig.subplots()

ax.plot(xtime, y48 ,label='CTIv48')
#ax.plot(xtime, y49,'or',label='CTIv49')
ax.plot(xtime, y49,label='CTIv49 orig')
ax.plot(xtime, y49a,'+',label='CTIv49 test0')
ax.plot(xtime, y49c,'o-',label='CTIv49 test2')
#ax.plot(xtime, y49c,label='CTIv49 test2')
#
ax.grid(True)
ax.set_ylabel("T_COEFF")
ax.set_xlabel("Year since 2000-01-01T00:00:00")
ax.legend()
#
#
 