#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov  9 10:36:31 2018

Modify CCF file:
    
    will use CCF v48 as baseline
    will replace MODE_ID=3 second entry with the newly derived on efor CalClosed

@author: ivaltchanov
"""
import os

from astropy.io import fits
import numpy as np

wdir = "/home/ivaltchanov/IVAN/PN_SW/ccfdev"
os.chdir(wdir)
#%%
ccf49_file = 'EPN_CTI_0049.CCF'
#
ccf48_file = '/xdata/ccf/valid/EPN_CTI_0048.CCF'

hdu48 = fits.open(f"{ccf48_file}")
ltc48 = hdu48['LONG_TERM_CTI']
ix48 = np.where(ltc48.data['MODE_ID'] == 3)[0]
tcoeff = hdu48['LONG_TERM_CTI'].data['T_COEFF']
aaa = tcoeff[ix48[1]]
#%%
hdu49 = fits.open(f"{ccf49_file}")
ltc49 = hdu49['LONG_TERM_CTI']
ix49 = np.where(ltc49.data['MODE_ID'] == 3)[0]
xtcoeff = hdu49['LONG_TERM_CTI'].data['T_COEFF']
bbb = xtcoeff[ix49[1]]
#
a0 = 4.22e-4 # empirical
print ("a0 = ", a0)
bbb[0] = a0
hdu48['LONG_TERM_CTI'].data['ENERGY'][ix48[1]] = 5.8967
hdu48['LONG_TERM_CTI'].data['T_COEFF'][ix48[1]] = bbb
#%%
hdu48.writeto("EPN_CTI_0049_test0.CCF",overwrite=True)
#
