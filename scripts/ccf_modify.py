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

wdir = "/home/ivaltchanov/IVAN/PN_SW/ccfdev"

os.chdir(wdir)
ccf49_file = 'EPN_CTI_0049.CCF_test25'
#
ccf48_file = '/xdata/ccf/valid/EPN_CTI_0048.CCF'

hdu48 = fits.open(f"{ccf48_file}")
ltc48 = hdu48['LONG_TERM_CTI']
ix48 = np.where(ltc48.data['MODE_ID'] == 3)[0]

with fits.open(f"{wdir}/{ccf49_file}", mode='update') as hdu49:
    ix49 = np.where(hdu49['LONG_TERM_CTI'].data['MODE_ID'] == 3)[0]
    tcoeff = hdu49['LONG_TERM_CTI'].data['T_COEFF']
    tcoeff[ix49[1]][0] = a0
    #tcoeff[ix49[2]][0] = 0.0004531 # taken from the spline at t=0
    tcoeff[ix49[2]][0] = a0x
    jjj = tcoeff[ix49][2][10:]
    jjj += 1.0e-5
    tcoeff[ix49[2]][10:] = jjj
    hdu49['LONG_TERM_CTI'].data['T_COEFF'] = tcoeff
    #
    print (hdu49['LONG_TERM_CTI'].data['T_COEFF'][ix49[1]][0])
    print (hdu49['LONG_TERM_CTI'].data['T_COEFF'][ix49][2][0])
    print (hdu49['LONG_TERM_CTI'].data['T_COEFF'][ix49][2])
    hdu49.flush()  # changes are written back to original.fits
#
#
