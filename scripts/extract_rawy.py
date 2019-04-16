#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 23 10:50:09 2018

Extract the RAWY coordinate for the PN SW AGN observations

The RAWY is stored in each pn_src_spectrum.fits primary header keyword SRCPOSY

@author: ivaltchanov
"""

import os
import numpy as np

from astropy.table import Table, Column
from astropy.io import fits

import matplotlib.pylab as plt
#from matplotlib.patches import Patch

targets = ["ngc5506","ngc5548", "ngc4151","ngc3516","ngc2992","ngc3227",'ngc3783','ngc1566']
#
wdir = "/xdata/xcaldata/XMM/IVAN/PN_SW"
os.chdir(wdir)

#%%
for target in targets:
    tabfile = f'{target}/output_xspec_nocti.csv'
    t = Table.read(tabfile)
    nt = len(t)
    # get the FITS header keyword SRCPOSY
    #
    obsids = t['obsid']
    rawy = []
    for obsid in obsids:
        specfile = f'{target}/{obsid:010}/no_cti/pn_src_spectrum_nocti.fits'
        if (os.path.isfile(specfile)):
            hdu = fits.open(specfile)
            rawy.append(hdu[0].header['SRCPOSY'])
        else:
            print (f"{specfile} not found.")
            raise FileNotFoundError
    #
    t.add_column(Column(rawy),index=7,name='rawy')
    t.write(tabfile,overwrite=True)
    print (f"Added RAWY column for {target}")
#
#%%
#
# now add RAWY for the CalClosed
#
wdir = "/xdata/xcaldata/XMM/IVAN/PN_calclosed"
os.chdir(wdir)
tabfile = "pn_sw_calclosed_fit_results_nocti.csv"
t = Table.read(tabfile)
nt = len(t)
obsids = t['obsid']
rawy = []
for obsid in obsids:
    specfile = f'{obsid:010}/no_cti/spectrum_bin5.fits'
    if (os.path.isfile(specfile)):
        hdu = fits.open(specfile)
        rawy.append(hdu[0].header['SRCPOSY'])
    else:
        print (f"{specfile} not found.")
        raise FileNotFoundError
#
t.add_column(Column(rawy),index=4,name='rawy')
t.write(tabfile,overwrite=True)
print (f"Added RAWY column for all calClosed in PN SW")
#
#%%
#
# Now check the average RAWY for all AGNs
for target in targets:
    tabfile = f'{target}/output_xspec_nocti.csv'
    t = Table.read(tabfile)
    nt = len(t)
    # get the FITS header keyword SRCPOSY
    #
    print ("{}: mean RAWY={}, st.dev.={}".format(target,np.mean(t['rawy']),np.std(t['rawy'])))
#
#%%
