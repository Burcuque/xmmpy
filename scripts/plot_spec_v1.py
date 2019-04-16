#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 12 15:43:05 2018

Compare spectra from the same physical region in PN FF and SW

@author: ivaltchanov
"""
import os
from astropy.io import fits
from astropy.table import Table

import matplotlib.pylab as plt

home = os.path.expanduser('~')
wdir = home + '/XMM/CAL/ODF'
#
obs_ff = '0120300501'
obs_sw = '0120300801'
#
file1 = '{}/{}/pn_ff_regionX.fits'.format(wdir,obs_ff)
file2 = '{}/{}/pn_sw_regionX.fits'.format(wdir,obs_sw)
#
#%%
# 
hdu = fits.open(file1)
chan1 = hdu['SPECTRUM'].data["CHANNEL"]
peak1 = hdu['SPECTRUM'].data["COUNTS"]
expo1 = hdu['SPECTRUM'].header['EXPOSURE']
hdu = fits.open(file2)
chan2 = hdu['SPECTRUM'].data["CHANNEL"]
peak2 = hdu['SPECTRUM'].data["COUNTS"]
expo2 = hdu['SPECTRUM'].header['EXPOSURE']
#
#%%
plt.plot(chan1,peak1/expo1,label='{}: PN FF calclosed'.format(obs_ff))
plt.plot(chan2,peak2/expo2,label='{}: PN SW calclosed'.format(obs_sw))
plt.xlim([0,300])
plt.xlabel("Channel")
plt.ylabel("Count rate")
plt.grid(True)
plt.legend()
#plt.show()
plt.savefig('{}/{}_vs_{}_spec.png'.format(wdir,obs_ff,obs_sw),dpi=100)
plt.close()
