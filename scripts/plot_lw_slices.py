#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 25 11:44:44 2019

plot the individual spectra for RAWY slices

@author: ivaltchanov
"""
import numpy as np
import os

from astropy.io import fits
from astropy.convolution import convolve, Box1DKernel

import matplotlib.pylab as plt
import seaborn as sns
sns.set(style="white")

plt.rc('text', usetex=False)
plt.rc('font', family='serif')
#%%

wdir = "/home/ivaltchanov/IVAN/PN_LW"
#obsid = "0727360201"
obsid = "0203720201"
expo = "S020"

specDir = f"{wdir}/{obsid}/PN_LW"
smooth = 7

fig, axs = plt.subplots(2,6,sharex=True,sharey=True,figsize=(15,10))
for j in np.arange(12):
    ccd = j+1
    #slices = glob.glob(f"{specDir}/pn_{ccd:02}_180_spec5.fits")
    slices = f"{specDir}/pn_{expo}_{ccd:02}_all_spec5.fits"
    if (not os.path.isfile(slices)):
        raise FileNotFoundError
    hdu = fits.open(slices)
    spec = hdu['SPECTRUM']
    channel = spec.data['CHANNEL']*5.0/1000.0
    counts = spec.data['COUNTS']
    y = convolve(counts, Box1DKernel(smooth))
    k = j
    kj = 0
    if (j >= 6):
        k = j - 6
        kj = 1
    axs[kj,k].plot(channel,counts,label=f'CCD #{ccd:02}')
    axs[kj,k].plot(channel,y,label='')
    axs[kj,k].set_xlim((1,9))
    axs[kj,k].set_ylim((0,300))
    axs[kj,k].grid(True)
    axs[kj,k].legend()
    #axs[kj,k].set_title(f"CCD #{ccd}")
#
plt.subplots_adjust(wspace=0, hspace=0)
plt.text(-13,-1,'Energy (keV)',ha='center', va='center')
plt.text(-36,10,'Counts',rotation='vertical',ha='center', va='center')
#plt.tight_layout()
#plt.savefig(f"{wdir}/{obsid}_{expo}_allccd_plot.png",dpi=100)
plt.show()
plt.close()
