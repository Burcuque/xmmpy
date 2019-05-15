#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 10 10:16:15 2019

Filter the XMM spectra to only keep those with more than XXX counts in [4,8] keV

For Fe Kalpha (6.4 keV) classification

@author: ivaltchanov
"""
import os
import glob
import numpy as np
import pandas as pd
from shutil import copy

from astropy.io import fits
from astropy.stats import histogram

import matplotlib.pylab as plt

#%%
wdir = '/lhome/ivaltchanov/XMM-data/ML_spectra'

df = pd.read_csv(f'{wdir}/xmm_spectra_for_ml.csv')

files = glob.glob(f'{wdir}/spectra01/*.fits')

print (f"Found {len(files)} spectra")
#%%
counts4_8 = []
os.chdir(f'{wdir}/spectra01')
for ifile in df['filename']:
    if (not os.path.isfile(ifile)):
        counts4_8.append(0)
        continue
    #if ('0651850501' not in ifile):
    #    continue
    hdu = fits.open(ifile)
    spec = hdu['SPECTRUM']
    binsize = spec.header['SPECDELT']
    x = spec.data['CHANNEL']*binsize/1000.0
    y = spec.data['COUNTS']
    ix = np.where((x<=8.0) & (x>=4.0))[0]
    tot = np.sum(y[ix])
    counts4_8.append(tot)
    print (f'Total counts in spectrum {os.path.basename(ifile)} in [4,8] keV: {tot}')
#
#%%
df['counts4_8'] = pd.Series(counts4_8)
df.to_csv(f'{wdir}/xmm_spectra_for_ml_v2.csv')
#%%
a = histogram(df.counts4_8,bins='blocks',range=(0,5000))
fig = plt.figure(figsize=(15,5))
ax = fig.subplots()
#ax.semilogx(x,y)
ax.step(a[1][1:],a[0],color='red')
#ax.step(a[0],a[1][1:],color='red')
ax.set_xlabel("Counts")
ax.set_ylabel("Frequency")
ax.grid()
#
#
#%%
#
# now select those with more than 1000 counts in [4,8] keV
#
out_dir = f'{wdir}/spectra01_select'
for i,qfile in enumerate(df.filename):
    if (df.counts4_8.iloc[i] > 1000):
        basename = os.path.basename(qfile).split('.')[0]
        pngfile = f'images/{basename}.png'
        if (os.path.isfile(qfile)):
            copy(qfile,out_dir)
        if (os.path.isfile(pngfile)):
            copy(pngfile,f'{out_dir}/images')
    pass
#
        
    



