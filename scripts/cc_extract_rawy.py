#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 23 10:50:09 2018

Extract the RAWY coordinate for the PN SW AGN observations

The RAWY is stored in each pn_src_spectrum.fits primary header keyword SRCPOSY

@author: ivaltchanov
"""

import os
import glob
import numpy as np

from astropy.table import Table, Column
from astropy.io import fits

import matplotlib.pylab as plt
#from matplotlib.patches import Patch

#%%
# now add RAWY for the CalClosed
#
wdir = "/xdata/xcaldata/XMM/IVAN/PN_calclosed"
os.chdir(wdir)
tabfile = "pn_sw_calclosed_fit_results_nocti.csv"
t = Table.read(tabfile)
nt = len(t)
obsids = t['obsid']

xc = 517
yc = 1430
xbox = 2500
ybox = 2250
#
all_wrawy = []
all_rawy = []
for i,obsid in enumerate(obsids):
    #specfile = '{:010}/no_cti/spectrum_bin5.fits'.format(int(obsid))
    evlist = glob.glob('{}/{:010}/no_cti/*PIEVLI*'.format(wdir,int(obsid)))
    hdu = fits.open(evlist[0])
    events = hdu['EVENTS']
    m1 = np.abs((events.data['DETX'] - xc)) <= xbox
    m2 = np.abs((events.data['DETY'] - yc)) <= ybox
    ix = np.where(m1*m2)[0]
    print ("Found {} events".format(len(ix)))
    xev = events.data[ix]
    #
    #what = "DET"
    #what = "RAW"
    #fig = plt.figure(figsize=(12,8))
    #ax = fig.subplots()
    #ax.plot(events.data[f'{what}X'],events.data[f'{what}Y'],'b.')
    #ax.plot(xev[f'{what}X'],xev[f'{what}Y'],'r.')
    #ax.set_xlabel(f"{what}X")
    #ax.set_ylabel(f"{what}Y")
    #
    # now weighted mean of RAWY, weighted by energy in [5.5-6.5] keV
    #
    mq1 = xev["PATTERN"] == 0
    mq2 = xev["FLAG"] == 0
    mq3 = xev["PAT_SEQ"] == 0
    # elect only Mn Kalpha at 5.89 keV
    mq4 = xev["PI"] <= 6200.0
    mq5 = xev["PI"] >= 5600.0
    iq = np.where(mq1 * mq2 * mq3 * mq4 *mq5)[0]
    #
    #ax.hist(xev['PI'][iq]/1000.0,bins=100)
    w = xev["PI"][iq]
    wrawy = np.average(xev["RAWY"][iq],weights=w)
    print (f"{obsid}: uweighted mean RAWY = {wrawy}")
    #
    rawy = np.average(xev["RAWY"][iq])
    print (f"{obsid}: unweighted mean RAWY = {rawy}")
#    t["rawy"][i] = rawy
    all_rawy.append(rawy)
    all_wrawy.append(wrawy)
    #ax.set_xlim([0.0,10.0])
    #plt.show()
#
#t.add_column(Column(rawy),index=4,name='rawy')
#tabfile_new = "pn_sw_calclosed_fit_results_nocti_rawy.csv"
#t.write(tabfile_new,overwrite=True)
#print (f"Added RAWY column for all calClosed in PN SW")
#
print ("Mean PI-weighted RAWY for all CalClosed: {} +/- {}".format(np.average(all_wrawy),np.std(all_wrawy)))
print ("Mean unweighted RAWY for all CalClosed: {} +/- {}".format(np.average(all_rawy),np.std(all_rawy)))
