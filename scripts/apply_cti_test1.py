#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 22 12:30:38 2018

PN long-term CTI calibration test

Will read the no_cti results on selected targets and apply the new CTI correction


@author: ivaltchanov
"""

#import os
import numpy as np
#import pandas as pd

import matplotlib.pylab as plt

#from astropy.io import fits
from astropy.table import Table

from scipy.interpolate import interp1d
#from scipy.interpolate import UnivariateSpline

#%%
#
def get_ltci_corr(ccf_file,relative_time,index=1,rawy=190,verbose=False):
    #
    # red the CTI CCF file, and return the correction to apply to E_obs/E_lab
    #
    # relative_time is in years since 2000-01-01T00:00:00 the output will be calculated on this grid
    #
    # index=1 is for 5.899 keV
    # index=2 is for 6.4 keV (new CCF only)
    # rawy is the average RAWY, default to boresight of 190
    #
    ltc = Table.read(ccf_file,hdu='LONG_TERM_CTI')
    xt = Table.read(ccf_file,hdu='LTC_TIMES')
    xtime = xt[0]['TIME'] # year since 01-01-2000T00:00:00
    ix = np.where(ltc['MODE_ID'] == 3)[0]
    tx = ltc[ix]
    if (len(tx) < 3 and index>=2):
        print (f"CCF file has no 3 entries for MODE_ID=3 (PN_SW) and index is {index}")
        return None
    tcoef = tx[index]['T_COEFF']
    if (verbose): 
        print (f'{ccf_file}, index={index}',tcoef)
    g = np.power((1.0 - tcoef)/(1.0 - tcoef[0]),rawy)
    #
    inter = interp1d(xtime,g)
    # and the correted observed ratios
    corr_out = inter(relative_time)
    return corr_out

#%%
#
# First we'll read the no_cti Fe K line fit results 
#
wdir = '/xdata/xcaldata/XMM/IVAN/PN_SW'
#
redshift = {'ngc4151': 0.003262, 'ngc3227': 0.00386, 'mrk1048': 0.0427, 'ngc3783': 0.009755,\
            'ngc4593': 0.008344, 'ngc5506': 0.00589, 'mcg-5-23-16': 0.008226, 'ngc3516': 0.008816,\
            'ngc5548': 0.01627, 'ngc2992':  0.007296, 'ngc1566': 0.005036}
feK = 6.399 # the weitghed mean of Kalpha_2 at 6.3908 (intensity 50) and Kalpha_1 at 6.40308 (intensity 100)
#
t = Table.read(f"{wdir}/selected_targets_nocti.csv",comment="\s*#")
targets = t['target']
nt = len(targets)

ratio = []
ratio_err = []
rtime = []
#
for i in np.arange(nt):
    target = targets[i]
    #
    lineX =  feK/(1.0 + redshift[target]) # redshifted line position
    line = t['lineE'].data[i]
    lineErrUp = t['lineE_low'].data[i]
    lineErrLo = t['lineE_up'].data[i]
    lineErr = (lineErrUp + lineErrLo)/2.0
    rx = line/lineX
    rxerr = lineErr/lineX
    # build up the arrays
    rtime = np.concatenate((rtime,t['delta_time'].data[i]),axis=None)
    ratio = np.concatenate((ratio,rx),axis=None)
    ratio_err = np.concatenate((ratio_err,rxerr),axis=None)
#
# sort on rtime for the interpolation
ix = np.argsort(rtime)
rtime = rtime[ix]
ratio = ratio[ix]
ratio_err = ratio_err[ix]
#%%
#
# Now the calClosed data
#
tcc = Table.read(f"{wdir}/../PN_calclosed/pn_sw_calclosed_fit_results_nocti_sorted.csv",
               comment= "\s*#")
cc_rtime = tcc['time'].data
#
MnKa = (0.162*5.888 + 0.6*5.899)/0.762 # keV probabilities from Wikipedia on Iron-55
line2 = tcc['mn1_cen'].data/1000.0
line2Err = tcc['mn1_cen_err'].data/1000.0
cc_ratio = line2/MnKa
cc_ratio_err = line2Err/MnKa

#%%
cal48 = '/xdata/ccf/pub/EPN_CTI_0048.CCF'
cal49a = '/xdata/xcaldata/XMM/IVAN/PN_SW/ccfdev/EPN_CTI_0049.CCF_orig'
cal49b = '/xdata/xcaldata/XMM/IVAN/PN_SW/ccfdev/EPN_CTI_0049.CCF_test3'

runyear = np.linspace(0.0,25.0,num=26)

corr48cc = get_ltci_corr(cal48,cc_rtime,index=1)
corr48 = get_ltci_corr(cal48,rtime)
xcorr48 = get_ltci_corr(cal48,runyear)
# CalClosed
cc_ratio48 = cc_ratio/corr48cc 
cc_ratio48err = cc_ratio_err/corr48cc
mnK48 = MnKa*(1-cc_ratio48)*1000 # in eV
mnK48err = MnKa*cc_ratio48err*1000
# AGNs
ratio48 = ratio/corr48 
ratio48err = ratio_err/corr48 
fe48 = feK*(1-ratio48)*1000 # in eV
fe48err = feK*ratio48err*1000
#
corr49acc = get_ltci_corr(cal49a,cc_rtime,index=1)
corr49a = get_ltci_corr(cal49a,rtime,index=2)
xcorr49a = get_ltci_corr(cal49a,runyear,index=2)
# CalClosed
cc_ratio49a = cc_ratio/corr49acc 
cc_ratio49a_err = cc_ratio_err/corr49acc
mnK49a = MnKa*(1-cc_ratio49a)*1000 # in eV
mnK49aerr = MnKa*cc_ratio49a_err*1000
# AGNs
ratio49a = ratio/corr49a 
ratio49a_err = ratio_err/corr49a 
fe49a = feK*(1-ratio49a)*1000
fe49a_err = feK*ratio49a_err*1000
#
corr49bcc = get_ltci_corr(cal49b,cc_rtime,verbose=True,index=1,rawy=170)
corr49b = get_ltci_corr(cal49b,rtime,index=2,verbose=True)
xcorr49b = get_ltci_corr(cal49b,runyear,index=2)
# CalClosed
cc_ratio49b = cc_ratio/corr49bcc 
cc_ratio49b_err = cc_ratio_err/corr49bcc
mnK49b = MnKa*(1-cc_ratio49b)*1000 # in eV
mnK49berr = MnKa*cc_ratio49b_err*1000
# AGNs
ratio49b = ratio/corr49b
ratio49b_err = ratio_err/corr49b 
fe49b = feK*(1-ratio49b)*1000
fe49b_err = feK*ratio49b_err*1000
#
#%%
fig = plt.figure(figsize=(12,8))
ax = fig.subplots(2,1,sharex=True)
fig.subplots_adjust(hspace=0)
#
# ratios and corrections
#
ax[0].errorbar(rtime, ratio, yerr=ratio_err,fmt='or',alpha=0.4,label='AGNs')
ax[0].errorbar(cc_rtime, cc_ratio, yerr=cc_ratio_err,fmt='^k',alpha=0.4,label='CalClosed')
ax[0].plot(runyear, xcorr48,'r-',label='v48')
ax[0].plot(runyear, xcorr49a,'g-',label='v49a')
ax[0].plot(runyear, xcorr49b,'b-',label='v49b')
ax[0].axhline(1.0,color='k',ls='dashed')
ax[0].grid(True)
ax[0].set_ylabel("Observed, E$_{obs}$/E$_{lab}$")
ax[0].legend()
#
#ax[1].errorbar(rtime, ratio48, yerr=ratio48err,fmt='or',alpha=0.4,label='v48')
#ax[1].errorbar(rtime, ratio49a, yerr=ratio49a_err,fmt='og',alpha=0.4,label='v49a')
#ax[1].errorbar(rtime, ratio49b, yerr=ratio49b_err,fmt='ob',alpha=0.4,label='v49b')
ax[1].errorbar(rtime, fe48, yerr=fe48err,fmt='or',alpha=0.4,label='v48')
ax[1].errorbar(rtime, fe49a, yerr=fe49a_err,fmt='og',alpha=0.4,label='v49a')
ax[1].errorbar(rtime, fe49b, yerr=fe49b_err,fmt='ob',alpha=0.4,label='v49b')
#
ax[1].errorbar(cc_rtime, mnK48, yerr=mnK48err,fmt='^r',alpha=0.4,label='v48')
ax[1].errorbar(cc_rtime, mnK49a, yerr=mnK49aerr,fmt='^g',alpha=0.4,markersize=8,label='v49a')
ax[1].errorbar(cc_rtime, mnK49b, yerr=mnK49berr,fmt='^b-',alpha=0.4,label='v49b')
ax[1].grid(True)
ax[1].set_ylabel("Corrected, E$_{obs}$-E$_{lab}$ (eV)")
ax[1].set_xlabel("Years since 2000-01-01T00:00:00")
ax[1].legend()
#
plt.show()
