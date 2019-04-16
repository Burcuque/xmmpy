#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 22 12:30:38 2018

Understanding the PN long-term CTI calibration file

work1 is for experiments
work2 is for the selected method (UnivariateSpline)

@author: ivaltchanov
"""

#import os
import numpy as np
#import pandas as pd

import matplotlib.pylab as plt

#from astropy.io import fits
from astropy.table import Table

from scipy.interpolate import interp1d
from scipy.interpolate import UnivariateSpline

#%%
# read the CAL file with the long-term CTI correction
#
wdir = '/xdata/xcaldata/XMM/IVAN'
#
calNew = 'EPN_CTI_0048.CCF'
calfile = f'/xdata/ccf/pub/{calNew}'
ltc = Table.read(calfile,hdu='LONG_TERM_CTI')
ltc.info()
xt = Table.read(calfile,hdu='LTC_TIMES')
xt.info()
xtime = xt[0]['TIME'] # year since 01-01-2000T00:00:00
#
#%%
#
# PN_SW is with MODE_ID=3
#
ix = np.where(ltc['MODE_ID'] == 3)[0]
tx = ltc[ix]
#
#%%
#
# now get the current Longe-term CTI curves for Mn Kalpha
#
#idx_line = 0 # Al Kaplpha at 1.4 keV
idx_line = 1 # Mn Kalpha at 5.8988 keV
ltci = tx[idx_line]['T_COEFF']
#ltci2 = tx2[idx_line]['T_COEFF']
#
rawy = 170 # derived as the weighted mean for each  CalClosed observation
#
# now the magic to convert to ratio
xcorr137 = np.power((1.0 - ltci)/(1.0 - ltci[0]),137.0)
xcorr170 = np.power((1.0 - ltci)/(1.0 - ltci[0]),rawy)
xcorr200 = np.power((1.0 - ltci)/(1.0 - ltci[0]),200.0)
# no interpolate, to be used later
#
f137 = interp1d(xtime, xcorr137)
f200 = interp1d(xtime, xcorr200)
f170 = interp1d(xtime, xcorr170)
#
#%%
#
# Now read the PN-SW Calclosed results with no long-term CTI correction applied
# the table is sorted on time, a couple of deviating results were also commented out
# 
t = Table.read(f"{wdir}/PN_calclosed/pn_sw_calclosed_fit_results_nocti_sorted.csv",
               comment= "\s*#")
nt = len(t)
#
rev = t['rev'].data
rtime = t['time'].data
#
line2_lab = (0.162*5.888 + 0.6*5.899)/0.762 #eV probabilities from Wikipedia on Iron-55
line2 = t['mn1_cen'].data/1000.0
line2Err = t['mn1_cen_err'].data/1000.0
ratio = line2/line2_lab
ratio_err = line2Err/line2_lab
# set the weights
weight = 1.0/line2Err**2
#%%
#
# this is the grid used in the CCF for the year
#
runyear = np.linspace(0.0,25.0,num=26)
#
# interpolate the current CCF on the grid
#
# now the magic to convert to ratio
ccf_int137 = f137(runyear)
ccf_int200 = f200(runyear)
ccf_int170 = f170(runyear)
#
#%%
#
# now use UnivariateSpline 
#
# latest pointw will have the max weight
i25 = np.where(rtime >= 2.5)[0]
weight[i25] = np.max(weight)
# normalize the weights
weight = weight/np.max(weight)
#
# the CTI at t=0
#a0 = 0.00043236
a0 = 0.0
#
iend = len(rtime) - 1
#
rtime0 = np.insert(rtime,iend+1,[20.0,22.0,24.0])
ratio0 = np.insert(ratio,iend+1,[ratio[iend]*0.998,ratio[iend]*0.998,ratio[iend]*0.998])
wxx = np.insert(weight,iend+1,[0.5,0.1,0.1])
# normalize the weights
wxx = wxx/np.max(wxx)
#
yy = 1.0 - (1.0-a0)*np.power(ratio,1.0/rawy)
yy0 = 1.0 - (1.0-a0)*np.power(ratio0,1.0/rawy)
#rtime0 = np.insert(rtime,0,0.0)
#yy0 = np.insert(yy,0,a0)


s = UnivariateSpline(rtime0, yy0, w=wxx,k=5,s=1)
#s = UnivariateSpline(rtime, yy, w=weight,s=1)
result = s(runyear)
#
# now everything after runyer = 20 is set to runyer=20
# runyer=20 is index 
result[21:] = result[20]
#result[0] = a0

#yn1 = np.power((1.0 - s(runyear))/(1.0-a0),190.0)
#yn200 = np.power((1.0 - s(runyear))/(1.0-a0),200.0)
#yn137 = np.power((1.0 - s(runyear))/(1.0-a0),137.0)
#
yn170 = np.power((1.0 - result)/(1.0-a0),rawy)
yn200 = np.power((1.0 - result)/(1.0-a0),200.0)
yn137 = np.power((1.0 - result)/(1.0-a0),137.0)

fig = plt.figure(figsize=(12,8))
ax = fig.subplots(2,1,sharex=True)
fig.subplots_adjust(hspace=0)

ax[0].errorbar(rtime, ratio, yerr=ratio_err,fmt='or',alpha=0.4,label='PN SW calclosed')
ax[0].plot(runyear, yn137,'--',label='RAWY=137')
ax[0].plot(runyear, yn170,color='black',linewidth=2.0,label=f'RAWY={rawy}')
ax[0].plot(runyear, yn200,'--',label='RAWY=200')
#
# the current CCF
#
ax[0].plot(runyear, ccf_int170,color='blue',label=f'CCFv48 at RAWY={rawy}')
#
#ax.set_xlim([0.6,1.0])
ax[0].set_ylim([0.96,1.0])
ax[0].grid(True)
ax[0].set_ylabel("Mn K$_\\alpha$, E$_{obs}$/E$_{lab}$")
ax[0].set_title("PN SmallWindow Energy Scale Analysis\n Long-term CTI at 5.8967 keV")
ax[0].legend()
#
# let's check the actual tables that wil go in the CCF
#
ax[1].plot(xtime,ltci*1.0e4,color='blue',label=f"{calNew}")
ax[1].plot(runyear,result*1.0e4,color='black',linewidth=2.0,label="Updated")
ax[1].set_xlabel("Years since 2000-01-01T00:00:00")
ax[1].set_ylabel(r"T_COEFF $\times 10^4$")
ax[1].grid(True)
ax[1].legend()

plt.savefig(f'{wdir}/pn_sw_ltci_ccf49.png',dpi=100)
plt.show()
plt.close()
#%%
# Save the new table to a file 
#
with open(f'{wdir}/PN_SW/ccfdev/ccf49_ltcti_table_mode3.txt',"w") as fout:
    print ("{:3} : {:3} :   {:6.4f} : ".format(3,4,5.8988),end=" ",file=fout)
    for j in result:
        print (" {:3.6e}".format(j),end=" ",file=fout)
    print (":   1",file=fout)
#
#%%
#
# now update the FITS file
#
from astropy.io import fits

calNew = 'EPN_CTI_0048.CCF'
calfile = f'/xdata/ccf/pub/{calNew}'
hdus = fits.open(calfile)
#
hdus[0].header['ISSUE'] = 49
ix = np.where(hdus['LONG_TERM_CTI'].data["MODE_ID"] == 3)[0]
hdus['LONG_TERM_CTI'].data["T_COEFF"][ix[1]] = result
hdus.writeto(f"{wdir}/EPN_CTI_0049_test.CCF",overwrite=True)

