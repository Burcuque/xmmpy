#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 22 12:30:38 2018

Understanding the PN long-term CTI calibration file

work1 is for experiments
work2 is for the selected method (UnivariateSpline) for CalClosed
work3 is for the AGNs
work4 is for the AGNs but only for the selected targets

@author: ivaltchanov
"""

import os
import numpy as np
#import pandas as pd

import matplotlib.pylab as plt

from astropy.io import fits
from astropy.table import Table

from scipy.interpolate import interp1d
from scipy.interpolate import UnivariateSpline

#%%
# read the CAL file with the long-term CTI correction
#
wdir = '/xdata/xcaldata/XMM/IVAN/PN_SW'
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
# now the magic to convert to ratio
xcorr137 = np.power((1.0 - ltci)/(1.0 - ltci[0]),137.0)
xcorr190 = np.power((1.0 - ltci)/(1.0 - ltci[0]),190.0)
xcorr200 = np.power((1.0 - ltci)/(1.0 - ltci[0]),200.0)
# no interpolate, to be used later
#
f137 = interp1d(xtime, xcorr137)
f200 = interp1d(xtime, xcorr200)
f190 = interp1d(xtime, xcorr190)
#
#%%
#
# Now read the PN-SW results for the sources with no long-term CTI correction applied
#

os.chdir(wdir)

#targets = ["ngc5506","ngc5548", "ngc4151","ngc3516","ngc2992","ngc3227",'ngc3783']
#
# redshifts
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
#
#%%
# now sort on rtime
# 
ix = np.argsort(rtime)
rtime = rtime[ix]
ratio = ratio[ix]
ratio_err = ratio_err[ix]
# set the weights
weight = 1.0/ratio_err**2
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
ccf_int190 = f190(runyear)
#
#%%
#
# now use UnivariateSpline 
#
# latest pointw will have the max weight
#m1 = rtime < 12
#m2 = ratio > 0.985
#iqq = np.where(m1*m2)[0]
#weight[i25] = np.max(weight)
# normalize the weights
weight = weight/np.max(weight)
#
# the CTI at t=0
a0 = 0.00043236
#
iend = len(rtime) - 1
#
rtime0 = np.insert(rtime,iend+1,[21.0,22.0])
ratio0 = np.insert(ratio,iend+1,[0.9735,0.972])
wxx = np.insert(weight,iend+1,[0.5,0.1])
rtime0 = np.insert(rtime0,0,0.0)
ratio0 = np.insert(ratio0,0,1.0)
wxx = np.insert(wxx,0,0.5)
# normalize the weights
#wxx[iqq] = 0.01
#wxx = wxx/np.max(wxx)
#
yy = 1.0 - (1.0-a0)*np.power(ratio,1.0/190.0)
yy0 = 1.0 - (1.0-a0)*np.power(ratio0,1.0/190.0)


s = UnivariateSpline(rtime0, yy0, w=wxx,k=5,s=1)
#s = UnivariateSpline(rtime, yy, w=weight,k=5,s=1)
result = s(runyear)
#
# now everything after runyer = 20 is set to runyer=20
# runyer=20 is index 
result[17:] = result[16]
result[0] = a0

#yn1 = np.power((1.0 - s(runyear))/(1.0-a0),190.0)
#yn200 = np.power((1.0 - s(runyear))/(1.0-a0),200.0)
#yn137 = np.power((1.0 - s(runyear))/(1.0-a0),137.0)
#
yn190 = np.power((1.0 - result)/(1.0-a0),190.0)
yn200 = np.power((1.0 - result)/(1.0-a0),200.0)
yn137 = np.power((1.0 - result)/(1.0-a0),137.0)

fig = plt.figure(figsize=(12,8))
ax = fig.subplots(2,1,sharex=True)
fig.subplots_adjust(hspace=0)

ax[0].errorbar(rtime, ratio, yerr=ratio_err,fmt='or',alpha=0.4,label='PN SW AGN')
ax[0].plot(runyear, yn137,'--',label='RAWY=137')
ax[0].plot(runyear, yn190,color='black',linewidth=2.0,label='RAWY=190')
ax[0].plot(runyear, yn200,'--',label='RAWY=200')
#
# the current CCF
#
ax[0].plot(runyear, ccf_int200,color='blue',label='CCFv48 at RAWY=190')
#
#ax.set_xlim([0.6,1.0])
ax[0].set_ylim([0.96,1.0])
ax[0].grid(True)
ax[0].set_ylabel("Fe K$_\\alpha$, E$_{obs}$/E$_{lab}$")
ax[0].set_title("PN SmallWindow Energy Scale Analysis\n Long-term CTI at 6.399 keV")
ax[0].legend()
#
# let's check the actual tables that wil go in the CCF
#
ax[1].plot(xtime,ltci*1.0e4,color='blue',label=f"{calNew}")
ax[1].plot(runyear,result*1.0e4,color='black',linewidth=2.0,label="Updated")
ax[1].set_xlabel("Years since 2000-01-01T00:00:00")
ax[1].set_ylabel(r"T_COEFF $\times 10^4$")
#
# now the CalClosed
# 
hdu = fits.open('ccfdev/EPN_CTI_0049.CCF')
ix = np.where(hdu['LONG_TERM_CTI'].data["MODE_ID"] == 3)[0]
ccres = hdu['LONG_TERM_CTI'].data["T_COEFF"][ix[1]]
#
ax[1].plot(xtime,ccres*1.0e4,color='green',label=f"New CalClosed")
ax[1].grid(True)
ax[1].legend()

plt.savefig(f'{wdir}/pn_sw_ltci_agn_ccf49.png',dpi=100)
plt.show()
plt.close()
#%%
# Save the new table to a file 
#
with open(f'{wdir}/ccf49_ltcti_table_mode3_agn.txt',"w") as fout:
    print ("{:3} : {:3} :   {:6.4f} : ".format(3,4,6.3990),end=" ",file=fout)
    for j in result:
        print (" {:3.6e}".format(j),end=" ",file=fout)
    print (":   1",file=fout)

#%%

hdus = fits.open(f'{wdir}/ccfdev/EPN_CTI_0049.CCF_test1')
#
#hdus[0].header['ISSUE'] = 49
ix = np.where(hdus['LONG_TERM_CTI'].data["MODE_ID"] == 3)[0]
#tab = hdus['LONG_TERM_CTI']
#nrows = tab.data.shape[0]
#xhdu = fits.BinTableHDU.from_columns(tab.columns, nrows=nrows+1)
#xhdu.data['MODE_ID'][nrows] = 3
#xhdu.data['CCD_ID'][nrows] = 4
#xhdu.data['ENERGY'][nrows] = 6.399
#xhdu.data['T_COEFF'][nrows] = result
#xhdu.data['SHIFT'][nrows] = 1
#idx_sorted = np.lexsort((xhdu.data['ENERGY'],xhdu.data['CCD_ID'],xhdu.data['MODE_ID']))
#xhdu.data = xhdu.data[idx_sorted]
hdus['LONG_TERM_CTI'].data[ix[2]]["T_COEFF"] = result
hdus['LONG_TERM_CTI'].data[ix[1]]["T_COEFF"][0] = 4.22e-4

#tab["T_COEFF"][ix[1]] = result
hdus.writeto(f"{wdir}/ccfdev/EPN_CTI_0049.CCF_test3",overwrite=True)




