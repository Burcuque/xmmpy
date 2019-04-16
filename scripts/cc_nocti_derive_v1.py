#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Derive the LTCTI only based on the Mn Ka and Kb lines

Will use the previous Mn Ka curve and derive a new one based on Mn Kb

Use the XSPEC fit results fot CalClosed 

Created: 15-01-2019

@author: ivaltchanov
"""
import os

import numpy as np
#from numpy.polynomial import polynomial
#import scipy

from astropy.table import Table
from astropy.io import fits
from scipy.interpolate import UnivariateSpline

import matplotlib
import matplotlib.pylab as plt
import seaborn as sns
sns.set(style="white")

matplotlib.rcParams['axes.formatter.useoffset'] = False
plt.rc('text', usetex=False)
plt.rc('font', family='serif')
#%%
#
wdir = "/xdata/xcaldata/XMM/IVAN/PN_SW/CalClosed"
os.chdir(wdir)
#
# READ the CCF, will only use the new LTCTI at 5.9 keV from CalClosed data
#
hdu = fits.open("/home/ivaltchanov/IVAN/PN_SW/ccfdev/EPN_CTI_0049.CCF")
# test4 is the previous correction with no extrapolation for the reference energy
tt = hdu['LTC_TIMES'].data['TIME'][0]
xx = hdu['LONG_TERM_CTI'].data
sel = np.where(xx['MODE_ID'] == 3)[0]
xsel = xx[sel]
e1 = xsel['ENERGY'][0] # AlKa at 1.486 keV
e2 = xsel['ENERGY'][1] # MnKa at 5.8988
tcoeff1 = xsel['T_COEFF'][0]
tcoeff2 = xsel['T_COEFF'][1]
tcoeff3 = xsel['T_COEFF'][2]
# now calculate the correction
corr1 = np.power((1.0 - tcoeff1)/(1.0 - tcoeff1[0]),190.0)
corr2x = np.power((1.0 - tcoeff2)/(1.0 - tcoeff2[0]),170.0)
corr3x = np.power((1.0 - tcoeff3)/(1.0 - tcoeff3[0]),170.0) # AGN derived curve at RAWY=170
#%%
msize=10
fig = plt.figure(figsize=(12,8))
ax = fig.add_subplot(111)
#
# plot the CalClosed data and curve
#
ccdir = '/xdata/xcaldata/XMM/IVAN/PN_SW/CalClosed'
tcc = Table.read(f"{ccdir}/pn_sw_calclosed_fit_results_nocti_sorted.csv",comment="\s*#")
rtime = tcc['time']
#
# Mn Ka
line2_lab = (0.162*5.888 + 0.6*5.899)/0.762 #keV probabilities from Wikipedia on Iron-55
line2 = tcc['mn1_cen'].data/1000.0
line2Err = tcc['mn1_cen_err'].data/1000.0
ax.plot(tt,corr2x,'b--',label='Correction for 5.9 keV, RAWY=170')
ax.errorbar(rtime,line2/line2_lab,line2Err/line2_lab,fmt='^',ms=10,label=r'MnK$\alpha$')
# Mn Kb
line3_lab = 6.49 #keV probabilities from Wikipedia on Iron-55
line3 = tcc['mn2_cen'].data/1000.0
line3Err = tcc['mn2_cen_err'].data/1000.0
ratio = line3/line3_lab
ratio_err = line3Err/line3_lab
#ax.errorbar(rtime,ratio,ratio_err,fmt='*',ms=10,label=r'MnK$\beta$')

#iq = np.where((rtime <= 3.0)*(ratio >= 0.992) + (rtime >= 3.0))[0]
#xsp = rtime[iq]
#ysp = ratio[iq]
#ysp_err = ratio_err[iq]
xsp = rtime
ysp = ratio
ysp_err = ratio_err
#
isort = np.argsort(xsp)
#
# force the curve to pass through (0.0,1.0), 
# need to put a weight of 1 and 0.5 to the rest
xsp = np.insert(xsp[isort],0,0.0)
ysp = np.insert(ysp[isort],0,1.0)
ysp_err = np.insert(ysp_err[isort],0,1.0e-6)
weights = 1.0/ysp_err

s = UnivariateSpline(xsp, ysp, w=weights, k=1,s=1)

result = s(tt)
result[18:] = result[18]
out = result.copy()
ixq = np.where(tt <= 10.0)[0]
out[ixq] = corr2x[ixq]
#
ax.errorbar(xsp,ysp,ysp_err,fmt='*',ms=10,label=r'MnK$\beta$')
ax.plot(tt,result,'g-', label=r'Correction for Mn K$\beta$, RAWY=170')
ax.plot(tt,out,'k-', label=r'Adopted for Mn K$\beta$, RAWY=170')
#
# AGN derived curve
ax.plot(tt,corr3x,'r--', label=r'Correction from AGN at RAWY=170')
#
# now derive g(t)
#
a0 = 0.00043236
gt = 1.0 - (1.0-a0)*np.power(out,1.0/170.0)
gt[0] = a0
#print (gt)
#ax.set_xlim((0.0,2.0))
ax.set_ylabel(r'$E_{obs}/E_{lab}$')
ax.set_xlabel("Years since 2000-01-01T00:00:00")
ax.legend(numpoints=1)
#ax.set_zlabel('E_obs/E_lab')
plt.grid(True)
plt.title('PN SW CalClosed without Long-term CTI correction')
plt.savefig(f'{wdir}/cc_curves_t11.png',dpi=100)
plt.show()
#plt.close()
#%%
#
# now replace T_COEFF in the calibration file
#
ccf49_file = '/home/ivaltchanov/IVAN/PN_SW/ccfdev/EPN_CTI_0049.CCF'
hdu49 = fits.open(f"{ccf49_file}")
#hdu49.close()
ltc49 = hdu49['LONG_TERM_CTI']
ix49 = np.where(ltc49.data['MODE_ID'] == 3)[0]
xtcoeff = hdu49['LONG_TERM_CTI'].data['T_COEFF']
##
hdu49['LONG_TERM_CTI'].data['ENERGY'][ix49[2]] = 6.49
hdu49['LONG_TERM_CTI'].data['T_COEFF'][ix49[2]] = gt
hdu49.writeto("/xdata/xcaldata/XMM/IVAN/PN_SW/ccfdev/EPN_CTI_0049.CCF_test11",overwrite=True)

