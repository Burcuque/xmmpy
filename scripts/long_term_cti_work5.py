#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Derive the LTCTI only based on the Mn Ka and Fe Ka

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
targets = ["ngc2992","ngc3227","ngc3783","ngc4151","ngc4593",'ngc5548']
#
# redshifts
#
redshift = {'ngc4151': 0.003262, 'ngc3227': 0.00386, 'mrk1048': 0.0427, 'ngc3783': 0.009755,\
            'ngc4593': 0.008344, 'ngc5506': 0.00589, 'mcg-5-23-16': 0.008226, 'ngc3516': 0.008816,\
            'ngc5548': 0.01627, 'ngc2992':  0.007296, 'ngc1566': 0.005036, 'iras09149': 0.057150,\
            "iras05078": 0.017879, 'ngc7213': 0.005869}
feK = 6.399 # the weitghed mean of Kalpha_2 at 6.3908 (intensity 50) and Kalpha_1 at 6.40308 (intensity 100)
#
wdir = "/xdata/xcaldata/XMM/IVAN/PN_SW/sources"
os.chdir(wdir)
#
# READ the CCF, will only use the new LTCTI at 5.9 keV from CalClosed data
#
#hdu = fits.open("/home/ivaltchanov/IVAN/PN_SW/ccfdev/EPN_CTI_0049.CCF_test14")
# test14 is the previous correction with good results
hdu = fits.open("ccfdev/EPN_CTI_0049.CCF_test27")
# test4 is the previous correction with no extrapolation for the reference energy
tt = hdu['LTC_TIMES'].data['TIME'][0]
xx = hdu['LONG_TERM_CTI'].data
sel = np.where(xx['MODE_ID'] == 3)[0]
xsel = xx[sel]
e1 = xsel['ENERGY'][0] # AlKa at 1.486 keV
e2 = xsel['ENERGY'][1] # MnKa at 5.8988
tcoeff1 = xsel['T_COEFF'][0]
tcoeff2 = xsel['T_COEFF'][1]
# now calculate the correction
corr1 = np.power((1.0 - tcoeff1)/(1.0 - tcoeff1[0]),170.0)
corr2x = np.power((1.0 - tcoeff2)/(1.0 - tcoeff2[0]),170.0)
#%%
msize=10
fig = plt.figure(figsize=(12,8))
ax = fig.add_subplot(111)
#
# plot the CalClosed data and curve
#,comment="\s*#"
ccdir = '/xdata/xcaldata/XMM/IVAN/PN_SW/CalClosed'
#tcc = Table.read(f"{ccdir}/pn_sw_calclosed_fit_results_nocti_sorted.csv",comment="\s*#")
tcc = Table.read(f"{ccdir}/cc_output_xspec_no_cti_sorted.csv",comment="\s*#")
tcc.sort('time')
#tcc = tcc[np.where(tcc['time'] >= 1.2)[0]]
rtime = tcc['time']
#
# Al Ka
line1_lab = 1.486 #keV probabilities from Wikipedia on Iron-55
line1 = tcc['al'].data
line1Err = tcc['al_err'].data
#ax.plot(tt,corr1,'c-',label='Correction for 1.486 keV, RAWY=170')
#ax.errorbar(rtime,line1/line1_lab,line1Err/line1_lab,fmt='h',ms=10,label=r'AlK$\alpha$')

# Mn Ka
#line2_lab = (0.162*5.888 + 0.6*5.899)/0.762 #keV probabilities from Wikipedia on Iron-55
line2_lab = 5.8988
line2 = tcc['mn1'].data
line2Err = tcc['mn1_err'].data
ax.errorbar(rtime,line2/line2_lab,line2Err/line2_lab,fmt='^',color='black',ms=10,label=r'MnK$\alpha$')

ratio = line2/line2_lab
ratio_err = line2Err/line2_lab

xsp = rtime
ysp = ratio
ysp_err = ratio_err
#
#isort = np.argsort(xsp)
#
# force the curve to pass through (0.0,1.0), 
# need to put a weight of 1 and 0.5 to the rest
#xsp = np.insert(xsp,0,0.0)
#ysp = np.insert(ysp,0,1.0)
#ysp_err = np.insert(ysp_err,0,1.0e-6)
weights = 1.0/ysp_err

#s = UnivariateSpline(xsp, ysp, w=weights, k=3,s=2)
s = UnivariateSpline(xsp, ysp, w=weights, k=1,s=1)

result = s(tt)
result[0] = 1.0
#result[1] = result[1]/1.001
#result[2] = result[2]/1.001
#result[3] = result[3]/1.0005
#result[16] = result[16]/1.001
#result[16:] = result[16]
#out = result.copy()
#ixq = np.where(tt <= 10.0)[0]
#out[ixq] = corr2x[ixq]
#
#ax.errorbar(xsp,ysp,ysp_err,fmt='*',ms=10,label=r'MnK$\alpha$')
ax.plot(tt,result,'r-', label=r'Correction for Mn K$\alpha$, RAWY=170')
#ax.plot(tt,out,'k-', label=r'Adopted for Mn K$\beta$, RAWY=170')

# now derive g(t)
a0 = 0.00043236
gt1 = 1.0 - (1.0-a0)*np.power(result,1.0/170.0)
#gt1 = 1.0 - (1.0-a0)*np.power(corr2x,1.0/170.0)
gt1[0] = a0

corr2z = np.power((1.0 - gt1)/(1.0 - gt1[0]),190.0)
ax.plot(tt,corr2z,'g--',label='Correction for 5.9 keV, RAWY=190')

#iqw = np.where(tt < 13.0)[0]
#corr2y = corr2x.copy()
#corr2y[iqw] = corr2y[iqw]/1.001
#corr2y = corr2y*1.001
#corr2y[0] = 1.0
#ax.plot(tt,corr2x,'g--',label='Correction for 5.9 keV, RAWY=170')
#ax.plot(tt,corr2y,'b--',label='Correction for 5.9 keV, RAWY=170')


# Mn Kb
line3_lab = 6.49 #keV probabilities from Wikipedia on Iron-55
line3 = tcc['mn2'].data
line3Err = tcc['mn2_err'].data
#ax.errorbar(rtime,ratio,ratio_err,fmt='*',ms=10,label=r'MnK$\beta$')
#
#
tab_file = 'selected_targets_nocti.csv'
tq = Table.read(tab_file,comment="\s*#" )
nq = len(tq)
lineX = []
for i in np.arange(nq):
    target = tq['target'].data[i]
    lineX.append(feK/(1.0 + redshift[target]))
#    
rev = tq['rev'].data
rtime = tq['delta_time'].data
line = tq['lineE'].data
lineErrUp = tq['lineE_low'].data
lineErrLo = tq['lineE_up'].data
lineErr = (lineErrLo + lineErrUp)/2.0
rchi2 = tq["chi2r"].data
ratio = line/lineX
ratio_err = lineErr/lineX

# AGN derived curve
#ax.plot(tt,corr3x,'r--', label=r'Correction from AGN at RAWY=170')

ax.errorbar(rtime,ratio,ratio_err,fmt='o',color='red',\
            label=fr'AGN FeK$\alpha$ at 6.4 keV')

#
ax.set_ylim((0.97,1.0,))
#ax.set_xlim((0.0,3.0,))
ax.set_ylabel(r'$E_{obs}/E_{lab}$')
ax.set_xlabel("Years since 2000-01-01T00:00:00")
ax.legend(numpoints=1)
#ax.set_zlabel('E_obs/E_lab')
plt.grid(True)
plt.title('PN SW CalClosed without Long-term CTI correction')
#plt.savefig(f'{wdir}/cc_curves_t18.png',dpi=100)
plt.show()
plt.close()
#%%
# now replace T_COEFF in the calibration file
#
ccf49_file = '/home/ivaltchanov/IVAN/PN_SW/ccfdev/EPN_CTI_0049.CCF_test24'
#
# will use the AGN derived curve, extrapolated to  at 6.4 keV
#
hdu49 = fits.open(f"{ccf49_file}")
#hdu49.close()
ltc49 = hdu49['LONG_TERM_CTI']
ix49 = np.where(ltc49.data['MODE_ID'] == 3)[0]
xtcoeff = hdu49['LONG_TERM_CTI'].data['T_COEFF']
##
hdu49['LONG_TERM_CTI'].data['ENERGY'][ix49[1]] = f"{line2_lab:.4f}"
hdu49['LONG_TERM_CTI'].data['T_COEFF'][ix49[1]] = gt1
hdu49.writeto("/xdata/xcaldata/XMM/IVAN/PN_SW/ccfdev/EPN_CTI_0049.CCF_test27",overwrite=True)

