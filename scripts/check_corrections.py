#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

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
##
#targets = ["ngc4593","ngc5506","ngc5548", "ngc4151","ngc3516","ngc2992","ngc3227",'ngc3783','ngc1566']
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
hdu = fits.open("/home/ivaltchanov/IVAN/PN_SW/ccfdev/EPN_CTI_0049.CCF_test27")
tt = hdu['LTC_TIMES'].data['TIME'][0]
xx = hdu['LONG_TERM_CTI'].data
sel = np.where(xx['MODE_ID'] == 3)[0]
xsel = xx[sel]
e1 = xsel['ENERGY'][0] # AlKa at 1.486 keV
e2 = xsel['ENERGY'][1] # MnKa at 5.8988
e3 = xsel['ENERGY'][2] # MnKa at 5.8988
tcoeff1 = xsel['T_COEFF'][0]
tcoeff2 = xsel['T_COEFF'][1]
tcoeff3 = xsel['T_COEFF'][2]
# now calculate the correction
corr1 = np.power((1.0 - tcoeff1)/(1.0 - tcoeff1[0]),190.0)
corr2 = np.power((1.0 - tcoeff2)/(1.0 - tcoeff2[0]),190.0) # needed for extrapolation
corr2x = np.power((1.0 - tcoeff2)/(1.0 - tcoeff2[0]),170.0)
corr3 = np.power((1.0 - tcoeff3)/(1.0 - tcoeff3[0]),190.0)
#%%
# This is the table with noCTI results on selected targets/OBSIDs
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
#
msize=10
fig = plt.figure(figsize=(12,8))
ax = fig.add_subplot(111)
#
ccdir = '/xdata/xcaldata/XMM/IVAN/PN_SW/CalClosed'
tcc = Table.read(f"{ccdir}/cc_output_xspec_no_cti_sorted.csv",comment="\s*#")
#
# Mn Ka
#line2_lab = (0.162*5.888 + 0.6*5.899)/0.762 #keV probabilities from Wikipedia on Iron-55
line2_lab = 5.8988
line2 = tcc['mn1'].data
line2Err = tcc['mn1_err'].data
ax.plot(tt,corr2x,'b--',label='Correction for 5.9 keV, RAWY=170')
ax.errorbar(tcc['time'],line2/line2_lab,line2Err/line2_lab,fmt='^',ms=10,label=r'MnK$\alpha$')
# Mn Kb
line3_lab = 6.49 #keV probabilities from Wikipedia on Iron-55
line3 = tcc['mn2'].data
line3Err = tcc['mn2_err'].data
#ax.plot(tt,corr2x,'b--',label='Correction for 5.9 keV, RAWY=170')
#ax.errorbar(tcc['time'],line3/line3_lab,line3Err/line3_lab,fmt='*',ms=10,label=r'MnK$\beta$')
#ax.plot(rtime,ratio,'bs',label='No CTI results')
#ax.plot(tt,corr37,'y-',label='Correction for 6.4 keV, RAWY=190')


ax.errorbar(rtime,ratio,ratio_err,fmt='o',color='red',\
            label=fr'AGN FeK$\alpha$ at 6.4 keV')


ax.plot(tt,corr3,'r--',label='Correction for 6.4 keV')
#
# try to fix it
#
corr3x = corr3.copy()
corr3x[8:] = corr3x[8:]*1.002

ax.plot(tt,corr3x,'g--',label='Correction for 6.4 keV')


ax.set_ylabel(r'$E_{obs}/E_{lab}$')
ax.set_xlabel("Years since 2000-01-01T00:00:00")
ax.legend(numpoints=1)
#ax.set_zlabel('E_obs/E_lab')
plt.grid(True)
plt.title('PN SW without Long-term CTI correction')
#plt.savefig(f'{wdir}/correction_with_extrapolation_t25.png',dpi=100)
plt.show()
plt.close()
#
#%%
#
