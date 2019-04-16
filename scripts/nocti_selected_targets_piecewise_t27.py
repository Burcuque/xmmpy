#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Implement energy extrapolation

The idea is the following:
    AGN observed at tx, FeK is at linex (due to redshift), with line/lab = ratiox
    then at tx we have two points, one at Energy 5.9 keV (from CC MnKa) and 
    another one at Energy linex with two corrections: cx and ratiox.
    The we extrapolate liearly in logE space the correction at energy 6.4 keV
    This will be the reference point in the CCF file

This is an implementation of an idea, from Norbert, to incorporate the energy dependence in the 
derivation of the correction curve.

Note: uses only selected targets and OBSIDs with no-cti results

Crafted: 10-01-2019

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
#target = "ngc4593"
#target = "ngc4151"
#target = "mrk1048"
#target = "ngc3227"
#target = "ngc3783"
#target = "ngc5506"
#target = "mcg-5-23-16"
#target = "ngc5548"
#target = "ngc3516"
#target = "WR140"
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
control = ["ngc1566","iras09149","ngc5506","iras05078","ngc7213"]
#
# subset of OBSIDs to include in the analysis
#

#if (target == 'WR140'):
#    lineX = 6.697 # see Sugawara et al. (2015)
#
wdir = "/xdata/xcaldata/XMM/IVAN/PN_SW/sources"
#wdir = "/home/ivaltchanov/XMM/CAL/PN_SW/{}".format(target)
#wdir = "/home/ivaltchanov/XMM/CAL/PN_SW/WR140"
os.chdir(wdir)
#
# READ the CCF, will only use the new LTCTI at 5.9 keV from CalClosed data
#
hdu = fits.open("/home/ivaltchanov/IVAN/PN_SW/ccfdev/EPN_CTI_0049.CCF_test27")
# test4 is the previous correction with no extrapolation for the reference energy
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
#%%
# Extrapolate and plot 
#
# log10 of the observed line energy, needed for the linear extrapolation 
xxe = np.log10(line)
nd = len(line)
ex = 6.4 # keV, the reference energy to use
logEx = np.log10(ex)
loge1 = np.log10(e1) # the Al Ka line at 1.48 keV, not used
loge2 = np.log10(e2) # the Mn Ka line at 5.9 keV, will be used
out = []
#
fig = plt.figure(figsize=(12,8))
ax = fig.add_subplot(111)
for i in np.arange(nd):
    # now find the line that passes through the two points
    # y - y2 = (y2-y1)/(x2-x1) * (x - x2)
    ccx = np.interp(rtime[i],tt,corr2)
    dx = xxe[i] - loge2  # x2-x1
    dy = ratio[i] - ccx # y2-y1
    yout = dy*(logEx - xxe[i])/dx + ratio[i]
    #ax.semilogx([e2,line[i],ex],[ccx,ratio[i],yout],'o-')
    ax.plot([loge2,xxe[i],logEx],[ccx,ratio[i],yout],'o-')
    out.append(yout)
    # now find the line that passes through the two points
    #xe = [np.log10(e1),np.log10(e2),xxe[i]]
    #ye = [np.interp(yy[i],tt,corr1),np.interp(yy[i],tt,corr2),zz[i]]
    #cc = np.polyfit(xe,ye,1)
    #p = np.poly1d(cc)
    #out.append(p(log64))
#ax.set_xlim([5.8,6.5])
ax.set_xlabel('Log10(Energy keV)')
ax.set_ylabel(r'$E_{obs}/E_{lab}$')
ax.grid(True)
plt.title('PN SW extrapolation')
plt.savefig(f'{wdir}/linear_extrapolation_t28.png',dpi=100)
plt.show()
plt.close()
#%%    
#
# The correction based on th extrapolated Eobs/Elab at Eref=6.4 keV
#
# just for comparison with test7 curve
hdu = fits.open("/home/ivaltchanov/IVAN/PN_SW/ccfdev/EPN_CTI_0049.CCF_test27")
# test4 is the previous correction with no extrapolation for the reference energy
tt7 = hdu['LTC_TIMES'].data['TIME'][0]
xx7 = hdu['LONG_TERM_CTI'].data
sel7 = np.where(xx7['MODE_ID'] == 3)[0]
xsel7 = xx7[sel7]
tcoeff37 = xsel['T_COEFF'][2]
# now calculate the correction
corr37 = np.power((1.0 - tcoeff37)/(1.0 - tcoeff37[0]),190.0)

#
out = np.asarray(out) 
msize=10
fig = plt.figure(figsize=(12,8))
ax = fig.add_subplot(111)
#
# plot the CalClosed data and curve
#
#tcc = Table.read("/xdata/xcaldata/XMM/IVAN/PN_SW/CalClosed/pn_sw_calclosed_fit_results_nocti.csv")
#tcc = Table.read("/xdata/xcaldata/XMM/IVAN/PN_SW/CalClosed/pn_sw_calclosed_fit_results_nocti.csv")
ccdir = '/xdata/xcaldata/XMM/IVAN/PN_SW/CalClosed'
tcc = Table.read(f"{ccdir}/cc_output_xspec_no_cti_sorted.csv",comment="\s*#")
#rtime = tcc['time']
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


ax.errorbar(rtime,out,ratio_err,fmt='o',color='red',\
            label=fr'AGN FeK$\alpha$ at {ex:.3} keV')

result  = corr37.copy()

#it10 = np.where(tt >= 12.0)[0]
#result[it10] = result[it10]/1.0025
#
#it10 = np.where((tt < 9.8)*(tt>=4.0))[0]
#result[it10] = corr2x[it10]/1.0005
#
#it10 = np.where(tt >= 9.8)[0]
#result[it10] = s(tt[it10])

result[13:] = result[13]

ax.plot(tt,result,'r--', label=f'Correction for {ex} keV, RAWY=190')

#ax.plot(tt7,corr37,'g--', label='Previous')
#
# now derive g(t)
#
a0 = 0.00043236
gt = 1.0 - (1.0-a0)*np.power(result,1.0/190.0)
gt[0] = a0
#print (gt)

ax.set_ylabel(r'$E_{obs}/E_{lab}$')
ax.set_xlabel("Years since 2000-01-01T00:00:00")
ax.legend(numpoints=1)
#ax.set_zlabel('E_obs/E_lab')
plt.grid(True)
plt.title('PN SW without Long-term CTI correction')
plt.savefig(f'{wdir}/correction_with_extrapolation_t28.png',dpi=100)
plt.show()
plt.close()
#%%
#
# now replace T_COEFF in the calibration file
#
ccf49_file = '/home/ivaltchanov/IVAN/PN_SW/ccfdev/EPN_CTI_0049.CCF_test27'
hdu49 = fits.open(f"{ccf49_file}")
#hdu49.close()
ltc49 = hdu49['LONG_TERM_CTI']
ix49 = np.where(ltc49.data['MODE_ID'] == 3)[0]
xtcoeff = hdu49['LONG_TERM_CTI'].data['T_COEFF']
##
hdu49['LONG_TERM_CTI'].data['ENERGY'][ix49[2]] = ex
hdu49['LONG_TERM_CTI'].data['T_COEFF'][ix49[2]] = gt
hdu49.writeto("/xdata/xcaldata/XMM/IVAN/PN_SW/ccfdev/EPN_CTI_0049.CCF_test28",overwrite=True)

