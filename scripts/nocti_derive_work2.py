#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Implement energy extrapolation

The CalClosed data give a point in the (time,energy,LTCTI) 3-d plane

For AGNs, each one has the FeKa line at differnet energy, depending on redshift

The idea is the following:
    AGN observed at tx, FeK is at linex (due to redshift), with line/lab = ratiox
    then at tx we have two points, one at Energy 5.9 keV (from CC MnKa) and 
    another one at Energy linex with two corrections: cx and ratiox.
    The we liearly (in logE space) extrapolate the correction at energy 6.4 keV
    This will be the reference point in the CCF file

This is an implementation of an idea, from Norbert, to incorporate the energy dependence in the 
derivation of the correction curve.

x = E_obs
y = time
z = E_obs/E_lab/(1+z)

Note #1: uses all AGNs, no selecttion

Note #2: this script is superseded by nocti_selected_targets.py


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
sns.set(style="darkgrid")
#import matplotlib as mpl
#mpl.style.use('classic')

#from matplotlib.patches import Patch
#from matplotlib.lines import Line2D
#from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 unused import
#%matplotlib

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
wdir = '/home/ivaltchanov/IVAN/PN_SW'
hdu = fits.open(f"{wdir}/ccfdev/EPN_CTI_0049.CCF")
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
corr2 = np.power((1.0 - tcoeff2)/(1.0 - tcoeff2[0]),190.0)
corr3 = np.power((1.0 - tcoeff3)/(1.0 - tcoeff3[0]),190.0)
#%%
xx = []
yy = []
zz = []
zz_err = []
#
for target in targets:
    print (f"Adding {target}")
    tab_file = f'{target}/output_xspec_nocti.csv'
    if (not os.path.isfile(tab_file)):
        print (f"No nocti results for {target}")
        continue
    tq = Table.read(tab_file,comment="\s*#" )
    iq = np.where(tq['chi2r'].data <= 1.5)[0]
    if (len(iq) < 1):
        print (f"All fit results for {target} have rChi2 > 1.3, discarding")
        continue
    t = tq[iq]
    nt = len(t)
    #
    lineX =  feK/(1.0 + redshift[target]) # redshifted line position
    rev = t['rev'].data
    rtime = t['delta_time'].data
    ff = t['full'].data
    inst = t['inst'].data
    line = t['lineE'].data
    lineErrUp = t['lineE_low'].data
    lineErrLo = t['lineE_up'].data
    lineErr = (lineErrLo + lineErrUp)/2.0
    rchi2 = t["chi2r"].data
    ratio = line/lineX
    ratio_err = lineErr/lineX
    #iqqx = np.where((ratio < 0.991)*(rtime < 2.0))[0]
    #if (len(iqqx) > 0):
    #    print ("Skipping {}, revoution {}".format(target,rev[iqqx]))
    print ("{}: min line energy {}, max line energy {} ".format(target,np.min(line),np.max(line)))
    xx = np.concatenate((xx,line),axis=None)
    yy = np.concatenate((yy,rtime),axis=None)
    zz = np.concatenate((zz,ratio),axis=None)
    zz_err = np.concatenate((zz_err,ratio_err),axis=None)
    #

#%%
#
e_mean = np.mean(line)
e_median = np.median(line)

xxe = np.log10(xx)
nd = len(xx)
ex = 6.5 # keV, the reference energy to use
#ex = e3
#ex = e_mean
logEx = np.log10(ex)
loge1 = np.log10(e1)
loge2 = np.log10(e2)
out = []
fig = plt.figure(figsize=(12,8))
ax = fig.add_subplot(111)
for i in np.arange(nd):
    # now find the line that passes through the two points
    # y - y2 = (y2-y1)/(x2-x1) * (x - x2)
    ccx = np.interp(yy[i],tt,corr2)
    dx = xxe[i] - loge2  # x2-x1
    dy = zz[i] - ccx # y2-y1
    yout = dy*(logEx - xxe[i])/dx + zz[i]
    ax.plot([loge2,xxe[i],logEx],[ccx,zz[i],yout],'o-')
    out.append(yout)
    # now find the line that passes through the two points
    #xe = [np.log10(e1),np.log10(e2),xxe[i]]
    #ye = [np.interp(yy[i],tt,corr1),np.interp(yy[i],tt,corr2),zz[i]]
    #cc = np.polyfit(xe,ye,1)
    #p = np.poly1d(cc)
    #out.append(p(log64))
ax.set_xlabel('Log10(Energy) keV')
ax.set_ylabel('Eobs/Elab')
plt.savefig(f'{wdir}/linear_extrapolation_t9.png',dpi=100)
plt.show()
plt.close()
#%%    
#
out = np.asarray(out) 
msize=10
fig = plt.figure(figsize=(12,8))
ax = fig.add_subplot(111)
ax.plot(yy,zz,'bs',label='Original')
ax.plot(tt,corr3,'b-',label='Original correction curve')
ax.plot(yy,out,'ro',label=f'Extrapolated to {ex:.3} keV')
#
# now the univariate spline
#
isort = np.argsort(yy)
#
# force the curve to pass through (0.0,1.0), 
# need to put a weight of 1 and 0.5 to the rest
xsp = np.insert(yy[isort],0,0.0)
ysp = np.insert(out[isort],0,1.0)
weights = np.zeros_like(xsp) + 0.5
weights[0] = 1.0
#isp = np.where(ysp <= 0.97)[0]
#weights[isp] = 0.3
s = UnivariateSpline(xsp, ysp, w=weights, k=5,s=1)

result = s(tt)
result[17:] = result[16]

ax.plot(tt,result,'r--', label='Updated correction curve')
#
# now derive g(t)
#
a0 = 0.00043236
gt = 1.0 - (1.0-a0)*np.power(result,1.0/190.0)
gt[0] = a0
print (gt)

#xcorr3 = np.power((1.0 - gt)/(1.0 - gt[0]),190.0)
#ax.plot(tt,xcorr3,'g-', label='Updated correction curve')

#qcorr3 = corr3.copy()
#ix = np.where(tt >= 10)[0]
#qcorr3[ix] = qcorr3[ix]*1.0015
#ix = np.where((tt >= 7)*(tt < 10))[0]
#qcorr3[ix] = qcorr3[ix]*1.00015
#ax.plot(tt,qcorr3,'r--')
ax.set_ylabel('Eobs/Elab')
ax.set_xlabel('Time')
ax.legend()
#ax.set_zlabel('E_obs/E_lab')
plt.grid(True)
plt.title('PN SW')
plt.savefig(f'{wdir}/correction_with_extrapolation_t9.png',dpi=100)
plt.show()
plt.close()
#%%
#
# now replace T_COEFF in the calibration file
#
ccf49_file = f'{wdir}/ccfdev/EPN_CTI_0049.CCF'
hdu49 = fits.open(f"{ccf49_file}")
hdu49.close()
ltc49 = hdu49['LONG_TERM_CTI']
ix49 = np.where(ltc49.data['MODE_ID'] == 3)[0]
xtcoeff = hdu49['LONG_TERM_CTI'].data['T_COEFF']
##
hdu49['LONG_TERM_CTI'].data['ENERGY'][ix49[2]] = ex
hdu49['LONG_TERM_CTI'].data['T_COEFF'][ix49[2]] = gt
#%%
#hdu49.writeto(f"{wdir}/ccfdev/EPN_CTI_0049.CCF_test7",overwrite=True)

