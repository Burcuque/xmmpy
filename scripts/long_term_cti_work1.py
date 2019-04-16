#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 22 12:30:38 2018

Understanding the PN long-term CTI calibration file

@author: ivaltchanov
"""

#import os
import numpy as np
#import pandas as pd

import matplotlib.pylab as plt

#from astropy.io import fits
from astropy.table import Table

from scipy.interpolate import interp1d
from scipy import optimize

#%%
# read the CAL file with the long-term CTI correction
#
wdir = '/xdata/xcaldata/XMM/IVAN'
#
calNew = 'EPN_CTI_0048.CCF'
calOld = 'EPN_CTI_0047.CCF'
calfile = f'/xdata/ccf/pub/{calNew}'
calfile_old = f'/xdata/ccf/pub/{calOld}'
ltc = Table.read(calfile,hdu='LONG_TERM_CTI')
ltc.info()
ltc_old = Table.read(calfile_old,hdu='LONG_TERM_CTI')
# read the CAL file with the long-term CTI correction, the extension with the times
#
xt = Table.read(calfile,hdu='LTC_TIMES')
xt.info()
xtime = xt[0]['TIME'] # year since 01-01-2000T00:00:00
#
oldcal = False
try:
    xt2 = Table.read(calfile_old,hdu='LTC_TIMES')
    xtime2 = xt2[0]['TIME'] # year since 01-01-2000T00:00:00
except:
    oldcal = True    
#
#%%
#
# PN_SW is with MODE_ID=3
#
ix = np.where(ltc['MODE_ID'] == 3)[0]
tx = ltc[ix]
ix2 = np.where(ltc_old['MODE_ID'] == 3)[0]
tx2 = ltc_old[ix]
#
#%%
# let's check
#
#idx_line = 0 # Al Kaplpha at 1.4 keV
idx_line = 1 # Mn Kalpha at 5.8988 keV
ltci = tx[idx_line]['T_COEFF']
#ltci2 = tx2[idx_line]['T_COEFF']
#
# now the magic to convert to ratio
xcorr1 = np.power((1.0 - ltci)/(1.0 - ltci[0]),137.0)
xcorr0 = np.power((1.0 - ltci)/(1.0 - ltci[0]),190.0)
xcorr2 = np.power((1.0 - ltci)/(1.0 - ltci[0]),200.0)
# no interpolate, to be used later
#
f1 = interp1d(xtime, xcorr1)
f2 = interp1d(xtime, xcorr2)
f0 = interp1d(xtime, xcorr0)
#
fig = plt.figure()
ax = fig.subplots()
ax.plot(xtime,xcorr1,'o-',label=f'{calNew}, RAWY=137')
ax.plot(xtime,xcorr2,'o-',label=f'{calNew}, RAWY=200')
#ax.plot(xtime,f(xtime),'o-',label='interpol')
#ax.plot(xtime2,ltci2,'+-',label=calOld)
ax.legend()
plt.show()
#
#
#%%
#
# 
t = Table.read(f"{wdir}/PN_calclosed/pn_sw_calclosed_fit_results_nocti_sorted.csv",
               comment= "\s*#")
#t.sort("time")
#t.write("/xdata/xcaldata/XMM/IVAN/PN_calclosed/pn_sw_calclosed_fit_results_nocti_sorted.csv")
nt = len(t)
#
rev = t['rev'].data
rtime = t['time'].data
#
line2_lab = (0.162*5.888 + 0.6*5.899)/0.762 #eV probabilities from Wikipedia on Iron-55
line2 = t['mn1_cen'].data/1000.0
line2Err = t['mn1_cen_err'].data/1000.0
ratio = line2/line2_lab
weight = 1.0/line2Err**2
#%%
# now the interpolated values from the CCF
#
runyear = np.linspace(0.0,25.0,num=26)
#
ccf_int1 = f1(runyear)
ccf_int2 = f2(runyear)
ccf_int0 = f0(runyear)
#
fig = plt.figure(figsize=(12,8))
msize=10
ax = fig.subplots()
#print (np.log(ratio)/np.log(ccf_interpol))
#
ax.errorbar(rtime,line2/line2_lab,yerr=([line2Err/line2_lab,line2Err/line2_lab]),fmt='^',\
           markersize=msize,color='black', label='PN SW calclosed')
ax.plot(runyear,ccf_int1, '--', label=f'{calNew}, RAWY=137')
ax.plot(runyear,ccf_int0, color='black', label=f'{calNew}, RAWY=190')
ax.plot(runyear,ccf_int2, '--', label=f'{calNew}, RAWY=200')
#
ax.axhline(1.0,color='k',ls='dashed')
ax.set_ylabel("Mn K$_\\alpha$, E$_{obs}$/E$_{lab}$")
ax.set_xlabel(r"Years since 2000-01-01T00:00:00")
ax.set_title("PN SmallWindow Energy Scale Analysis\n Long-term CTI off")
#ax.set_title("{} Analysis".format(target.capitalize()))
ax.grid(True)
#
plt.legend()
#plt.savefig('pn_sw_ltci_check.png',dpi=100)
plt.show()
#plt.close()
#%%
# Now the fitting part
# We aregoing to fit the following function:
# ratio = [a + bt + ct^2 + dt^3 + et^4]^190.0
# then we have to encode [1 - (a + b't + c't^2 + d't^3 + e't^4)]/(1-a),
# where a is the CTI at t=0 ==> ltci[0]
# and store the values of this function:
# (a + b't + c't^2 + d't^3 + e't^4)

def f1(x, b, c, d):
    # polynomial function, will work ni LOG space to avoid overflows
    a0 = 0.00043236
    n0 = 1-a0
    return np.power((1.0 - (a0 + b*x + c*x**2 + d*x**3))/n0,190.0)

def f2(x, b, c, d):
    # polynomial function, will work ni LOG space to avoid overflows
    a0 = 0.00043236
    n0 = 1-a0
    return np.power((1.0 - (a0 + b*x + c*x**2 + d*x**3))/n0,190.0)

def residual1(p, x, y):
    return y - f1(x, *p)

def residual2(p, x, y):
    return y - f2(x, *p)

#
# now do this on pieces
#
i1 = np.where(rtime <= 3)[0]
i2 = np.where(rtime > 3)[0]

#p0 = [1.0e-3]
p0 = [1.0e-3, 1.0e-4, 1.0e-5]
popt1, pcov1 = optimize.leastsq(residual1, p0, args=(rtime[i1], ratio[i1]))

yn1 = f1(runyear, *popt1)

p0 = [1.0e-3, 1.0e-4, 1.0e-5]
popt2, pcov2 = optimize.leastsq(residual2, p0, args=(rtime[i2], ratio[i2]))

yn2 = f2(runyear, *popt2)

fig = plt.figure(figsize=(12,8))
ax = fig.subplots()
ax.plot(rtime, ratio, 'or')
ax.plot(runyear, yn1)
ax.plot(runyear, yn2)
ax.set_ylim([0.96,1.0])
ax.grid(True)

plt.show()

#%%
# try with piecewise function
#
def piecewise_linear(x, x0, y0, k1, k2):
    x0 = 1
    return np.piecewise(x, [x < x0], [lambda x:k1*x + y0-k1*x0, lambda x:k2*x + y0-k2*x0])

#def piecewise_linear(x, x0, y0, k1, k2):
#    x0 = 1.5
#    return np.piecewise(x, [x < x0], [lambda x:1.0+k1*x, lambda x:y0-k2*x0+k2*x])

#def piecewise_linear(x, x0, y0, k1, k2):
#    x0 = 3.0
#    return np.piecewise(x, [x < x0], [lambda x:1 + k2*x, lambda x:k2*x0-1 + k2*x])

p , e = optimize.curve_fit(piecewise_linear, rtime, ratio)
#
yn = piecewise_linear(runyear, *p)

fig = plt.figure(figsize=(12,8))
ax = fig.subplots()
ax.plot(rtime, ratio, 'or')
ax.plot(runyear, yn)
ax.set_ylim([0.96,1.0])
ax.grid(True)

plt.show()
#%%
#
# another idea, 
# 
#
def func1(x, b, c, d, e, f):
    a0 = 0.00043236
    return a0 + b*x + c*x**2 + d*x**3 + e*x**4 + f*x**5
#
def func2(x, b, c, d):
    a0 = 0.00043236
    return a0 + b*x + c*x**2 + d*x**3
#
a0 = 0.00043236
yy = 1.0 - (1.0-a0)*np.power(ratio,1.0/190.0)
#p0 = [1.0e-3]
ix1 = np.where(rtime > 3)[0]
p1 = [1.0e-3, 1.0e-4, 1.0e-5]
popt, pcov = optimize.curve_fit(func2, rtime[ix1], yy[ix1])

yn1 = np.power((1.0 - func2(runyear,*popt))/(1.0-a0),190.0)
yn200 = np.power((1.0 - func2(runyear,*popt))/(1.0-a0),200.0)
yn137 = np.power((1.0 - func2(runyear,*popt))/(1.0-a0),137.0)

fig = plt.figure(figsize=(12,8))
ax = fig.subplots()
ax.plot(rtime, ratio, 'or',label='PN SW calclosed')
ax.plot(runyear, yn137,'--',label='RAWY=137')
ax.plot(runyear, yn1,color='black',label='RAWY=190')
ax.plot(runyear, yn200,'--',label='RAWY=200')
ax.set_ylim([0.96,1.0])
ax.grid(True)
ax.set_ylabel("Mn K$_\\alpha$, E$_{obs}$/E$_{lab}$")
ax.set_xlabel(r"Years since 2000-01-01T00:00:00")
ax.set_title("PN SmallWindow Energy Scale Analysis\n Long-term CTI off")
plt.legend()
plt.savefig('pn_sw_ltci_newfit.png',dpi=100)
plt.show()
plt.close()
#
#%%
#
# now try with Splines
#
from scipy.interpolate import UnivariateSpline

i25 = np.where(rtime >= 2.5)[0]
weight[i25] = np.max(weight)
# normalize the weights
weight = weight/np.sum(weight)
#
a0 = 0.00043236
iend = len(rtime) - 1
#
rtime0 = np.insert(rtime,iend+1,[18.0,20.0,24.0])
ratio0 = np.insert(ratio,iend+1,[ratio[iend],ratio[iend],ratio[iend]])
wxx = np.insert(weight,iend+1,[0.5,0.5,0.5])
# normalize the weights
wxx = wxx/np.sum(wxx)
#
yy = 1.0 - (1.0-a0)*np.power(ratio,1.0/190.0)
yy0 = 1.0 - (1.0-a0)*np.power(ratio0,1.0/190.0)
#rtime0 = np.insert(rtime,0,0.0)
#yy0 = np.insert(yy,0,a0)


s = UnivariateSpline(rtime0, yy0, w=wxx,k=4,s=1)
#s = UnivariateSpline(rtime, yy, w=weight,s=1)
result = s(runyear)
#result[0] = a0

#yn1 = np.power((1.0 - s(runyear))/(1.0-a0),190.0)
#yn200 = np.power((1.0 - s(runyear))/(1.0-a0),200.0)
#yn137 = np.power((1.0 - s(runyear))/(1.0-a0),137.0)
#
yn1 = np.power((1.0 - result)/(1.0-a0),190.0)
yn200 = np.power((1.0 - result)/(1.0-a0),200.0)
yn137 = np.power((1.0 - result)/(1.0-a0),137.0)

fig = plt.figure(figsize=(12,8))
ax = fig.subplots()
ax.plot(rtime0, ratio0, 'sg',label='PN SW calclosed')
ax.plot(rtime, ratio, 'or',label='PN SW calclosed')
ax.plot(runyear, yn137,'--',label='RAWY=137')
ax.plot(runyear, yn1,color='black',label='RAWY=190')
ax.plot(runyear, yn200,'x--',label='RAWY=200')
#ax.set_xlim([0.6,1.0])
ax.set_ylim([0.96,1.0])
ax.grid(True)
ax.set_ylabel("Mn K$_\\alpha$, E$_{obs}$/E$_{lab}$")
ax.set_xlabel(r"Years since 2000-01-01T00:00:00")
ax.set_title("PN SmallWindow Energy Scale Analysis\n Long-term CTI off")
plt.legend()
plt.savefig(f'{wdir}/pn_sw_ltci_spline.png',dpi=100)
plt.show()
plt.close()
