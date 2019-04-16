#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 20 11:29:20 2018

@author: ivaltchanov
"""
import os
import numpy as np 
from datetime import datetime
import matplotlib.pyplot as plt

from astropy.io import fits
from astropy.table import Table
from astropy.convolution import convolve, Box1DKernel

from lmfit.models import GaussianModel, PolynomialModel

def fit_cal_lines(fitsfile, minMax=(500,10000), plotIt=True, savePng=False, \
    pngName = "fitted.png", verbose=False, smooth=None):
    #
    # read a spectrum form CalClosed observation and with the Al and Mn lines
    # minMax is the tuple for the min and max channel to include in the fit
    # good values are (20,400)
    hdu = fits.open(fitsfile)
    x = hdu['SPECTRUM'].data["CHANNEL"]*5.0 # in eV
    peak1 = hdu['SPECTRUM'].data["COUNTS"]
    if ('EXPOSURE' in hdu['SPECTRUM'].header.keys()):
        expo1 = hdu['SPECTRUM'].header['EXPOSURE']
    elif ('EXPOSURE' in hdu[0].header.keys()):
        expo1 = hdu[0].header['EXPOSURE']
    elif ('DURATION' in hdu[0].header.keys()):
        expo1 = hdu[0].header['DURATION']
    else:
        print ("Warning: cannot find the exposure, stting it to 1.0")
        expo1 = 1.0
    #
    obsid = hdu[0].header["OBS_ID"]
    revol = hdu[0].header["REVOLUT"]
    window = hdu[0].header["SUBMODE"]
    start_time = hdu[0].header['DATE-OBS']
    #
    # Times will be relative to 2000-01-01T00:00:00
    #
    time0 = datetime.strptime("2000-01-01T00:00:00","%Y-%m-%dT%H:%M:%S")
    stime = datetime.strptime(start_time,"%Y-%m-%dT%H:%M:%S")
    delta_time = (stime-time0).total_seconds()/(365.0*24.0*3600.0) # in years
    
    if (smooth == None):
        y = peak1
    else:
        y = convolve(peak1, Box1DKernel(smooth))
    #
    m1 = x >= minMax[0]
    m2 = x <= minMax[1]
    x = x[m1*m2]
    y = y[m1*m2]
    ig1 = (x >= 750) * (x <= 1600)
    xq = x[ig1]
    yq = y[ig1]
    i1max = np.argmax(yq)
    y1max = yq[i1max]
    x1max = xq[i1max]
    ig2 = (x >= 5600) * (x <= 6100)
    xq = x[ig2]
    yq = y[ig2]
    i2max = np.argmax(yq)
    y2max = yq[i2max]
    x2max = xq[i2max]
    ig3 = (x >= 6200) * (x <= 6700)
    xq = x[ig3]
    yq = y[ig3]
    i3max = np.argmax(yq)
    y3max = yq[i3max]
    x3max = xq[i3max]
    #print (x1max,x2max,x3max)
    #
    poly_mod = PolynomialModel(2,prefix='poly_')
    pars = poly_mod.guess(y, x=x)
    #
    gauss1  = GaussianModel(prefix='g1_')
    pars.update( gauss1.make_params())
    pars['g1_center'].set(x1max,min=x1max-100.0,max=x1max+100.0)
    pars['g1_sigma'].set(80.0,min=40,max=150)
    pars['g1_amplitude'].set(y1max)
    #
    gauss2  = GaussianModel(prefix='g2_')
    pars.update(gauss2.make_params())
    pars['g2_center'].set(x2max,min=x2max-100.0,max=x2max+100.0)
    pars['g2_sigma'].set(100.0,min=60,max=150)
    pars['g2_amplitude'].set(y2max)
    #
    gauss3  = GaussianModel(prefix='g3_')
    pars.update(gauss3.make_params())
    pars['g3_center'].set(x3max,min=x3max-100.0,max=x3max+100.0)
    pars['g3_sigma'].set(100.0,min=60,max=160)
    pars['g3_amplitude'].set(y3max)
    #
    mod = gauss1 + gauss2 + gauss3 + poly_mod
    #init = mod.eval(pars, x=x)
    out = mod.fit(y, pars, x=x)    
    #comps = out.eval_components(x=x)
    #
    # prepare the output
    aaa = "{},{},{},{},".format(obsid, revol, delta_time,expo1)
    for models in ['g1','g2', 'g3']:
        cen = out.params['%s_center'%models].value
        cen_err = out.params['%s_center'%models].stderr
        fwhm = out.params['%s_fwhm'%models].value
        fwhm_err = out.params['%s_fwhm'%models].stderr
        peak = out.params['%s_amplitude'%models].value/expo1
        peak_err = out.params['%s_amplitude'%models].stderr/expo1
        aaa += "{},{},{},{},{},{},".format(cen,cen_err,fwhm,fwhm_err,peak,peak_err)
    #
    aaa += "%f"%out.chisqr
    #
    if (verbose):
        print(out.fit_report(min_correl=0.5))
    if (plotIt):
        fig = plt.figure(figsize=(10,10))
        ax1 = fig.add_subplot(211)
        ax1.plot(x, y, label="{}_{}".format(revol,obsid))
        #ax1.plot(x, init, 'k--')
        ax1.plot(x, out.best_fit, 'r--')
        #ax1.plot(x, comps['g1_'], 'b--')
        #ax1.plot(x, comps['g2_'], 'b--')
        #ax1.plot(x, comps['g3_'], 'b--')
        #ax1.plot(x, comps['poly_'], 'k--')
        ax1.set_ylabel("Counts")
        ax1.grid(True)
        ax1.legend()
        plt.title("PN %s"%window)
        #
        ax2 = fig.add_subplot(212)
        ax2.plot(x,out.best_fit - y, 'k-')
        ax2.set_ylabel("Residual (counts)")
        ax2.set_xlabel("Channel")
        ax2.grid(True)
        if (savePng):
            plt.savefig(pngName,dpi=100)
            plt.close()
        else:
            plt.show()
    return aaa

#
# 
#%%
wdir = '/xdata/xcaldata/XMM/IVAN/PN_SW/CalClosed'

ver = 't11'

t = Table.read(f'{wdir}/PN_CalClosed_SW.csv',format='ascii.csv')
obsid_list = t['observation_id']
nt = len(obsid_list)
#
#fout = open(f"{wdir}/pn_sw_calclosed_fit_results_nocti.csv","a")
fout = open(f"{wdir}/pn_sw_calclosed_fit_results_cti49_{ver}.csv","w")
print ("obsid,rev,time,exp,al_cen,al_cen_err,al_fwhm,al_fwhm_err,al_peak,al_peak_err,mn1_cen,mn1_cen_err,mn1_fwhm,mn1_fwhm_err,mn1_peak,mn1_peak_err,mn2_cen,mn2_cen_err,mn2_fwhm,mn2_fwhm_err,mn2_peak,mn2_peak_err,chisqr",file=fout)

for i,iobs in enumerate(obsid_list):
    print (f"Fitting CC OBSID {iobs} ({i+1}/{nt})")
    sobs = "%010i"%iobs
    in_spec = f'{wdir}/{sobs}/cti49/spectrum_bin5.fits'
    if (not os.path.isfile(in_spec)):
        print (f"No spectrum found for {sobs}")
        continue
    result = fit_cal_lines(in_spec,smooth=9, savePng=True, pngName=f"{wdir}/spectral_fit/{sobs}_cti49_fitted_{ver}.png")
    print (result,file=fout)
fout.close()
print ("all done")

