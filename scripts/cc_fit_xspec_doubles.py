#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 15 18:23:44 2019

Fit the CalClosed lines with XSPEC, using spectra with only double events

@author: ivaltchanov
"""

import os
import xspec

import argparse

#import numpy as np 
from datetime import datetime

#from astropy.table import Table
from astropy.io import fits
#from astropy.stats import median_absolute_deviation as mad
#from astropy.stats import mad_std

#import matplotlib as mpl
#mpl.use('Agg')

#import matplotlib.pylab as plt

#import logging
#
# global folders
#
os.environ['SAS_CCFPATH'] = '/xdata/ccf/pub'
#%%
#
# get the arguments
parser = argparse.ArgumentParser(description='Read, fit and plot spectra for PNSW in CC')
parser.add_argument('obsid', type=str,
                    help='The OBSID to process')
parser.add_argument('subfolder', type=str, default="no_cti",
                    help='The variant \"no_cti\" or \"cti49\"')
parser.add_argument('--output', type=str, default="output_xspec_doubles.csv",
                    help='Output file, will be appended')
parser.add_argument('--wdir', type=str, default=os.getcwd(),
                    help='Where the ODF files are')
#
args = parser.parse_args()
#
#%%
#
wdir = os.path.abspath(args.wdir)
inst = "pn"
obsid = args.obsid
var = args.subfolder
file_out = os.path.join(wdir,args.output)
#
odfDir = os.path.join(wdir,obsid)

if (not odfDir):
    print ("The ODF folder {} does not exist! Cannot continue".format(odfDir))
    raise FileNotFoundError
#    
ppsDir = os.path.join(wdir,obsid,var)
#
if (not os.path.isdir(ppsDir)):
    print ("PPS folder {} does not exist.".format(ppsDir))
    raise FileNotFoundError

os.chdir(ppsDir)
#
# we'll use th eunbinned spectrum for fitting in XSPEC, using C-stat
# the binned version is used only for plotting
#
spec_file = 'pn_spectrum_doubles_grp0.fits'
#
try:
    s = xspec.Spectrum(spec_file)
except:
    print ("Cannot read {} in XSPEC".format(spec_file))
    raise Exception

#spec_file = "{}_spectrum_grp25.fits".format(inst)
#
# Times will be relative to 2000-01-01T00:00:00
#
time0 = datetime.strptime("2000-01-01T00:00:00","%Y-%m-%dT%H:%M:%S")
hdu = fits.open(spec_file)
rev = hdu[0].header['REVOLUT']
submode = hdu[0].header['SUBMODE']
xfilt = hdu[0].header['FILTER']
ontime = hdu[1].header['EXPOSURE']/1000.0 # in ksec
start_time = hdu[0].header['DATE-OBS']

stime = datetime.strptime(start_time,"%Y-%m-%dT%H:%M:%S")
delta_time = (stime-time0).total_seconds()/(365.0*24.0*3600.0) # in years
#
title = "PN-SW CalClosed ({}_{})".format(obsid,rev)
#
logfile = xspec.Xset.openLog(f"cc_xspec_doubles_{var}.log")
#
#file_out = os.path.join(wdir,args.output)
#file_out = os.path.join(wdir,f'cc_output_xspec_{var}.csv')
#
#%%
#
# now fit with Xspec
#
#s = xspec.Spectrum('{}_src_spectrum.fits'.format(inst))
#s.background = '{}_bkg_spectrum.fits'.format(inst)
#s.response.arf = "{}.arf".format(inst)
#s.response.rmf = "{}.rmf".format(inst)
s.notice("all")
#s.ignore("bad")
xspec.Plot.xAxis = 'keV'
s.ignore("**-1.0 8.0-**")
#xspec.Fit.statMethod = "cstat"
#xspec.Xset.abund = "wilm"
#
mod2 = xspec.Model("power + gauss + gauss + gauss")
#
par1 = mod2.powerlaw.PhoIndex
par1.values = 1.1
par2 = mod2.powerlaw.norm
par2.values = 1.0
#
par3 = mod2.gaussian.LineE
par3.values = [1.486,0.01,1.45,1.45,1.53,1.53]
par4 = mod2.gaussian.Sigma
par4.values = 1.0e-3
par4.frozen = True
par5 = mod2.gaussian.norm
par5.values = 3.0e-4
#
par6 = mod2.gaussian_3.LineE
par6.values = [5.8988,0.01,5.7,5.7,6.1,6.1]
par7 = mod2.gaussian_3.Sigma
par7.values = 1.0e-3
par7.frozen = True
par8 = mod2.gaussian_3.norm
par8.values = 3.0e-5
#
par9 = mod2.gaussian_4.LineE
par9.values = [6.49,0.01,6.3,6.3,6.6,6.6]
par10 = mod2.gaussian_4.Sigma
par10.values = 1.0e-3
par10.frozen = True
par11 = mod2.gaussian_4.norm
par11.values = 3.0e-5
#
xspec.Fit.nIterations = 100
xspec.Fit.query = 'yes'
xspec.Fit.perform()
xspec.Fit.show()
#
xspec.Fit.error("2.706 3 6 9")
#
# saving the results
#
cstat = xspec.Fit.statistic
dof = xspec.Fit.dof
chi2r = xspec.Fit.testStatistic/dof
#
psfile = f"{wdir}/spectral_fit/pn_{obsid}_xspec_doubles_{var}.ps"
xspec.Plot.device = '{}/cps'.format(psfile)
#xspec.Plot.device = '/xs'
xspec.Plot.setRebin(minSig=5.0,maxBins=10,groupNum=1)
xspec.Plot.addCommand("csize 1.2")
xspec.Plot.addCommand("lwidth 3 on 1")
xspec.Plot.addCommand("lab top \"{}\"".format(title))
#xspec.Plot.addCommand("rescale y 0.1 3.0")
xspec.Plot.addCommand("color 1 on 1")
xspec.Plot.addCommand("color 2 on 2")
xspec.Plot.addCommand("color 1 on 3")
#xspec.Plot()    
xspec.Plot('data','ratio') 
xspec.Plot.device = '/xs'
xspec.Plot()
xspec.AllData.clear()
xspec.Xset.closeLog()
#
p3error = (par3.error[1] - par3.error[0])/2.0
p6error = (par6.error[1] - par6.error[0])/2.0
p9error = (par9.error[1] - par9.error[0])/2.0
#

with open(file_out,'a') as fout:
    full = 0
    if ('Full' in submode or 'Large' in submode):
        full = 1
    output = "{},{},{:.5f},{},{},{},{},{:.2f},{:.4f},{:.4f},{:.4f},{:.4f},{:.4f},{:.4f},{:.1f},{:.3f},{}".format(\
          obsid,rev,delta_time,submode,full,xfilt,inst,ontime,\
          par3.values[0],p3error,\
          par6.values[0],p6error,\
          par9.values[0],p9error,\
          cstat,chi2r,dof)
    print (output,file=fout)
#
print (f"All done, results saved to {file_out}")