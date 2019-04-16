#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 15 18:23:44 2019

Fit the CalClosed lines with XSPEC

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
parser.add_argument('--output', type=str, default="output_xspec_fit2_cti49.csv",
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
#
odfDir = os.path.join(wdir,obsid)

if (not odfDir):
    print ("The ODF folder {} does not exist! Cannot continue".format(odfDir))
    raise FileNotFoundError
#    
ppsDir = os.path.join(wdir,obsid,'no_cti')
if (not os.path.isdir(ppsDir)):
    print ("PPS folder {} does not exist.".format(ppsDir))
    raise FileNotFoundError

os.chdir(ppsDir)
#
# we'll use th eunbinned spectrum for fitting in XSPEC, using C-stat
# the binned version is used only for plotting
#
spec_file = "spectrum_bin5_grp0.fits"
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
logfile = xspec.Xset.openLog("cc_xspec_nocti.log")
#
#file_out = os.path.join(wdir,args.output)
file_out = os.path.join(wdir,'cc_output_xspec_nocti.csv')
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
mod2 = xspec.Model("power + gauss + gauss + gauss)")
#
par1 = mod2.powerlaw.PhoIndex
par1.values = 1.5
par1 = mod2.powerlaw.norm
par1.values = 1.0e-3
#
par2 = mod2.gaussian.LineE
par2.values = [1.486,0.01,1.48,1.48,1.49,1.49]
par2 = mod2.gaussian.Sigma
par2.values = 1.0e-3
par2.frozen = True
par2 = mod2.gaussian.norm
par2.values = 3.0e-4
#
par3 = mod2.gaussian_3.LineE
par3.values = [5.8988,0.01,5.8,5.8,6.0,6.0]
par3 = mod2.gaussian_3.Sigma
par3.values = 1.0e-3
par3.frozen = True
par3 = mod2.gaussian_3.norm
par3.values = 3.0e-5
#
par4 = mod2.gaussian_4.LineE
par4.values = [6.49,0.01,6.45,6.45,6.55,6.55]
par4 = mod2.gaussian_4.Sigma
par4.values = 1.0e-3
par4.frozen = True
par4 = mod2.gaussian_4.norm
par4.values = 3.0e-5
#
xspec.Fit.nIterations = 100
xspec.Fit.query = 'yes'
xspec.Fit.perform()
xspec.Fit.show()
#xspec.Fit.error("2.706 6")
#
# saving the results
#
#lineE = par6.values[0]
#lineE_low = lineE - par6.error[0]
#lineE_up = par6.error[1] - lineE
cstat = xspec.Fit.statistic
dof = xspec.Fit.dof
chi2r = xspec.Fit.testStatistic/dof
#
#psfile = "{}/spec_fit/{}_{}_fit2.ps".format(wdir,obsid,inst)
#xspec.Plot.device = '{}/cps'.format(psfile)
#xspec.Plot.device = '/xs'
xspec.Plot.setRebin(minSig=30.0,maxBins=4)
xspec.Plot.addCommand("csize 1.2")
xspec.Plot.addCommand("lwidth 3 on 1")
xspec.Plot.addCommand("lab top \"{}\"".format(title))
xspec.Plot.addCommand("rescale y 0.1 3.0")
xspec.Plot.addCommand("color 1 on 1")
xspec.Plot.addCommand("color 1 on 3")
#xspec.Plot()    
xspec.Plot('ldata','ratio') 
xspec.Plot.device = '/xs'
xspec.Plot()
xspec.AllData.clear()
xspec.Xset.closeLog()
