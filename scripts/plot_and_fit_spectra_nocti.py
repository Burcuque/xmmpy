#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 18 16:38:44 2018

Plot the extracted spectra and fit a simple model

@author: ivaltchanov
"""

import os
import xspec

import argparse

import numpy as np 

from datetime import datetime

#from astropy.table import Table
from astropy.io import fits
#from astropy.stats import median_absolute_deviation as mad
#from astropy.stats import mad_std

#import matplotlib as mpl
#mpl.use('Agg')

import matplotlib.pylab as plt

#import logging
#
# global folders
#
os.environ['SAS_CCFPATH'] = '/xdata/ccf/pub'
#%%
#
# get the arguments
parser = argparse.ArgumentParser(description='Read, fit and plot spectra, saving the results')
parser.add_argument('obsid', type=str,
                    help='The OBSID to process')
parser.add_argument('--line', type=float, default=6.3,
                    help='Initial line enery in keV')
parser.add_argument('--output', type=str, default="output_xspec_nocti.csv",
                    help='Output file, will be appended')
parser.add_argument('--wdir', type=str, default=os.getcwd(),
                    help='Where the spectral files are')
#
args = parser.parse_args()
#args = parser.parse_args("0555471101 pn --line 6.7 --output output_test.csv --wdir /home/ivaltchanov/XMM/CAL/PN_SW/WR140")
#
#%%
#
wdir = args.wdir
obsid = args.obsid
inst = "pn"
specDir = os.path.join(os.path.abspath(wdir),obsid,"no_cti")
lineIn = args.line
#

#wdir = "/home/ivaltchanov/XMM/CAL/PN_SW/WR140"
#wdir = "/home/ivaltchanov/XMM/CAL/PN_SW/mrk1048"
#inst = "pn"
#obsid = "0690870501"

#lineIn = 6.1
#file_out = os.path.join(wdir,'output_test.csv')
#
if (not os.path.isdir(specDir)):
    print ("The folder {} does not exist! Spectra should be in it. Cannot continue".format(wdir))
    raise FileNotFoundError
#    
file_out = os.path.join(wdir,args.output)
os.chdir(specDir)
#
#
spec_file = "{}_spectrum_nocti_grp0.fits".format(inst)
#spec_file = "{}_spectrum_grp25.fits".format(inst)
if (not os.path.isfile(spec_file)):
    print ("No spectrum file found {}".format(spec_file))
    raise FileNotFoundError
#
try:
    s = xspec.Spectrum(spec_file)
except:
    print ("Cannot read {} in XSPEC".format(spec_file))
    raise Exception
#
# Times will be relative to 2000-01-01T00:00:00
#
time0 = datetime.strptime("2000-01-01T00:00:00","%Y-%m-%dT%H:%M:%S")
hdu = fits.open(spec_file)
target = hdu[0].header['OBJECT']
rev = hdu[0].header['REVOLUT']
submode = hdu[0].header['SUBMODE']
xfilt = hdu[0].header['FILTER']
ontime = hdu[1].header['EXPOSURE']/1000.0 # in ksec
start_time = hdu[0].header['DATE-OBS']

stime = datetime.strptime(start_time,"%Y-%m-%dT%H:%M:%S")
delta_time = (stime-time0).total_seconds()/(365.0*24.0*3600.0) # in years
#
title = "{} for {} ({}_{})".format(inst,target,obsid,rev)
#
logfile = xspec.Xset.openLog("{}_{}_xspec.log".format(target,inst))
#
#
s.notice("all")
xspec.Plot.setRebin(minSig=40.0,maxBins=4)
#%%
#
# XPSEC fitting in [5,8] keV
#
s.ignore("**-5.0 8.0-**")
xspec.Fit.statMethod = "cstat"
xspec.Xset.abund = "wilm"
#
mod1 = xspec.Model("po + ga")
mod1.setPars(2.5,4.0e-2,lineIn,1.0e-2,7.0e-5)
par4 = mod1.gaussian.LineE
par4.values = [lineIn, .01, lineIn-0.5, lineIn-0.5, lineIn+0.5, lineIn+0.5]
mod1.gaussian.Sigma.frozen=True
#
xspec.Fit.nIterations = 100
xspec.Fit.query = 'yes'
xspec.Fit.perform()
xspec.Fit.show()
xspec.Fit.error("2.706 3")
#
# saving the results
#
par3 = mod1(3)
lineE = par3.values[0]
lineE_low = lineE - par3.error[0]
lineE_up = par3.error[1] - lineE
cstat = xspec.Fit.statistic
dof = xspec.Fit.dof
chi2r = xspec.Fit.testStatistic/dof
#
psfile = "{}/spec_fit/{}_{}_nocti_fit.ps".format(wdir,obsid,inst)
xspec.Plot.device = '{}/cps'.format(psfile)
xspec.Plot.xAxis = "keV"
xspec.Plot.setRebin(minSig=30.0,maxBins=4)
xspec.Plot.addCommand("csize 1.2")
xspec.Plot.addCommand("lwidth 4")
xspec.Plot.addCommand("lab top \"{}\"".format(title))
#xspec.Plot()    
xspec.Plot('data','ratio') 
xspec.Plot.device = '/xs'
xspec.Plot()
#
with open(file_out,'a') as fout:
    full = 0
    if ('Full' in submode or 'Large' in submode):
        full = 1
    output = "{},{},{},{},{},{},{},{:.2f},{:.3f},{:.3f},{:.3f},{:.1f},{:.4f},{}".format(\
          obsid,rev,delta_time,submode,full,xfilt,inst,ontime,lineE,lineE_low,lineE_up,cstat,chi2r,dof)
    print (output,file=fout)
#
#%%
xspec.AllData.clear()
xspec.Xset.closeLog()
