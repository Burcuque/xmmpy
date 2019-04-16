#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 18 16:38:44 2018

Pot the extracted spectra and fit a complex model,
applicable for NGC 4151


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
parser = argparse.ArgumentParser(description='Read, fit and plot spectra for NGC4151 with a complex model, saving the results')
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
ppsDir = os.path.join(wdir,obsid,'cti49')
if (not os.path.isdir(ppsDir)):
    print ("PPS folder {} does not exist.".format(ppsDir))
    raise FileNotFoundError

os.chdir(ppsDir)
#
# we'll use th eunbinned spectrum for fitting in XSPEC, using C-stat
# the binned version is used only for plotting
#
spec_file = "{}_spectrum_cti49_grp0.fits".format(inst)
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
target = hdu[0].header['OBJECT']
target = "".join(target.split()) # remove whitespace in name
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
logfile = xspec.Xset.openLog("{}_{}_xspec_fit2_cti49.log".format(target,inst))
#
#file_out = os.path.join(wdir,args.output)
file_out = os.path.join(wdir,'output_xspec_cti49_fit2.csv')
#
if (not os.path.isfile(spec_file)):
    print ("No spectrum file found {}".format(spec_file))
    raise FileNotFoundError
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
s.ignore("**-5.0 8.0-**")
xspec.Fit.statMethod = "cstat"
xspec.Xset.abund = "wilm"
#
mod2 = xspec.Model("phabs*(power + bbody + gauss + gauss)")
#
#mod2.setPars(36.0,1.3,2.5,4.0e-2,6.3,1.0e-2,7.0e-5)
par1 = mod2.phabs.nH
par1.values = 30.0 # 10^22 cm^-2
#
par2 = mod2.powerlaw.PhoIndex
par2.values = 1.5
par3 = mod2.powerlaw.norm
par3.values = 1.0e-3
#
par4 = mod2.bbody.kT
par4.values = 2.0
par5 = mod2.bbody.norm
par5.values = 3.0e-3
#
par6 = mod2.gaussian.LineE
par6.values = [6.3,0.01,6.1,6.1,6.5,6.5]
par7 = mod2.gaussian.Sigma
par7.values = 1.0e-3
par7.frozen = True
par8 = mod2.gaussian.norm
par8.values = 3.0e-4
#
par9 = mod2.gaussian_5.LineE
par9.values = [6.9,0.01,6.8,6.8,7.2,7.2]
par10 = mod2.gaussian_5.Sigma
par10.values = 1.0e-3
par10.frozen = True
par11 = mod2.gaussian_5.norm
par11.values = 3.0e-5
#
xspec.Fit.nIterations = 100
xspec.Fit.query = 'yes'
xspec.Fit.perform()
xspec.Fit.show()
xspec.Fit.error("2.706 6")
#
# saving the results
#
lineE = par6.values[0]
lineE_low = lineE - par6.error[0]
lineE_up = par6.error[1] - lineE
cstat = xspec.Fit.statistic
dof = xspec.Fit.dof
chi2r = xspec.Fit.testStatistic/dof
#
psfile = "{}/spec_fit/{}_{}_fit2.ps".format(wdir,obsid,inst)
xspec.Plot.device = '{}/cps'.format(psfile)
#xspec.Plot.device = '/xs'
xspec.Plot.setRebin(minSig=30.0,maxBins=4)
xspec.Plot.addCommand("csize 1.2")
xspec.Plot.addCommand("lwidth 3 on 1")
xspec.Plot.addCommand("lab top \"{}\"".format(title))
if (inst == 'pn'):
    xspec.Plot.addCommand("rescale y 0.1 3.0")
else:
    xspec.Plot.addCommand("rescale y 0.02 0.7")
xspec.Plot.addCommand("color 1 on 1")
xspec.Plot.addCommand("color 1 on 3")
#xspec.Plot()    
xspec.Plot('ldata','ratio') 
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
xspec.AllData.clear()
xspec.Xset.closeLog()
