#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 29 12:30:37 2019

Function to fit an input spectrum with XSPEC

@author: ivaltchanov
"""

import os
import xspec

import numpy as np 
from datetime import datetime

from astropy.io import fits

import matplotlib.pylab as plt

#%%
wdir = os.path.join(os.path.expanduser('~'),"IVAN/PN_LW")

spec_file = f"{wdir}/0727360201/PN_LW/pn_S003_06_all_spec5.fits"
#spec_file = f"{wdir}/0727360201/PN_LW/pn_S003_04_180_spec5.fits"
#
if (not os.path.isfile(spec_file)):
    print (f"No spectrum file found {spec_file}")
    raise FileNotFoundError
#
hdu = fits.open(spec_file)
time0 = datetime.strptime("2000-01-01T00:00:00","%Y-%m-%dT%H:%M:%S")
obsid = hdu[0].header['OBS_ID']
expo = hdu[0].header['EXPIDSTR']
target = hdu[0].header['OBJECT']
rev = hdu[0].header['REVOLUT']
submode = hdu[0].header['SUBMODE']
xfilt = hdu[0].header['FILTER']
inst = hdu[0].header['INSTRUME']
ontime = hdu[1].header['EXPOSURE']/1000.0 # in ksec
start_time = hdu[0].header['DATE-OBS']

stime = datetime.strptime(start_time,"%Y-%m-%dT%H:%M:%S")
delta_time = (stime-time0).total_seconds()/(365.0*24.0*3600.0) # in years
#
xfile = os.path.basename(spec_file)
ccd = xfile.split('_')[2]
#
title = f"{xfile} ({rev}_{obsid}_{expo})"
#title = f"{inst} in {submode} for {target} ({rev}_{obsid}_{expo})"
hdu.close()
#        
try:
    s = xspec.Spectrum(spec_file)
except:
    print (f"Cannot read {spec_file} in XSPEC")
    raise Exception
#
# Times will be relative to 2000-01-01T00:00:00
#
#
logfile = xspec.Xset.openLog(f"xspec.log")
#
file_out = os.path.join(wdir,f"{obsid}_{expo}_xspec_output.csv")
#
#s.response.rmf = f"{wdir}/epn_e3_lw20_sY9_v17.0.rmf"
s.notice("all")
#
# XPSEC fitting in [1,10] keV
#
s.ignore("**-1.0 10.0-**")
xspec.Fit.statMethod = "cstat"
xspec.Xset.abund = "wilm"

#%%epn_e3_lw20_sY9_v17.0.rmf
#
# For PN Large Window CalClosed, will fit 4 lines:
#    Al Ka (1.5 keV), Mn Ka (5.9 keV), Mn Kb (6.4 keV) and Cu Ka (8.0 keV)
#
model = xspec.Model("po + ga + ga + ga + ga",\
                    setPars={3:1.5, 4:1.0e-2, 6:5.9, 5:1.0e-2, 9:6.4, 10:1.0e-2, 12: 8.0, 13:1.0e-2})
#
# freeze the sigmas
#
model.gaussian.Sigma.frozen = True
model.gaussian_3.Sigma.frozen = True
model.gaussian_4.Sigma.frozen = True
model.gaussian_5.Sigma.frozen = True
#
xspec.Fit.nIterations = 100
xspec.Fit.query = 'yes'
xspec.Fit.perform()
xspec.Fit.show()
#
# saving the results
#
#par3 = model(3)
#lineE = par3.values[0]
#lineE_low = lineE - par3.error[0]
#lineE_up = par3.error[1] - lineE
cstat = xspec.Fit.statistic
dof = xspec.Fit.dof
chi2r = xspec.Fit.testStatistic/dof
xspec.Fit.error("2.706 3 6 9 12")
#
par3 = model.gaussian.LineE
par6 = model.gaussian_3.LineE
par9 = model.gaussian_4.LineE
par12 = model.gaussian_5.LineE

p3err = (par3.error[1] - par3.error[0])/2.0
p6err = (par6.error[1] - par6.error[0])/2.0
p9err = (par9.error[1] - par9.error[0])/2.0
p12err = (par12.error[1] - par12.error[0])/2.0
#
#psfile = "{}/spec_fit/{}_{}_fit_cti49.ps".format(wdir,obsid,inst)
#xspec.Plot.device = '{}/cps'.format(psfile)
xspec.Plot.xAxis = "keV"
xspec.Plot.setRebin(minSig=30.0,maxBins=10)
xspec.Plot.addCommand("csize 1.2")
xspec.Plot.addCommand("lwidth 2 on 2")
xspec.Plot.addCommand("lab top \"{}\"".format(title))
xspec.Plot.addCommand("color 1 on 1")
xspec.Plot.addCommand("color 2 on 2")
xspec.Plot.addCommand("color 1 on 3")
xspec.Plot.addCommand("color 3 on 4")
#xspec.Plot()
psfile = f"{wdir}/spectral_fit/pn_{obsid}_{expo}_{ccd}_xspec.ps"
xspec.Plot.device = '{}/cps'.format(psfile)
xspec.Plot("data","ratio") 
#
#
xspec.Plot.device = '/xs'
xspec.Plot() 
#
#xspec.Plot.addCommand("vp 0.15 0.15 0.9 0.9")
#
# output the XSPEC data
#xvals = xspec.Plot.x()
#yvals = xspec.Plot.y()
#xErrs = xspec.Plot.xErr()
#yErrs = xspec.Plot.yErr()
#modvals = xspec.Plot.model()
#
#tout = Table([xvals, xErrs,yvals,yErrs,modvals], names=('energy', 'energyErr', 'norm','normErr','model'))
#tout.write(f"{wdir}/{obsid}_xspec_dump.csv",overwrite=True)
#
with open(file_out,'a') as fout:
    output = "{},{},{},{:.4f},{},{},{},{:.1f},{},{:.3f},{:.3f},{:.3f},{:.3f},{:.3f},{:.3f},{:.3f},{:.3f},{:.1f},{:.4f},{}".format(\
          obsid,expo,rev,delta_time,submode,xfilt,inst,ontime,ccd,\
          par3.values[0],p3err,\
          par6.values[0],p6err,\
          par9.values[0],p9err,\
          par12.values[0],p12err,\
          cstat,chi2r,dof)
    print (output,file=fout)
#
#%%
xspec.Plot.commands = ()
xspec.AllData.clear()
xspec.Xset.closeLog()
