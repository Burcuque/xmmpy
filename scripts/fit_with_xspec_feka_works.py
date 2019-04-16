#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 29 12:30:37 2019

Function to fit an input spectrum for the Fe Kalpha line with XSPEC

@author: ivaltchanov
"""

import os
import xspec

from datetime import datetime

from astropy.io import fits

import argparse

#%%
# get the arguments
parser = argparse.ArgumentParser(description='Read, fit FeKa line and plot for PNLW')
parser.add_argument('spec_file', type=str,
                    help='The spectrum file to process with XSPEC')
parser.add_argument('--output', type=str, default="output_xspec.csv",
                    help='Output file, will be appended')
parser.add_argument('--plot_file', type=str, default="",
                    help='XSPEC plot file')
#
# check with one spectrum
#
#args = parser.parse_args(["ngc6814/0764620101/nocti/pn_spectrum_grp0.fits"])
args = parser.parse_args()
#

#%%

wdir = os.path.join(os.path.expanduser('~'),'IVAN','PN_LW','sources')
os.chdir(wdir)

spec_file = args.spec_file
file_out = f"{wdir}/{args.output}"
psfile = args.plot_file
#spec_file = f"{wdir}/0727360201/PN_LW/pn_S003_06_all_spec5.fits"
#spec_file = f"{wdir}/0727360201/PN_LW/pn_S003_04_180_spec5.fits"
#
if (not os.path.isfile(spec_file)):
    print (f"No spectrum file found {spec_file}")
    raise FileNotFoundError
#
# Need to change to the folder with the spectrum in order to find the RMF and ARF files
#
spec_dir = os.path.dirname(spec_file)
os.chdir(spec_dir)
#
hdu = fits.open(os.path.basename(spec_file))
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
#ccd = xfile.split('_')[2]
#
title = f"{target} ({rev}_{obsid}_{expo})"
#title = f"{inst} in {submode} for {target} ({rev}_{obsid}_{expo})"
hdu.close()
#        
try:
    s = xspec.Spectrum(os.path.basename(spec_file))
except:
    print (f"Cannot read {spec_file} in XSPEC")
    raise Exception
#
# Times will be relative to 2000-01-01T00:00:00
#
#
logfile = xspec.Xset.openLog(f"xspec.log")
#
#file_out = os.path.join(wdir,f"{obsid}_{expo}_xspec_output.csv")
#
#s.response.rmf = f"{wdir}/epn_e3_lw20_sY9_v17.0.rmf"
s.notice("all")
#
# XPSEC fitting in [4,8] keV
#
s.ignore("**-4.0 8.0-**")
xspec.Fit.statMethod = "cstat"
xspec.Xset.abund = "wilm"

#%%epn_e3_lw20_sY9_v17.0.rmf
#
# For PN Large Window CalClosed, will fit 4 lines:
#    Al Ka (1.5 keV), Mn Ka (5.9 keV), Mn Kb (6.4 keV) and Cu Ka (8.0 keV)
#
model = xspec.Model("po + ga")
#
#
model.gaussian.Sigma = [5.0e-2,0.001,1.0e-3,1.0e-3,1.0e-1,1.0e-1]
#
model.gaussian.LineE = [6.4,0.001,6.1,6.1,6.5,6.5]
#
xspec.Fit.nIterations = 100
xspec.Fit.query = 'yes'
xspec.Fit.perform()
xspec.Fit.show()
cstat = xspec.Fit.statistic
dof = xspec.Fit.dof
chi2r = xspec.Fit.testStatistic/dof
#
# calculate the 1-sigma errors
#
xspec.Fit.error("2.706 3")
#xspec.Fit.error("3")
par3 = model.gaussian.LineE
p3err = (par3.error[1] - par3.error[0])/2.0
p3sig = par3.sigma # the fit error (one parameter)
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
psfile_out = f"{wdir}/{psfile}"
if (psfile != ""):
    print (f"Saving the plot to {psfile} ({psfile_out})")
    xspec.Plot.device = f'{psfile_out}/cps'
    xspec.Plot("data","ratio") 
#
xspec.Plot.device = '/xs'
xspec.Plot("data","ratio") 
#
#
with open(file_out,'a') as fout:
    output = "{},{},{},{:.4f},{},{},{},{:.1f},{:.4f},{:.4f},{:.4f},{:.1f},{:.4f},{}".format(\
          obsid,expo,rev,delta_time,submode,xfilt,inst,ontime,\
          par3.values[0],p3err,p3sig,\
          cstat,chi2r,dof)
    print (output,file=fout)
#
#%%
xspec.Plot.commands = ()
xspec.AllData.clear()
xspec.Xset.closeLog()
