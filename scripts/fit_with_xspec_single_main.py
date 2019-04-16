#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 29 12:30:37 2019

Function to fit an input spectrum with XSPEC

Each line is individually fit

@author: ivaltchanov
"""

import os
import xspec

from datetime import datetime

from astropy.io import fits

import argparse

#%%
# get the arguments
parser = argparse.ArgumentParser(description='Read, fit and plot CalClosed spectra for PNLW')
parser.add_argument('spec_file', type=str,
                    help='The spectrum file to process with XSPEC')
parser.add_argument('--output', type=str, default="output_xspec.csv",
                    help='Output file, will be appended')
parser.add_argument('--plot_file', type=str, default="",
                    help='XSPEC plot file')
#
args = parser.parse_args()
#

#%%
wdir = os.getcwd()

spec_file = args.spec_file
file_out = args.output
psfile = args.plot_file
#spec_file = f"{wdir}/0727360201/PN_LW/pn_S003_06_all_spec5.fits"
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
#file_out = os.path.join(wdir,f"{obsid}_{expo}_xspec_output.csv")
#
#s.response.rmf = f"{wdir}/epn_e3_lw20_sY9_v17.0.rmf"

fout = open(file_out,'a')
output = f"{obsid},{expo},{rev},{delta_time:.4f},{submode},{xfilt},{inst},{ontime:.1f},{ccd},"
#          par3.values[0],p3err,\
#          par6.values[0],p6err,\
#          par9.values[0],p9err,\
#          par12.values[0],p12err,\
#          cstat,chi2r,dof)
#    print (output,file=fout)

#
xspec.Fit.statMethod = "cstat"
xspec.Xset.abund = "wilm"

#%%
lines = {}
lines['AlKa'] = {'cen': 1.468, 'range': [1.0,2.0]}
lines['MnKa'] = {'cen': 5.8988, 'range': [5.0,7.0]}
lines['CuKa'] = {'cen': 8.04, 'range': [7.0,9.0]}

for line in lines.keys():
    s.notice("all")
    s.ignore("**-{} {}-**".format(lines[line]['range'][0],lines[line]['range'][1]))
    #
    model = xspec.Model("po + ga")
    model.gaussian.Sigma = [5.0e-2,0.001,1.0e-3,1.0e-3,1.0e-1,1.0e-1]
    line_cen = lines[line]['cen']
    model.gaussian.LineE = [line_cen,0.001,line_cen-0.2,line_cen-0.1,line_cen+0.1,line_cen+0.2]
    #
    # will also fit MnKb
    if (line == 'MnKa'):
        model = xspec.Model("po + ga + ga")
        model.gaussian.Sigma = [5.0e-2,0.001,1.0e-3,1.0e-3,1.0e-1,1.0e-1]
        model.gaussian_3.Sigma = [5.0e-2,0.001,1.0e-3,1.0e-3,1.0e-1,1.0e-1]
        line_cen = lines[line]['cen']
        model.gaussian.LineE = [line_cen,0.001,line_cen-0.2,line_cen-0.1,line_cen+0.1,line_cen+0.2]
        model.gaussian_3.LineE = [6.49,0.001,6.49-0.2,6.49-0.1,6.49+0.1,6.49+0.2]
        
    # freeze the sigmas
    #model.gaussian.Sigma.frozen = True
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
    if (line == 'MnKa'):
        xspec.Fit.error("2.706 3 6")
        par3 = model.gaussian.LineE
        p3err = (par3.error[1] - par3.error[0])/2.0
        par6 = model.gaussian_3.LineE
        p6err = (par6.error[1] - par6.error[0])/2.0
        output += f"{par3.values[0]:.4f},{p3err:.4f},{par6.values[0]:.4f},{p6err:.4f},{cstat:.2f},{chi2r:.3f},{dof},"
    else:
        xspec.Fit.error("2.706 3")
        par3 = model.gaussian.LineE
        p3err = (par3.error[1] - par3.error[0])/2.0
        output += f"{par3.values[0]:.4f},{p3err:.4f},{cstat:.2f},{chi2r:.3f},{dof},"
    #
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
    psfile = f"{wdir}/spectral_fit/pn_{obsid}_{expo}_{ccd}_{line}_xspec.ps"
    if (psfile != ""):
        xspec.Plot.device = '{}/cps'.format(psfile)
        xspec.Plot("data","ratio") 
    #
    xspec.Plot.device = '/xs'
    xspec.Plot("data","ratio") 
    #
    #v = input('Continue? (y/n) ')
    
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
##
#%%
xspec.Plot.commands = ()
xspec.AllData.clear()
xspec.Xset.closeLog()
print (output,file=fout)
fout.close()