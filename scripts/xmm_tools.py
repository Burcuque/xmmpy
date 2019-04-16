#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 11 11:09:13 2019

@author: ivaltchanov
"""
import numpy as np
from astropy.table import Table
from astropy.io import fits
from astropy.convolution import convolve, Box1DKernel
from lmfit.models import GaussianModel, PolynomialModel
from datetime import datetime

import matplotlib.pylab as plt

import subprocess
import sys
import logging

logging.basicConfig(level=logging.DEBUG,
                    format='%(asctime)s %(levelname)s %(message)s',
                    filename='output.log',
                    filemode='w')

def run_command(command,verbose=True):
    #
    # Execute a shell command with the stdout and stderr being redirected to a log file 
    #
    try:
        result = subprocess.run(command, shell=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
        retcode=result.returncode
        if retcode < 0:
            if (verbose):
                print(f"Execution of {command} was terminated by signal", -retcode, file=sys.stderr)
            logging.warning("Execution of {} was terminated by signal: {} \n {}".format(command,-retcode,result.stdout.decode()))
        else:
            if (verbose):
                print(f"Execution of {command} returned", retcode, file=sys.stderr)
            logging.info("Execution of {} returned {}, \n {}".format(command,retcode,result.stdout.decode()))
    except OSError as e:
        print(f"Execution of {command} failed:", e, file=sys.stderr)
        logging.error("Execution of {} failed: {}".format(command,e))
    return retcode

def read_pn_obstable():
    #
    # load the observations from EPICMON account
    #
    
    t = Table.read('/home/epicmon/monitoring/mjss/BadPix/db/pn_all.dat',delimiter='|',
                   format='ascii.no_header',data_start=1,
                   names=('rev','ac_obs_id','ac_exp_id','exp_start_time','exp_end_time','perf_dur',
                          'mode','filter','target','s_obs_id','s_exp_id','ro','rx','rs'))
    return t
#
#
def fit_cu_line(xin,yin, minMax=(7.0,9.0),line_c=8.0):
    #
    # PURPOSE:
    #   Fit the Cu Ka line (8.04 keV), the model is a polynomial(2) + a Gaussian line.
    #
    # INPUTS:
    #   xin is the energy channel (in keV)
    #   yin is the counts
    #   line_c is the initial energy of the line (in keV)
    #   minMax is a tuple (in keV) of the energy limits to consider for the model and the fit
    #
    # OUTPUTS:
    #  a tuple of the full fit output class and the results line in ascii.
    #
    # NOTES:
    #   the Gaussian sigma of the line is only allowed within a certain range: 80 to 250 eV
    #
    #
    m1 = xin >= minMax[0]
    m2 = xin <= minMax[1]
    xw = xin[m1*m2]
    yw = yin[m1*m2]
    i1max = np.argmax(yw)
    y1max = yw[i1max]
    #x1max = xw[i1max]
    #
    poly_mod = PolynomialModel(1,prefix='poly_')
    pars = poly_mod.guess(yw, x=xw)
    #
    gauss1  = GaussianModel(prefix='g1_')
    pars.update( gauss1.make_params())
    pars['g1_center'].set(line_c,min=line_c-0.2,max=line_c+0.1)
    pars['g1_sigma'].set(0.08,min=0.04,max=0.250)
    pars['g1_amplitude'].set(y1max)
    #
    mod = poly_mod + gauss1
    #init = mod.eval(pars, x=x)
    out = mod.fit(yw, pars, x=xw)
    #
    cen = out.params['g1_center'].value
    cen_err = out.params['g1_center'].stderr
    fwhm = out.params['g1_fwhm'].value
    fwhm_err = out.params['g1_fwhm'].stderr
    #peak = out.params['g1_amplitude'].value/expo1
    #peak_err = out.params['g1_amplitude'].stderr/expo1
    chi2 = out.chisqr
    df = len(xw)
    try:
        results  = f"{cen:.3f},{cen_err:.3f},{fwhm:.5f},{fwhm_err:.5f},{chi2:.3f},{df}"
    except:
        results = None
    #
    return (out,results)

def fit_cu_line2(fitsfile, minMax=(7000,9000), plotIt=True, savePng=False, \
    pngName = "fitted.png", verbose=False, smooth=None):
    #
    # PURPOSE:
    #   Fit the Cu Ka line (8.04 keV), the model is a polynomial(2) + a Gaussian line.
    #
    # INPUTS:
    #   fitsfile: the specrtum FITS file, an output from XMM-SAS evselect. 
    #           Must have extension 'SPECTRUM' and columns 'CHANNEL' and 'COUNTS'.
    #   minMax is a tuple (in eV) of the energy limits to consider for the model and the fit
    #   plotIt: if True will plot the results
    #   savePng: if True will save the plot to a PNG file
    #   pgnName: the name of the PNG file, default fitted.png in the current working directory
    #   verbose: if True will print some verbose information
    #   smooth: if an integer number then it will be used for Boxcar 1D smoothing box. 
    #       If None then the fitting will be performed on the raw input spectrum (no smoothing).
    #
    # OUTPUTS:
    #  a tuple of the full fit output class and the results line in ascii.
    #
    # NOTES:
    #   the Gaussian sigma of the line is only allowed within a certain range: 80 to 250 eV
    #
    #
    # read a spectrum and fit for the Cu Ka line at 8.04 keV
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
        print ("Warning: cannot find the exposure, setting it to 1.0")
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
    i1max = np.argmax(y)
    y1max = y[i1max]
    x1max = x[i1max]
    #
    poly_mod = PolynomialModel(1,prefix='poly_')
    pars = poly_mod.guess(y, x=x)
    #
    gauss1  = GaussianModel(prefix='g1_')
    pars.update( gauss1.make_params())
    pars['g1_center'].set(x1max,min=x1max-100.0,max=x1max+100.0)
    pars['g1_sigma'].set(80.0,min=40,max=250)
    pars['g1_amplitude'].set(y1max)
    #
    mod = poly_mod + gauss1
    #init = mod.eval(pars, x=x)
    out = mod.fit(y, pars, x=x)    
    #comps = out.eval_components(x=x)
    #
    # prepare the output
    cen = out.params['g1_center'].value
    cen_err = out.params['g1_center'].stderr
    fwhm = out.params['g1_fwhm'].value
    fwhm_err = out.params['g1_fwhm'].stderr
    peak = out.params['g1_amplitude'].value/expo1
    peak_err = out.params['g1_amplitude'].stderr/expo1
    chi2 = out.chisqr
    df = len(x)
    aaa  = f"{obsid},{revol},{delta_time},{expo1},{cen},{cen_err},{fwhm},{fwhm_err},{peak},{peak_err},{chi2},{df}"
    #
    if (verbose):
        print(out.fit_report(min_correl=0.5))
    if (plotIt):
        fig = plt.figure(figsize=(10,10))
        ax1 = fig.add_subplot(211)
        ax1.plot(x, y, label="{}_{}".format(revol,obsid))
        #ax1.plot(x, init, 'k--')
        ax1.plot(x, out.best_fit, 'r--')
        ax1.axvline(8040.0, color='red',ls='dashed')
        #ax1.plot(x, comps['g1_'], 'b--')
        #ax1.plot(x, comps['g2_'], 'b--')
        #ax1.plot(x, comps['g3_'], 'b--')
        #ax1.plot(x, comps['poly_'], 'k--')
        ax1.set_ylabel("Counts")
        ax1.grid(True)
        ax1.legend()
        plt.title(f"PN {window}")
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
    return (out,aaa)

def read_cti_ccf(ccf_file,mode_id=2, energy_index=0):
    #
    # Read the EPN_CTI CCF file
    #
    #
    # mode_id = 0 for PRIME_FULL_WINDOW
    # mode_id = 1 for PRIME_FULL_WINDOW_EXTENDED
    # mode_id = 2 for PRIME_LARGE_WINDOW
    # mode_id = 3 for PRIME_SMALL_WINDOW
    #
    # energy_index = 0 is for Al Ka at 1.486
    # energy_index = 1 is for Mn Ka at 5.8988
    # energy_index = 2 is for Fe Ka at 6.4 keV
    #
    try:
        t = Table.read(f'{ccf_file}',hdu='LONG_TERM_CTI')
        tmp = Table.read(f'{ccf_file}',hdu='LTC_TIMES')
        #times = tmp.data['TIME'][0].flatten()                           
    except:
        print (f"Cannot read CCF file {ccf_file}, it must contain extension with name \"LONG_TERM_CTI\"")
        return None
    times = tmp['TIME'][0]
    tx = t.group_by('MODE_ID')
    mmask = tx.groups.keys['MODE_ID'] == mode_id
    tenergy = tx.groups[mmask].group_by('ENERGY')
    ixxx = np.arange(len(tenergy.groups))
    if (energy_index not in ixxx):
        print (f"The index for the energy {energy_index} is not available for this mode. Only indices: ", ixxx)
        return None
    grp = tenergy.groups[energy_index]
    #print (grp)
    tout = {}
    if (mode_id != 3):
        for iccd in np.arange(1,13):
            xmask = grp['CCD_ID'] == iccd
            tcoeff = grp[xmask]['T_COEFF']
            tout[iccd] = {'mode_id': mode_id, 'ccd': iccd, 'energy': grp['ENERGY'][0],
                'times': times, 'tcoeff': tcoeff[0]}
    else:
        tcoeff = grp['T_COEFF']
        tout[0] = {'mode_id': mode_id, 'ccd': 4, 'energy': grp['ENERGY'][0],
            'times': times, 'tcoeff': tcoeff[0]}
    return tout