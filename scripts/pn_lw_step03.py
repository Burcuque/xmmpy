#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 09 12:00:00 2018

Step 03 of PN LW processing:  GTI selection and extracting filtered events 

Wil only process PN LargeWindow mode.

Will use the PN GTI and region files from the standard procesing
 
The output will be in a foder called PN_LW

Using epchain

@author: ivaltchanov
"""
import os
import subprocess
import sys
import logging
import glob

import numpy as np

from astropy.io import fits
from astropy.stats import mad_std

import matplotlib as mpl
mpl.use('Agg')

import matplotlib.pylab as plt

import argparse

os.environ['SAS_CCFPATH'] = '/xdata/ccf/pub'
#%%
#
# get the arguments
parser = argparse.ArgumentParser(description='XMM SAS step 03 for PN LW, GTI filtering')
parser.add_argument('obsid', type=str,
                    help='The OBSID to process')
#parser.add_argument('expo', type=str, default="all",
#                    help='The OBSID Exposure to process, default all')
#
args = parser.parse_args()
#
#%%
#
obsid = args.obsid
odfdir = os.path.join(os.getcwd(),obsid)
#
if (not os.path.isdir(odfdir)):
    print (f"The ODF folder {odfdir} does not exist! Cannot continue")
    raise FileNotFoundError
#
ppsDir = os.path.join(odfdir,"PN_LW")
#
if (not os.path.isdir(ppsDir)):
    print (f"PPS folder {ppsDir} does not exist. This means there are no event lists produced")
    print (f"Have you run pn_lw_step02.py ?")
    raise FileNotFoundError
#
logging.basicConfig(level=logging.DEBUG,
                    format='%(asctime)s %(levelname)s %(message)s',
                    filename=os.path.join(ppsDir,'pn_lw_step03.log'),
                    filemode='w')
#
os.environ["SAS_ODF"] = odfdir
ccffile = os.path.join(odfdir,"ccf.cif")
if (not os.path.isfile(ccffile)):
    print (f"No ccf.cif file found in ODF folder {odfdir}")
    print ("     ===> This means pn_sw_step01.py was not run?")
    logging.error(f"No CCF file ccf.cif found in {odfdir}")
    raise FileNotFoundError
#
os.environ["SAS_CCF"] = ccffile
logging.info("Set SAS_ODF to {}".format(odfdir))
logging.info("Set SAS_CCF to {}".format(ccffile))

os.chdir(ppsDir)

#%%
# Select the event lists
#
haveEvlist = False
evfind = glob.glob(f"{ppsDir}/*PIEVLI*")
if (len(evfind) >= 1):
    logging.info("Found {} event lists".format(len(evfind)))
    haveEvlist = True
else:
    logging.error("Event list not produced")
    raise FileNotFoundError
#
# now loop over all exposures 

for evlist in evfind:
    nexpo = os.path.basename(evlist)[13:17]
    #  ========== GTI filtering ============
    #
    #
    print (f"Generating rate curve for expo {nexpo}")
    rate = f'rate_pn_{nexpo}.fits'
    expr = ' #XMMEA_EP && (PI > 10000 && PI < 12000) && (PATTERN == 0)'
    #    
    command = 'evselect table={} withrateset=Y rateset={}'.format(evlist,rate) + \
    ' maketimecolumn=Y timebinsize=100 makeratecolumn=Y' + \
    ' expression=\'{}\''.format(expr)
    try:
        result = subprocess.run(command, shell=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
        retcode=result.returncode
        if retcode < 0:
            logging.warning("evselect was terminated by signal {} \n {}".format(-retcode,result.stdout.decode()))
            print("evselect was terminated by signal {} \n {}".format(-retcode,result.stdout.decode()))
        else:
            print("evselect returned {}, \n {}".format(retcode,result.stdout.decode()))
            logging.info("evselect returned {}".format(retcode,result.stdout.decode()))
    except OSError as e:
        print("Execution of evselect failed: {}".format(e))
        logging.error("Execution of evselect failed: {}".format(e))
    #
    # now selet the GTI using the lightcurve
    # the count-rate limit for filtering = median + 3*MAD
    #
    tpn = fits.open(f'rate_pn_{nexpo}.fits')
    time_min = np.min(tpn[1].data['TIME'])
    expo_time = tpn['RATE'].header['EXPOSURE'] 
    #
    maxrate = np.max(tpn[1].data['RATE'])
    medrate = np.median(tpn[1].data['RATE'])
    xmad = mad_std(tpn[1].data['RATE'])
    ulimit = medrate + 3*xmad
    rate_lim = 0.4
    #
    use_limit = max(ulimit,rate_lim)
    logging.info("PN: actual rate limit to use for GTI filtering: {:.3f}".format(use_limit))
    #
    #%%
    print ("Running tabgtigen")
    gti_command = f'tabgtigen table=rate_pn_{nexpo}.fits' + \
    ' expression=\'RATE<={}\' gtiset=gti_pn_{}.fits'.format(use_limit,nexpo)
    try:
        result = subprocess.run(gti_command, shell=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
        retcode=result.returncode
        if retcode < 0:
            print("tabgtigen was terminated by signal", -retcode, file=sys.stderr)
            logging.warning("tabgtigen was terminated by signal: {}".format(-retcode))
        else:
            print("tabgtigen returned", retcode, file=sys.stderr)
            logging.info("tabgtigen returned {}, \n {}".format(retcode,result.stdout.decode()))
    except OSError as e:
        print("Execution of tabgtigen failed:", e, file=sys.stderr)
        logging.error("Execution of tabgtigen failed: {}".format(e))
    #
    #%%
    # now filter the event lists with the merged GTI
    print (f"Filtering the calibrated event lists with the GTI for expo {nexpo}")
    gti_file = f"gti_pn_{nexpo}.fits"
    if (not os.path.isfile(gti_file)):
        logging.error ("GTI file {} not found".format(gti_file))
        raise FileNotFoundError
    #
    # check the duration of the GTI
    #
    hdu = fits.open(gti_file)
    gti_time = hdu[1].header['ONTIME']
    t_ratio = gti_time/expo_time
    print ("GTI filtered time is {}, the original is {}, fraction {}".format(gti_time,expo_time,t_ratio))
    if (t_ratio <= 0.5):
        print ("Warning! {} time fraction is discarded for high background".format(t_ratio))
        logging.warning("Time fraction is discarded for high background: {}".format(t_ratio))
    #
    filt_evlist = f"{ppsDir}/pn_evlist_clean_{nexpo}.fits"
    
    expr = "#XMMEA_EP && gti({},TIME) && (PI>150)".format(gti_file)
    #
    ev_command = 'evselect table={} withfilteredset=Y filteredset={}'.format(os.path.basename(evlist),os.path.basename(filt_evlist)) +  \
       ' destruct=Y keepfilteroutput=T expression=\'{}\''.format(expr)
    #
    try:
        result = subprocess.run(ev_command, shell=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
        retcode=result.returncode
        if retcode < 0:
            print("evselect was terminated by signal", -retcode, file=sys.stderr)
            logging.warning("evselect was terminated by signal {}".format(-retcode))
        else:
            print("evselect returned", retcode, file=sys.stderr)
            logging.info("evselect returned {}, \n {}".format(retcode,result.stdout.decode()))
    except OSError as e:
        print("Execution of evselect failed:", e, file=sys.stderr)
        logging.error("Execution of evselect failed: {}".format(e))
        #
    #
    #%%
    #
    # plot the merged GTI on lightcurve
    
    #
    f, ax = plt.subplots()
    start = hdu[1].data["START"]
    end = hdu[1].data["STOP"]
    ngti = len(start)
    #
    reltime = tpn[1].data['TIME'] - time_min
    ax.plot(reltime,tpn[1].data['RATE'],label="PN")
    rate_lim = 0.4
    #
    ax.axhline(y=rate_lim,color='magenta',ls='solid',label='SAS')
    ax.axhline(y=use_limit,color='k',ls='dashed',label='Adopted'.format(use_limit))
    ax.set_ylim((0.0,maxrate))
    ax.grid(True)
    ax.legend(loc=2)
    ax.set_ylabel("Count rate (cts/s)")
    #
    # now the GTI bands
    #
    for jj in np.arange(ngti):
        xx = (start[jj]-time_min,end[jj]-time_min)
        yy1 = (0.0,0.0)
        yy2 = (maxrate,maxrate)
        ax.fill_between(xx,yy1,yy2,facecolor='yellow',alpha=0.3)
    f.subplots_adjust(hspace=0)
    plt.setp([a.get_xticklabels() for a in f.axes[:-1]], visible=False)
    plt.xlabel("Relative time (s)")
    ax.set_title("GTI for PN {} {}, ONTIME={} s".format(obsid,nexpo, gti_time))
    #
    png_file = "{}/{}_rates_gti_{}.png".format(ppsDir,obsid,nexpo)
    plt.savefig(png_file,dpi=100)
    #plt.show()
    logging.info("Saved GTI limits in figure: {}".format(png_file))
    plt.close()
