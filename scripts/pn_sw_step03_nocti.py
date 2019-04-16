#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 09 12:00:00 2018

Step 03 of PN SW processing:  GTI selection and extracting filtered events 

Wil only process PN SmallWindow mode.

Will use the PN GTI and region files from the standard procesing
 
The output will be in a foder called PN_SW

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
parser = argparse.ArgumentParser(description='XMM SAS step 03 for PN SW, GTI filtering')
parser.add_argument('obsid', type=str,
                    help='The OBSID to process')
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
ppsDir = os.path.join(odfdir,"no_cti")
#
if (not os.path.isdir(ppsDir)):
    print (f"PPS folder {ppsDir} does not exist. This means there are no event lists produced")
    print (f"Have you run pn_sw_step02.py ?")
    raise FileNotFoundError
#
logging.basicConfig(level=logging.DEBUG,
                    format='%(asctime)s %(levelname)s %(message)s',
                    filename=os.path.join(ppsDir,'pn_sw_step03.log'),
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
    lt_max = 0.0
    for xev in evfind:
        hdu = fits.open(xev)
        lt = hdu[1].header['ONTIME']
        if (lt >= lt_max):
            evlist = xev
            lt_max = lt
        #
    haveEvlist = True
else:
    logging.error("Event list not produced")
    raise FileNotFoundError
#
#%%
#  ========== GTI filtering not needed will use the one from the standard epchain ============
#
#%%
# now filter the event lists with the GTI already available 
#
print ("Filtering the calibrated event lists with the GTI")
gti_file = "../PN_SW/gti_pn.fits"
if (not os.path.isfile(gti_file)):
    logging.error ("GTI file {} not found".format(gti_file))
    raise FileNotFoundError
#
# check the duration of the GTI
#
hdu = fits.open(gti_file)
gti_time = hdu[1].header['ONTIME']
t_ratio = gti_time/lt_max
print ("GTI filtered time is {}, the original is {}, fraction {}".format(gti_time,lt_max,t_ratio))
if (t_ratio <= 0.5):
    print ("Warning! {} time fraction is discarded for high background".format(t_ratio))
    logging.warning("Time fraction is discarded for high background: {}".format(t_ratio))
#
filt_evlist = f"{ppsDir}/pn_evlist_clean_nocti.fits"

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
