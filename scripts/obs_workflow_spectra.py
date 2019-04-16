#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 14 16:58:36 2018

Create a spectrum after finding the source and background regions

@author: ivaltchanov
"""

import os
import subprocess
import sys
#
import numpy as np 

#from astropy.table import Table
#from astropy.io import fits
#from astropy.stats import median_absolute_deviation as mad
#from astropy.stats import mad_std

#import matplotlib.pylab as plt

import logging
import argparse
#
# global folders
#
home = os.path.expanduser('~')
sasversion = 'xmmsas_20180504_1731'
os.environ['SAS_CCFPATH'] = '/xdata/ccf/pub'
#%%
#%%
#
#parser = argparse.ArgumentParser(description='XMMSAS workflow step 4 - spectra')
#parser.add_argument('obsid', type=str,
#                    help='The OBSID to process')
#parser.add_argument('inst', type=str,default='pn',
#                    help='Which instrument: pn, mos1 or mos2')
#parser.add_argument('odfdir', type=str, default=os.getcwd(),
#                    help='ODF folder')

#parser.add_argument('task', type=str,default='emproc',
#                    help='Which proc to run: emproc or epproc')
#
#args = parser.parse_args()
#
#obsid = args.obsid
#odfdir = args.odfdir
obsid = '0656580601'
odfdir = home + "/XMM/CAL/" + obsid

# check if ODF exists, if not fail
if (not os.path.isdir(odfdir)):
    raise FileNotFoundError
os.environ['SAS_ODF'] = odfdir
#
wdir = '{}/{}'.format(odfdir,sasversion)
if (not os.path.isdir(wdir)):
    raise FileNotFoundError
#
#obsid = '0656580601'
#obsid = '0111240101'
# check if CCF file exists in ODF dir, if not then fail
ccf_file = "{}/ccf.cif".format(odfdir)
if (not os.path.isfile(ccf_file)):
    logging.error("CCF filenot found in SAS_ODF: {}".format(odfdir))
    raise FileNotFoundError
os.environ['SAS_CCF'] = ccf_file
os.chdir(wdir)
#
logging.basicConfig(level=logging.DEBUG,
                    format='%(asctime)s %(levelname)s %(message)s',
                    filename='{}/spectra.log'.format(wdir),
                    filemode='w')
#
#%%
#
inst = 'pn'
#
evfile = "{}_evlist_clean.fits".format(inst)
if (not os.path.isfile(evfile)):
    raise FileNotFoundError
#
# read the region file it must have a name source_region_[inst].txt
# and the background region must have a name bg_region_[inst].txt
#
with open('{}/regions/source_region_{}.txt'.format(wdir,inst)) as regfile:
    src_reg = regfile.readline()
with open('{}/regions/bg_region_{}.txt'.format(wdir,inst)) as regfile:
    bg_reg = regfile.readline()
#%%
#
# now select source photons form the GTI cleaned event list
#
src_spec = "{}_src_spectrum.fits".format(inst)
command = "evselect table={} withspectrumset=yes spectrumset=".format(evfile,src_spec) +  \
" energycolumn=PI spectralbinsize=5 withspecranges=yes specchannelmin=0 specchannelmax=20479" +  \
" expression='(FLAG==0) && (PATTERN<=4) && ((X,Y) IN {})'".format(src_reg.strip())
try:
    result = subprocess.run(command, shell=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
    retcode=result.returncode
    if retcode < 0:
        print("evselect was terminated by signal", -retcode, file=sys.stderr)
        logging.warning("evselect was terminated by signal: {}".format(-retcode))
    else:
        print("evselect returned", retcode, file=sys.stderr)
        logging.info("evselect returned {}, \n {}".format(retcode,result.stdout.decode()))
except OSError as e:
    print("Execution of evselect failed:", e, file=sys.stderr)
    logging.error("Execution of evselect failed: {}".format(e))
#%%
# background selection
#

bg_spec = "{}_bkg_spectrum.fits".format(inst)
command = "evselect table={} withspectrumset=yes spectrumset={}".format(evfile,bg_spec) +  \
" energycolumn=PI spectralbinsize=5 withspecranges=yes specchannelmin=0 specchannelmax=20479" +  \
" expression='(FLAG==0) && (PATTERN<=4) && ((X,Y) IN {})'".format(bg_reg.strip())
try:
    result = subprocess.run(command, shell=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
    retcode=result.returncode
    if retcode < 0:
        print("evselect was terminated by signal", -retcode, file=sys.stderr)
        logging.warning("evselect was terminated by signal: {}".format(-retcode))
    else:
        print("evselect returned", retcode, file=sys.stderr)
        logging.info("evselect returned {}, \n {}".format(retcode,result.stdout.decode()))
except OSError as e:
    print("Execution of evselect failed:", e, file=sys.stderr)
    logging.error("Execution of evselect failed: {}".format(e))
#%%
# Calculate the area of source and background region used to make the spectral files. 
#The area is written into the header of the SPECTRUM table of the file as keyword BACKSCAL 
# (if the spectrum is created via xmmselect, backscale will run automatically).
command = "backscale spectrumset={} badpixlocation={}".format(src_spec,evfile)
try:
    result = subprocess.run(command, shell=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
    retcode=result.returncode
    if retcode < 0:
        print("backscale was terminated by signal", -retcode, file=sys.stderr)
        logging.warning("backscale was terminated by signal: {}".format(-retcode))
    else:
        print("backscale returned", retcode, file=sys.stderr)
        logging.info("backscale returned {}, \n {}".format(retcode,result.stdout.decode()))
except OSError as e:
    print("Execution of backscale failed:", e, file=sys.stderr)
    logging.error("Execution of backscale failed: {}".format(e))
#%%
# background
#
command = "backscale spectrumset={} badpixlocation={}".format(bg_spec,evfile)
try:
    result = subprocess.run(command, shell=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
    retcode=result.returncode
    if retcode < 0:
        print("backscale was terminated by signal", -retcode, file=sys.stderr)
        logging.warning("backscale was terminated by signal: {}".format(-retcode))
    else:
        print("backscale returned", retcode, file=sys.stderr)
        logging.info("backscale returned {}, \n {}".format(retcode,result.stdout.decode()))
except OSError as e:
    print("Execution of backscale failed:", e, file=sys.stderr)
    logging.error("Execution of backscale failed: {}".format(e))
#%%
# now generate the RMF (can take time)
command = "rmfgen spectrumset={} rmfset=pn.rmf".format(src_spec)
try:
    result = subprocess.run(command, shell=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
    retcode=result.returncode
    if retcode < 0:
        print("rmfgen was terminated by signal", -retcode, file=sys.stderr)
        logging.warning("rmfgen was terminated by signal: {}".format(-retcode))
    else:
        print("rmfgen returned", retcode, file=sys.stderr)
        logging.info("rmfgen returned {}, \n {}".format(retcode,result.stdout.decode()))
except OSError as e:
    print("Execution of rmfgen failed:", e, file=sys.stderr)
    logging.error("Execution of rmfgen failed: {}".format(e))
#%%
# now generate the ARF
command = "arfgen spectrumset={} arfset=pn.arf withrmfset=yes rmfset=pn.rmf".format(src_spec) + \
    " badpixlocation={} detmaptype=psf".format(evfile)
try:
    result = subprocess.run(command, shell=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
    retcode=result.returncode
    if retcode < 0:
        print("arfgen was terminated by signal", -retcode, file=sys.stderr)
        logging.warning("arfgen was terminated by signal: {}".format(-retcode))
    else:
        print("arfgen returned", retcode, file=sys.stderr)
        logging.info("arfgen returned {}, \n {}".format(retcode,result.stdout.decode()))
except OSError as e:
    print("Execution of arfgen failed:", e, file=sys.stderr)
    logging.error("Execution of arfgen failed: {}".format(e))
#
#%%
# now group/bin the spectrum
#
command = "specgroup spectrumset={} mincounts=25 oversample=3 rmfset=pn.rmf".format(src_spec) + \
    " arfset=pn.arf backgndset={} groupedset=pn_spectrum_grp25.fits".format(bg_spec)
try:
    result = subprocess.run(command, shell=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
    retcode=result.returncode
    if retcode < 0:
        print("specgroup was terminated by signal", -retcode, file=sys.stderr)
        logging.warning("specgroup was terminated by signal: {}".format(-retcode))
    else:
        print("specgroup returned", retcode, file=sys.stderr)
        logging.info("specgroup returned {}, \n {}".format(retcode,result.stdout.decode()))
except OSError as e:
    print("Execution of specgroup failed:", e, file=sys.stderr)
    logging.error("Execution of specgroup failed: {}".format(e))

print ("All done")