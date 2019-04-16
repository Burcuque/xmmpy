#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 10 10:43:19 2018

Step 04 of the workflow for EPIC processing an XMM observation

Use the GTI filtered event lists to make images in different energy bands

@author: ivaltchanov
"""

import os
import subprocess
import sys
import glob

#import argparse


import logging

import argparse
#
# global folders
#
os.environ['SAS_CCFPATH'] = '/xdata/ccf/pub'
#
#%%
#
parser = argparse.ArgumentParser(description='XMMSAS workflow step 4 - images')
parser.add_argument('obsid', type=str,
                    help='The OBSID to process')
parser.add_argument('lowe', type=int,
                    help='Low energy of the band in eV (integer)')
parser.add_argument('hie', type=int,
                    help='High energy of the band in eV (integer)')
parser.add_argument('--expo', type=str, default="S003",
                    help='The exposure to use')
parser.add_argument('--binSize', type=int, default=80,
                    help='The image bin size (integer), default 80 ==> 4 arcsec pixel')
parser.add_argument('--odfdir', type=str, default=os.getcwd(),
                    help='Where ODF files are')
parser.add_argument('--detmask', default=True, action='store_true',\
                    help='Generate detector mask? Will generate exposure map too.')

args = parser.parse_args()
#
obsid = args.obsid
odfdir = os.path.join(os.path.abspath(args.odfdir),obsid)
expo = args.expo
pi0 = args.lowe
pi1 = args.hie
bin_size = args.binSize
make_detmask = args.detmask

#
if (not os.path.isdir(odfdir)):
    print ("The ODF folder {} does not exist! Cannot continue".format(odfdir))
    raise FileNotFoundError
#
#
# crate output folder for PPS and logging
#
# get the SAS version, will be used for the output folder 
try:
    result = subprocess.run(["sasversion -v"], shell=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
    retcode = result.returncode
    if retcode < 0:
        print("sasversion was terminated by signal {} \n {}".format(-retcode,result.stdout.decode()))
    else:
        print("sasversion returned {} \n {}".format(retcode,result.stdout.decode()))
except OSError as e:
    print("Execution of sasversion failed:", e, file=sys.stderr)
#
sasversion = result.stdout.decode().split('[')[1].split(']')[0]
#
ppsDir = os.path.join(odfdir,sasversion)
if (not os.path.isdir(ppsDir)):
    print ("PPS folder {} does not exist. This means step #01 was not run?".format(ppsDir))
    raise FileNotFoundError
#
logging.basicConfig(level=logging.DEBUG,
                    format='%(asctime)s %(levelname)s %(message)s',
                    filename=os.path.join(ppsDir,'proc_step04.log'),
                    filemode='w')
#
os.environ["SAS_ODF"] = odfdir
ccffile = os.path.join(odfdir,"ccf.cif")
if (not os.path.isfile(ccffile)):
    print ("No ccf.cif file found in ODF folder {}. This means step #01 was not run?".format(odfdir))
    raise FileNotFoundError
#
os.environ["SAS_CCF"] = ccffile
logging.info("Set SAS_ODF to {}".format(odfdir))
logging.info("Set SAS_CCF to {}".format(ccffile))
#
os.chdir(ppsDir)
#
task = "evselect"
#
print ("*** Generating PN image for exposure {} in band [{},{}] eV".format(expo,pi0,pi1))
#
evlist = f"pn_{expo}_evlist_clean.fits"
if (not os.path.isfile(evlist)):
    logging.error("Cannot find cleaned event lists")
image_name = f'pn_{expo}_image_{pi0}_{pi1}.fits'
#
expr = f'PI in [{pi0}:{pi1}] &&  FLAG==0 && PATTERN in [0:4]'
#    
command = 'evselect table={} xcolumn=X ycolumn=Y imagebinning=binSize'.format(evlist) +  \
     ' ximagebinsize={} yimagebinsize={}'.format(bin_size,bin_size) + \
     ' expression=\'{}\''.format(expr) +  \
     ' withimageset=true imageset={}'.format(image_name)
try:
    result = subprocess.run(command, shell=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
    retcode=result.returncode
    if retcode < 0:
        print("evselect was terminated by signal", -retcode, file=sys.stderr)
        logging.warning("evselect was terminated by signal: {} \n {}".format(-retcode,result.stdout.decode()))
    else:
        print("evselect returned", retcode, file=sys.stderr)
        logging.info("evselect returned {}, \n {}".format(retcode,result.stdout.decode()))
except OSError as e:
    print("Execution of evselect failed:", e, file=sys.stderr)
    logging.error("Execution of evselect failed: {}".format(e))
#
if (make_detmask):
    print ("*** Generating exposure maps and detector masks ")
    #
    atthk = glob.glob("*_AttHk.ds")
    if (not os.path.isfile(atthk[0])):
        print("Cannot find Attitude HK file")
        logging.error("Cannot find Attitude HK file")
        raise FileNotFoundError
    expfile = "pn_{}_expimage_{}_{}.fits".format(expo,pi0,pi1)
    detfile = "pn_{}_detmask_{}_{}.fits".format(expo,pi0,pi1)
    command = 'eexpmap imageset={} attitudeset={}'.format(image_name,atthk[0]) + \
    ' eventset={} expimageset={} pimin={} pimax={}'.format(evlist,expfile,pi0,pi1)
    try:
        result = subprocess.run(command, shell=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
        retcode=result.returncode
        if retcode < 0:
            print("eexpmap was terminated by signal", -retcode, file=sys.stderr)
            logging.warning("eexpmap was terminated by signal: {} \n {}".format(-retcode,result.stdout.decode()))
        else:
            print("eexpmap returned", retcode, file=sys.stderr)
            logging.info("eexpmap returned {}, \n {}".format(retcode,result.stdout.decode()))
    except OSError as e:
        print("Execution of eexpmap failed:", e, file=sys.stderr)
        logging.error("Execution of eexpmap failed: {}".format(e))
    #
    # now the detector map
    #
    command = 'emask expimageset={} detmaskset={}'.format(expfile,detfile) + \
    ' threshold1=0.3 threshold2=0.5'
    try:
        result = subprocess.run(command, shell=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
        retcode=result.returncode
        if retcode < 0:
            print("emask was terminated by signal", -retcode, file=sys.stderr)
            logging.warning("emask was terminated by signal: {} \n {}".format(-retcode,result.stdout.decode()))
        else:
            print("emask returned", retcode, file=sys.stderr)
            logging.info("emask returned {}, \n {}".format(retcode,result.stdout.decode()))
    except OSError as e:
        print("Execution of emask failed:", e, file=sys.stderr)
        logging.error("Execution of emask failed: {}".format(e))
        

logging.info ("All done")
#
