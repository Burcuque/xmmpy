#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 09 12:00:00 2018

Step 04 of PN SW processing:  build images in user supplied bands

Wil only process PN Large Window mode.

Will use the PN GTI and region files from the standard procesing
 
The output will be in a foder called PN_LW

Using epchain

@author: ivaltchanov
"""
import os
import subprocess
import sys
import logging

import argparse

os.environ['SAS_CCFPATH'] = '/xdata/ccf/pub'
#%%
#
# get the arguments
parser = argparse.ArgumentParser(description='XMM SAS step 04 for PN SW, making image')
parser.add_argument('obsid', type=str,
                    help='The OBSID to process')
parser.add_argument('expo_name', type=str,
                    help='The OBSID exposure to process')
parser.add_argument('lowe', type=int,
                    help='Low energy of the band in eV (integer)')
parser.add_argument('hie', type=int,
                    help='High energy of the band in eV (integer)')
#
args = parser.parse_args()
#
#%%
#
obsid = args.obsid
nexpo = args.expo_name
odfdir = os.path.join(os.getcwd(),obsid)
pi0 = args.lowe
pi1 = args.hie
#
if (not os.path.isdir(odfdir)):
    print (f"The ODF folder {odfdir} does not exist! Cannot continue")
    raise FileNotFoundError
#
ppsDir = os.path.join(odfdir,"PN_LW")
#
if (not os.path.isdir(ppsDir)):
    print (f"PPS folder {ppsDir} does not exist. This means there are no event lists produced")
    print (f"Have you run pn_sw_step02.py ?")
    raise FileNotFoundError
#
logging.basicConfig(level=logging.DEBUG,
                    format='%(asctime)s %(levelname)s %(message)s',
                    filename=os.path.join(ppsDir,'pn_lw_step04.log'),
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
#
os.chdir(ppsDir)

#%%
# Select the cleaned event lists
#
evlist = os.path.join(ppsDir,f'pn_evlist_clean_{nexpo}.fits')
if (not os.path.isfile(evlist)):
    print(f"GTI filtered event list {evlist} not available. Have you run pn_lw_step03.py ?")
    logging.error(f"GTI filtered event list {evlist} not available.")
    raise FileNotFoundError
#
#%%
#
print ("*** Generating images in band [{},{}] eV".format(pi0,pi1))
image_name = f'pn_image_{pi0}_{pi1}_{nexpo}.fits'
#
expr = f'PI in [{pi0}:{pi1}] &&  FLAG==0 && PATTERN in [0:4]'
#    
command = 'evselect table={} xcolumn=X ycolumn=Y imagebinning=binSize'.format(evlist) +  \
     ' ximagebinsize=80 yimagebinsize=80' + \
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
logging.info ("pn_lw_step04 done")
#
