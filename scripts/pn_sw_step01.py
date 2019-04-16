#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 09 12:00:00 2018

Step #1 of processing PN SmallWindow observations:
 1. odfingest
 2. cifbuild

Will only process PN SmallWindow mode.
 
The output will be in the ODF_DIR folder

Using epchain

@author: ivaltchanov
"""
import os
import subprocess
import logging
import glob


import argparse

# global place for XMM current calibration files (CCF)
#
os.environ['SAS_CCFPATH'] = '/xdata/ccf/pub'
#%%
#
# get the arguments
parser = argparse.ArgumentParser(description='XMM SAS PN SW step01')
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
    print ("The ODF folder {} does not exist! Cannot continue".format(odfdir))
    raise FileNotFoundError
#
logging.basicConfig(level=logging.DEBUG,
                    format='%(asctime)s %(levelname)s %(message)s',
                    filename=os.path.join(odfdir,'pn_sw_step01.log'),
                    filemode='w')
#
os.environ["SAS_ODF"] = odfdir
#
os.chdir(odfdir)
#%%
#
# ======== CIFBUILD ============
#
ccffile = os.path.join(odfdir,"ccf.cif")
if (not os.path.isfile(ccffile)):
    try:
        result = subprocess.run("cifbuild", shell=True,stderr=subprocess.STDOUT,stdout=subprocess.PIPE)
        retcode=result.returncode
        if retcode < 0:
            print("cifbuild was terminated by signal {} \n {}".format(-retcode))
            logging.warning("cifbuild was terminated by signal {}".format(-retcode))
        else:
            print("cifbuild returned {}".format(retcode))
            logging.info("cifbuild returned {} \n {}".format(retcode,result.stdout.decode()))
    except OSError as e:
        print("Execution of cifbuild failed: {}".format(e))
        logging.error("Execution of cifbuild failed: {}".format(e))
#
os.environ["SAS_CCF"] = ccffile
logging.info("Set SAS_ODF to {}".format(odfdir))
logging.info("Set SAS_CCF to {}".format(ccffile))

#%%
#
# ======== ODFINGEST
try:
    result = subprocess.run("odfingest withodfdir=yes odfdir=\"{}\"".format(odfdir), shell=True,stderr=subprocess.STDOUT,stdout=subprocess.PIPE)
    retcode = result.returncode
    if retcode < 0:
        print("odfingest was terminated by signal {}".format(-retcode))
        logging.warning("odfingest was terminated by signal {}".format(-retcode))
    else:
        print("odfingest returned {} \n {}".format(retcode,result.stdout.decode()))
        logging.info("odfingest returned {} \n {}".format(retcode,result.stdout.decode()))
except OSError as e:
    print("Execution of odfingest failed: {}".format(e))
    logging.error("Execution of odfingest failed: {}".format(e))
#
#%%
# parse the ODF Summary SAS file to only select observations with SmallWindow mode for PN
#
pnSW = False
sumFile = glob.glob("{}/*.SAS".format(odfdir))

if (len(sumFile) != 1):
    print (f"No ODF summary file *.SAS found in ODF folder {odfdir}")
    raise FileNotFoundError
#
with open(sumFile[0],'r') as sas:
    lines = sas.readlines()
for qline in lines:
    if ("MODE = PrimeSmallWindow" in qline):
        pnSW = True
        break
#
if (not pnSW):
    print (f"Not processing {obsid} as it is not with PN Small Window Mode")
    raise Exception
else:
    print (f"Processing {obsid} as it is with PN Small Window Mode")
#
#
print ("PN SW Step 01 done")