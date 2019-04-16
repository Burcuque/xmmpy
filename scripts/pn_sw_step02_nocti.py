#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 09 12:00:00 2018

Step 02 of PN SW processing:  running epchain 

Will only process PN SmallWindow mode.
 
The output will be in a foder called PN_SW

Using epchain

@author: ivaltchanov
"""
import os
import subprocess
import sys
import logging
import glob

import argparse

os.environ['SAS_CCFPATH'] = '/xdata/ccf/pub'
#%%
#
# get the arguments
parser = argparse.ArgumentParser(description='XMM SAS step02 for PN SW, running epchain')
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
ppsDir = os.path.join(odfdir,"no_cti")
#
if (not os.path.isdir(ppsDir)):
    print ("PPS folder {} does not exist. Will create it".format(ppsDir))
    os.mkdir(ppsDir)
#
logging.basicConfig(level=logging.DEBUG,
                    format='%(asctime)s %(levelname)s %(message)s',
                    filename=os.path.join(ppsDir,'pn_sw_step2.log'),
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
#%%
# parse the ODF Summary SAS file to only select observations with SmallWindow mode for PN
#
pnSW = False
sumFile = glob.glob("{}/*.SAS".format(odfdir))

if (len(sumFile) != 1):
    print (f"No ODF summary file *.SAS found in ODF folder {odfdir}")
    print ("     ===> This means pn_sw_step01.py was not run?")
    logging.error(f"No ODF summary file *.SAS found in {odfdir}")
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
    logging.error(f"Observation {obsid} not in Small Window mode")
    raise Exception
else:
    print (f"Processing {obsid} as it is with PN Small Window Mode")
    logging.info(f"Observation {obsid} is in Small Window mode")
#
#
#%%
#
# ======== EPCHAIN =============
#
os.chdir(ppsDir)
#
task = "epchain"
xtask = f"{task} withctilongterm=N odf={odfdir} odfaccess=all"
print ("*** RUNNING %s"%xtask)
try:
    result = subprocess.run(xtask, shell=True,stderr=subprocess.STDOUT,stdout=subprocess.PIPE)
    retcode = result.returncode
    if retcode < 0:
        print("%s was terminated by signal"%task, -retcode, file=sys.stderr)
        logging.warning("{} was terminated by signal {} \n {}".format(task, -retcode,result.stdout.decode()))
    else:
        print("%s returned"%task, retcode, file=sys.stderr)
        logging.info("{} returned {} \n {}".format(task, retcode, result.stdout.decode()))
except OSError as e:
    print("Execution of %s failed:"%task, e, file=sys.stderr)
    logging.error("Execution of {} failed: {}".format(task, e))
#
print ("*** all done for {} and {}".format(obsid,task))
logging.info ("*** all done for {} and {}".format(obsid,task))
#
