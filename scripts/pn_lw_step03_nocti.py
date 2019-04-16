#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 09 12:00:00 2018

Step 02 of PN SW processing:  running epchain with the long-term CTI turned off.

Will only process PN LargeWindow mode.
 
The output will be in a foder called nocti

Using epchain

@author: ivaltchanov
"""
import os
import subprocess
import sys
import logging
import glob

import argparse

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

os.environ['SAS_CCFPATH'] = '/xdata/ccf/pub'
#%%
#
# get the arguments
parser = argparse.ArgumentParser(description='XMM SAS step03 with no CTI for PN LW, running epchain')
parser.add_argument('obsid', type=str,
                    help='The OBSID to process')
parser.add_argument('--wdir', type=str,
                    help='The working directory (root folder)')
#
args = parser.parse_args()
#
#%%
#
obsid = args.obsid
wdir = os.path.abspath(args.wdir)
odfdir = os.path.join(wdir,obsid)
#
if (not os.path.isdir(odfdir)):
    print (f"The ODF folder {odfdir} does not exist! Cannot continue.")
    raise FileNotFoundError
#
ppsDir = os.path.join(odfdir,"nocti")
#
if (not os.path.isdir(ppsDir)):
    print (f"PPS folder {ppsDir} does not exist. Will create it")
    os.mkdir(ppsDir)
#
logging.basicConfig(level=logging.DEBUG,
                    format='%(asctime)s %(levelname)s %(message)s',
                    filename=os.path.join(ppsDir,'pn_lw_step03_nocti.log'),
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
logging.info(f"Set SAS_ODF to {odfdir}")
logging.info(f"Set SAS_CCF to {odfdir}")
#%%
# parse the ODF Summary SAS file to only select observations with SmallWindow mode for PN
#
pnLW = False
sumFile = glob.glob(f"{odfdir}/*.SAS")

if (len(sumFile) != 1):
    print (f"No ODF summary file *.SAS found in ODF folder {odfdir}")
    print ("     ===> This means pn_sw_step01.py was not run?")
    logging.error(f"No ODF summary file *.SAS found in {odfdir}")
    raise FileNotFoundError
#
with open(sumFile[0],'r') as sas:
    lines = sas.readlines()
for qline in lines:
    if ("MODE = PrimeLargeWindow" in qline):
        pnLW = True
        break
#
if (not pnLW):
    print (f"Not processing {obsid} as it is not with PN Large Window Mode")
    logging.error(f"Observation {obsid} not in Large Window mode")
    raise Exception
else:
    print (f"Processing {obsid} as it is with PN Large Window Mode")
    logging.info(f"Observation {obsid} is in Large Window mode")
#
#
#%%
#
# ======== EPCHAIN =============
#
os.chdir(ppsDir)
#
task = "epchain"
xtask = f"{task} odf={odfdir} odfaccess=all backgroundtres=N "
print (f"*** RUNNING {xtask}")
status = run_command(xtask)
if (status != 0):
    print("Execution of {xtask} failed.")
    logging.error(f"Execution of {xtask} failed.")
    raise Exception
#
print (f"*** step02 done for {obsid}".format(obsid,task))
logging.info (f"*** step02 done for {obsid}".format(obsid,task))
#
