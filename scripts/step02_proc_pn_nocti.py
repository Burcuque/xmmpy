#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 09 12:00:00 2018

Step 02 of the workflow for EPIC processing an XMM observation

This is a version of step02_proc.py, only for PN, only for SW and
will not apply the long-term CTI on the results.

The output will be in a foder called no_cti

no long-term CTI parameter is only available in epchain script.

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
parser = argparse.ArgumentParser(description='XMM SAS workflow step 02 em/ep proc')
parser.add_argument('obsid', type=str,
                    help='The OBSID to process')
parser.add_argument('--odfdir', type=str, default=os.getcwd(),
                    help='Where the ODF files are')
#
args = parser.parse_args()
#
#%%
#
odfdir = os.path.abspath(args.odfdir)
obsid = args.obsid
#
if (not os.path.isdir(odfdir)):
    print ("The ODF folder {} does not exist! Cannot continue".format(odfdir))
    raise FileNotFoundError
#
# crate output folder for PPS and logging
#
# get the SAS version, will be used for the output folder 
#try:
#    result = subprocess.run(["sasversion -v"], shell=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
#    retcode = result.returncode
#    if retcode < 0:
#        print("sasversion was terminated by signal {} \n {}".format(-retcode,result.stdout.decode()))
#    else:
#        print("sasversion returned {} \n {}".format(retcode,result.stdout.decode()))
#except OSError as e:
#    print("Execution of sasversion failed:", e, file=sys.stderr)
#
#sasversion = result.stdout.decode().split('[')[1].split(']')[0]
#
#ppsDir = os.path.join(odfdir,sasversion)
#
ppsDir = os.path.join(odfdir,"no_cti")
#
if (not os.path.isdir(ppsDir)):
    print ("PPS folder {} does not exist. Will create it".format(ppsDir))
    os.mkdir(ppsDir)
#
task = 'epchain'
#
logging.basicConfig(level=logging.DEBUG,
                    format='%(asctime)s %(levelname)s %(message)s',
                    filename=os.path.join(ppsDir,'proc_{}_step02_nocti.log'.format(task)),
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
# Have to run ODFINGEST as I moved the files
#
print ("*** ODFINGEST")
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
# parse the ODF SUmmary SAS file to only select observations with SmallWindow mode for PN
#
sumFile = glob.glob("{}/*.SAS".format(odfdir))
if (len(sumFile) != 1):
    raise FileNotFoundError
#
pnSW = False
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
os.chdir(ppsDir)
#
#
task += " withctilongterm=N odf={}".format(odfdir)
print ("*** RUNNING %s"%task)
try:
    result = subprocess.run(task, shell=True,stderr=subprocess.STDOUT,stdout=subprocess.PIPE)
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
