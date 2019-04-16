#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 09 12:00:00 2018

Step 02 of the workflow for EPIC processing an XMM observation

emproc and epproc tasks

version 1

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
parser = argparse.ArgumentParser(description='XMM SAS workflow step 02 em/ep proc')
parser.add_argument('obsid', type=str,
                    help='The OBSID to process')
parser.add_argument('inst', type=str, default='pn',
                    help='Which instrument to process')
parser.add_argument('--odfdir', type=str, default=os.getcwd(),
                    help='Where the ODF files are')
#
args = parser.parse_args()
#
#%%
#
obsid = args.obsid
odfdir = os.path.join(os.path.abspath(args.odfdir),obsid)
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
if ('pn' in args.inst):
    task = 'epproc'
else:
    task = 'emproc'
#
logging.basicConfig(level=logging.DEBUG,
                    format='%(asctime)s %(levelname)s %(message)s',
                    filename=os.path.join(ppsDir,'proc_{}_step02.log'.format(task)),
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
#
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
