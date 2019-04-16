#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 09 12:00:00 2018

Full workflow for EPIC processing an XMM observation with no CTI correction applied.

Wil only process PN SmallWindow mode in CALCLOSED
 
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
parser = argparse.ArgumentParser(description='XMM SAS workflow for PN SW in CALCLOSED')
parser.add_argument('obsid', type=str,
                    help='The OBSID to process')
#
args = parser.parse_args()
#
#%%
#
obsid = args.obsid
odfdir = f'/xdata/xcaldata/XMM/IVAN/PN_calclosed/{obsid}'
#
#skipDone = args.skip
#obsid = "0771000201"
#odfdir = f'/xdata/xcaldata/XMM/IVAN/PN_SW/ngc5548/{obsid}'
#skipDone = True
#
if (not os.path.isdir(odfdir)):
    print ("The ODF folder {} does not exist! Cannot continue".format(odfdir))
    raise FileNotFoundError

ppsDir = os.path.join(odfdir,'no_cti')
#
if (not os.path.isdir(ppsDir)):
    print ("PPS folder {} does not exist. Will create it".format(ppsDir))
    os.mkdir(ppsDir)
#
os.environ["SAS_ODF"] = odfdir
ccffile = os.path.join(odfdir,"ccf.cif")
if (not os.path.isfile(ccffile)):
    print ("No ccf.cif file found in ODF folder {}. Run ccfbuild".format(odfdir))
    raise FileNotFoundError
#
#%%
#
task = 'epchain'
#
logging.basicConfig(level=logging.DEBUG,
                    format='%(asctime)s %(levelname)s %(message)s',
                    filename=os.path.join(ppsDir,'pn_proc_nocti.log'.format(task)),
                    filemode='w')
#
os.environ["SAS_CCF"] = ccffile
logging.info("Set SAS_ODF to {}".format(odfdir))
logging.info("Set SAS_CCF to {}".format(ccffile))
#
#
#%%
# Have to run ODFINGEST as I moved the files
#
print ("*** ODFINGEST")
os.chdir(odfdir)
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
#
#%%
os.chdir(ppsDir)
#
xtask = task + f" withctilongterm=N odf={odfdir} withphagaincolumn=yes propagatecolumns=all "
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
evfind = glob.glob(f"{ppsDir}/*PIEVLI*")
if (len(evfind) == 1):
    evlist = evfind[0]
    haveEvlist = True
else:
    logging.error("Event list not produced")
    raise FileNotFoundError

#%%
# 
# For CALCLOSED no need to filter for GTI, directly extract a spectrum
#
spec_bin_size = 5
specFile = 'spectrum_bin{}.fits'.format(spec_bin_size)

evfile = glob.glob("%s/*PIEVLI*"%ppsDir)
if (len(evfile) == 0):
    print ("No event list file was created. Check the processing, cannot continue!")
    raise FileNotFoundError
#
evselect_command = "evselect " + \
        "table=%s"%evfile[0] + " energycolumn='PI' withspectrumset=yes" + \
        " expression='(((DETX,DETY) in BOX(517,1430,2500,2250,0)) && (FLAG==0) " + \
        " && (PATTERN==0) && (PAT_SEQ==0))' " + \
        " xcolumn=DETX ycolumn=DETY" +  \
        " withspecranges=yes specchannelmin=0 specchannelmax=20479" + \
        " spectrumset={} spectralbinsize={}".format(specFile,spec_bin_size)
try:
    result = subprocess.run(evselect_command, shell=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
    retcode = result.returncode
    if retcode < 0:
        print("evselect was terminated by signal", -retcode, file=sys.stderr)
        logging.warning("evselect was terminated by signal: {}".format(-retcode))
    else:
        print("evselect returned", retcode, file=sys.stderr)
        logging.info("evselect returned {}, \n {}".format(retcode,result.stdout.decode()))
except OSError as e:
    print("Execution of evselect failed:", e, file=sys.stderr)
    logging.error("Execution of evselect failed: {}".format(e))
#