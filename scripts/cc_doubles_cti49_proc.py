#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 09 12:00:00 2018

Full workflow for EPIC processing an XMM observation with no CTI correction applied.

Wil only process PN SmallWindow mode in CALCLOSED
 
Processing only double events (added 25 Feb 2019)

The output will be in a foder called cti49

Will use the new EPN_CTI_0049.CCF

@author: ivaltchanov
"""
import os
import subprocess
import sys
import logging
import glob

import argparse

from astropy.io import fits

os.environ['SAS_CCFPATH'] = '/xdata/ccf/pub'
#%%
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
odfdir = f'/xdata/xcaldata/XMM/IVAN/PN_SW/CalClosed/{obsid}'
#
#skipDone = args.skip
#obsid = "0771000201"
#odfdir = f'/xdata/xcaldata/XMM/IVAN/PN_SW/ngc5548/{obsid}'
#skipDone = True
#
if (not os.path.isdir(odfdir)):
    print ("The ODF folder {} does not exist! Cannot continue".format(odfdir))
    raise FileNotFoundError
#
ppsDir = os.path.join(odfdir,'cti49')
#
if (not os.path.isdir(ppsDir)):
    print ("PPS folder {} does not exist. Will create it".format(ppsDir))
    os.mkdir(ppsDir)
#
os.environ["SAS_ODF"] = odfdir
#
# Now rerun cifbuild with the new CCf file
#
logging.basicConfig(level=logging.DEBUG,
                    format='%(asctime)s %(levelname)s %(message)s',
                    filename=os.path.join(ppsDir,'pn_proc_doubles_cti49.log'),
                    filemode='w')

os.environ["SAS_CCFPATH"] = "/ccf/valid:/xdata/xcaldata/XMM/IVAN/ccfdev"
#
os.chdir(odfdir)
#
cif_file = "ccf_test.cif"

cifbuild = f"cifbuild calindexset={cif_file} withccfpath=yes ccfpath=$SAS_CCFPATH fullpath=yes -w 1"
status = run_command(cifbuild)
if (status != 0):
    raise Exception
#
#
ccffile = os.path.join(odfdir,cif_file)
if (not os.path.isfile(ccffile)):
    print ("No ccf.cif file found in ODF folder {}. Run ccfbuild".format(odfdir))
    raise FileNotFoundError
#
#%%
#
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

odf_comm = "odfingest withodfdir=yes odfdir=\"{}\"".format(odfdir)
status = run_command(odf_comm)
if (status != 0):
    raise Exception
#
#%%
os.chdir(ppsDir)
#

xtask = f"epchain odf={odfdir} withphagaincolumn=yes propagatecolumns=all backgroundtres=N"
print ("*** RUNNING %s"%xtask)
status = run_command(xtask)
if (status != 0):
    raise Exception
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
specFile = 'spectrum_doubles_bin{}.fits'.format(spec_bin_size)

evfile = glob.glob("%s/*PIEVLI*"%ppsDir)
if (len(evfile) == 0):
    print ("No event list file was created. Check the processing, cannot continue!")
    raise FileNotFoundError
#
expr1 = "((DETX,DETY) in BOX(517,1430,2500,2250,0))"
expr2 = "(FLAG == 0) && (PATTERN >= 1 && PATTERN <= 4)"
ev_comm = "evselect " + \
        f"table={evfile[0]} energycolumn='PI' withspectrumset=yes" + \
        f" expression='{expr1} && {expr2}'" + \
        " xcolumn=DETX ycolumn=DETY" +  \
        " withspecranges=yes specchannelmin=0 specchannelmax=20479" + \
        f" spectrumset={specFile} spectralbinsize={spec_bin_size}"

status = run_command(ev_comm)
if (status != 0):
    raise Exception

#%%
# Now preapring for the XSPEC
#
outfile = 'pn_spectrum_doubles_grp0.fits'
spec = 'spectrum_doubles_bin5.fits'
hdu = fits.open(spec)
hdu[1].header['BACKFILE'] = 'NONE'
hdu.writeto(spec,overwrite=True)
hdu.close()
#%%
#
# Wil use already available canned RMF and ARF
#
rmfset = '/home/ivaltchanov/IVAN/PN_SW/epn_e3_sw20_sY9_v17.0.rmf'
arfset = '/home/ivaltchanov/IVAN/PN_SW/cc_ancillary.ds'
#
command = "specgroup spectrumset={} addfilenames=yes rmfset={}".format(spec,rmfset) + \
    " arfset={} groupedset={}".format(arfset,outfile)
print ("*** Link bkg, ARF and RMF filenames to the {obsid} spectrum")

status = run_command(command)
if (status != 0):
    raise Exception

print ("All done")