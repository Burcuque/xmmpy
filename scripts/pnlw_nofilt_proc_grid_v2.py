#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 09 12:00:00 2018

PN LW processing:  running epchain with the new testing CCF.

If <subfolder> contains "nocti" then it will not apply the long-term CTI.

Will only process PN Large Window mode.
 
The output will be in a foder called <subfolder>

If there are more than one exposures will pick up the one wit hthe largest ONTIME

Version adapted to use the local disk space on each grid node. This relies on existing 
$TMPDIR environment variable. If it does not exist then it will fall back to the normal 
version, which uses /xdata/xcaldata disk.

If $TMPDIR is available, then it will copy the content of the ODF_DIR and will create subfolder
on the node and at the end of the script it will copy it back to /xdata/xcaldata

Modification history:
    v1: will use the ODF from /xdata/xcaldata/XMM/IVAN but will create the products in $TMPDIR
    v2 (this): will copy the ODF to $TMPDIR and will create the products there

Using epchain

@author: ivaltchanov
"""
import os
import shutil
import subprocess
import sys
import logging
import glob
from astropy.io import fits


import argparse
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
#
#%%
#
# get the arguments
parser = argparse.ArgumentParser(description='XMM SAS full processing chain for PN LW')
parser.add_argument('obsid', type=str,
                    help='The OBSID to process')
parser.add_argument('subfolder', type=str,default="cti49x",
                    help='The subfolder where to save the products')
parser.add_argument('--wdir', type=str,default=os.getcwd(),
                    help='The working directory (root folder)')
parser.add_argument('-t','--test_ccf', action='store_true',
                    help='If used will use a testing CCF in IVAN/ccfdev')
parser.add_argument('-p','--patterns', action='store_true',
                    help='If set will create spectra per pattern from 0 to 4, otherwise only 0.')
parser.add_argument('-c','--clean',  action='store_true',
                    help='If set will clean the subfolder before processing')
#
args = parser.parse_args()
#
#%%
#
home = os.path.expanduser('~')
copy_to_node = False
if ('TMPDIR' in os.environ.keys()):
    copy_to_node = True
    tmpdir = os.environ["TMPDIR"]
#
#%%
obsid = args.obsid
wdir = os.path.abspath(args.wdir)
odfdir = os.path.join(wdir,obsid)
subfolder = args.subfolder
test_ccf = args.test_ccf
patts = args.patterns
clean = args.clean
#
logging.basicConfig(level=logging.DEBUG,
                    format='%(asctime)s %(levelname)s %(message)s',
                    filename=os.path.join(wdir,f'pnlw_nofilt_proc_{obsid}_{subfolder}.log'),
                    filemode='w')

if (not os.path.isdir(odfdir)):
    print (f"The ODF folder {odfdir} does not exist! Cannot continue.")
    logging.error (f"The ODF folder {odfdir} does not exist! Cannot continue.")
    raise FileNotFoundError
#
# now copy the ODF folder to tge grid node
#
if (copy_to_node):
    print ("Copying the ODF products to the grid node folder")
    logging.info ("Copying the ODF products to the grid node")
    os.chdir(wdir)
    destdir = os.path.join(tmpdir,obsid)
    shutil.copytree(obsid,destdir)
    #
    odfdir = os.path.join(tmpdir,obsid)
#
ppsDir = os.path.join(odfdir,subfolder)
if (copy_to_node):
    print ("Will save PPS products to the grid node")
    logging.info ("Will save PPS products to the grid node")
#
# Now clean the folder if asked to
#
if (os.path.isdir(ppsDir) and clean):
    print (f"Found an already existing folder {ppsDir}, removing the files as clean is {clean}.")
    logging.info (f"Found an already existing folder {ppsDir}, removing the files as clean is {clean}.")
    shutil.rmtree(ppsDir)
#
if (not os.path.isdir(ppsDir)):
    print (f"PPS folder {ppsDir} does not exist. Will create it")
    logging.info (f"PPS folder {ppsDir} does not exist. Will create it")
    os.mkdir(ppsDir)
#
#
#%%
# === CIFBUILD ====
#
# adding the path where the testing CCF file is:
#
os.environ["SAS_ODF"] = odfdir
os.environ["SAS_CCFPATH"] = "/ccf/pub"
cif_file = "ccf.cif"
#
if (test_ccf):
    logging.info(f"Will use testing CCF files found in /xdata/xcaldata/XMM/IVAN/ccfdev.")    
    os.environ["SAS_CCFPATH"] = "/ccf/valid:/xdata/xcaldata/XMM/IVAN/ccfdev"
    cif_file = f"ccf_test_{subfolder}.cif"
os.chdir(odfdir)
#
cifbuild = f"cifbuild calindexset={cif_file} withccfpath=yes ccfpath=$SAS_CCFPATH fullpath=yes -w 1"
#
status = run_command(cifbuild)
if (status != 0):
    #
    # remove ODF_DIR/*.SAS and try again
    #
    sumFile = glob.glob(f"{odfdir}/*.SAS")
    if (len(sumFile) > 0):
        logging.warning(f"Removing {sumFile} as CIFBUILD failed. Trying again.")
        print(f"Removing {sumFile} as CIFBUILD failed. Trying again.")
        os.remove(sumFile[0])
    #
    # let's try again
    #
    status = run_command(cifbuild)
    if (status != 0):
        logging.error(f"Execution of cifbuild failed.")
        raise Exception
#
ccffile = os.path.join(odfdir,cif_file)
if (not os.path.isfile(ccffile)):
    print (f"No ccf.cif file found in ODF folder {odfdir}.")
    logging.error (f"No ccf.cif file found in ODF folder {odfdir}.")
    raise FileNotFoundError
#
os.environ["SAS_CCF"] = ccffile
logging.info(f"Set SAS_CCF to {ccffile}.")
#%%
# Have to run ODFINGEST as I moved the files
#
print ("*** ODFINGEST")
odfcomm = f"odfingest withodfdir=yes odfdir=\"{odfdir}\""
status = run_command(odfcomm)
if (status != 0):
    logging.error (f"ODFINGEST failed.")
    raise Exception
#
# parse the ODF Summary SAS file to only select observations with SmallWindow mode for PN
#
pnLW = False
sumFile = glob.glob(f"{odfdir}/*.SAS")

if (len(sumFile) != 1):
    print (f"No ODF summary file *.SAS found in ODF folder {odfdir}. Cannot continue!")
    logging.error(f"No ODF summary file *.SAS found in {odfdir}. Cannot continue!")
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
logging.info(f"Set SAS_ODF to {odfdir}.")
#%%
# ======== EPCHAIN =============
#
os.chdir(ppsDir)
#
task = "epchain"
if ("nocti" in subfolder):
    logging.warning(f"Epchain will be run with no long-term CTI correction!")
    xtask = f"{task} odf={odfdir} withctilongterm=N odfaccess=all backgroundtres=N "
else:
    logging.warning(f"Epchain will be run with the long-term CTI correction!")
    xtask = f"{task} odf={odfdir} withctilongterm=Y odfaccess=all backgroundtres=N "
#
print (f"*** RUNNING {xtask}")
status = run_command(xtask)
if (status != 0):
    print(f"Execution of {xtask} failed.")
    logging.error(f"Execution of {xtask} failed.")
    raise Exception
#
#%%
haveEvlist = False
evfind = glob.glob(f"{ppsDir}/*PIEVLI*")
tmax = 0.0
for xevl in evfind:
    hdu = fits.open(xevl)
    ontime = hdu['EVENTS'].header['ONTIME']
    if (ontime >= tmax):
        haveEvlist = True
        tmax = ontime
        evlist = xevl
    hdu.close()
#        
if not haveEvlist:
    logging.error(f"Event list for {obsid} not found.")
    raise FileNotFoundError
#
#%%
#
# Now generate the spectra, will reuse the regions from previous run
#
if (not os.path.isfile(evlist)):
    print(f"Event list {evlist} not available!")
    logging.error(f"Event list {evlist} not available.")
    raise FileNotFoundError
else:
    print(f"Will use event list {evlist}.")
    logging.info(f"Will use event list {evlist}.")    
#
print ("*** Generating spectra")
#
# now loop over all 12 CCDs in EPN
#
bore_ccds = ["03","04","09","10"]
#
for i in range(12):
    ccd = f"{i+1:02}"
    print (f"Processing CCD # {ccd}")
    #
    # PATTERN 0
    #
    evsel = "(PATTERN==0) && (PAT_SEQ==0) && #XMMEA_EP"
    sel1 = f"(CCDNR == {i+1})"
    if (ccd in  bore_ccds):
        sel1 = f"(CCDNR == {i+1}) && (RAWY >= 100) && (RAWY <= 120)"
        logging.info (f"Extracting a spectrum from a restricted RAWY in [100,120] for CCD {ccd}")
    spec_file = f"pn_{ccd}_pat00_spec5.fits"
    evt_file = f"pn_{ccd}_pat00_events.fits"
    command = "evselect " + \
        f"table={evlist} energycolumn='PI' withspectrumset=yes" + \
        f" expression='{evsel} && {sel1}'" + \
        " withspecranges=yes specchannelmin=0 specchannelmax=20479" + \
        f" spectrumset={spec_file} spectralbinsize=5" + \
        f" withfilteredset=yes filteredset={evt_file}"
    #    
    status = run_command(command)
    if (status != 0):
        print(f"evselect for pattern 0 failed for CCD {ccd}")
        logging.error(f"evselect for pattern 0 failed for CCD {ccd}")
        raise Exception
    #
    # PATTERN in [1:4] combined
    #
    evsel = "(PATTERN <= 4) && (PATTERN >= 1) && #XMMEA_EP"
    sel1 = f"(CCDNR == {i+1})"
    if (ccd in  bore_ccds):
        sel1 = f"(CCDNR == {i+1}) && (RAWY >= 100) && (RAWY <= 120)"
        logging.info (f"Extracting a spectrum from a restricted RAWY in [100,120] for CCD {ccd}")
    spec_file = f"pn_{ccd}_pat14_spec5.fits"
    evt_file = f"pn_{ccd}_pat14_events.fits"
    command = "evselect " + \
        f"table={evlist} energycolumn='PI' withspectrumset=yes" + \
        f" expression='{evsel} && {sel1}'" + \
        " withspecranges=yes specchannelmin=0 specchannelmax=20479" + \
        f" spectrumset={spec_file} spectralbinsize=5" + \
        f" withfilteredset=yes filteredset={evt_file}"
    #    
    status = run_command(command)
    if (status != 0):
        print(f"evselect for patterns 1 to 4 failed for CCD {ccd}")
        logging.error(f"evselect for patterns 1 to 4 failed for CCD {ccd}")
        raise Exception
    if (patts):
        #
        for xpat in range(1,5):
            # split per pattern
            evsel = f"(PATTERN=={xpat}) && #XMMEA_EP"
            spec_file = f"pn_{ccd}_pat{xpat:02}_spec5.fits"
            evt_file = f"pn_{ccd}_pat{xpat:02}_events.fits"
            command = "evselect " + \
                f"table={evlist} energycolumn='PI' withspectrumset=yes" + \
                f" expression='{evsel} && {sel1}'" + \
                " withspecranges=yes specchannelmin=0 specchannelmax=20479" + \
                f" spectrumset={spec_file} spectralbinsize=5" + \
                f" withfilteredset=yes filteredset={evt_file}"
                #    
            status = run_command(command)
            if (status != 0):
                print(f"evselect for pattern {xpat} failed for CCD {ccd}")
                logging.error(f"evselect for pattern {xpat} failed for CCD {ccd}")
                raise Exception
#
#
# now copy back the files from the node
#
if (copy_to_node):
    dest_dir = os.path.join(wdir,obsid,subfolder)
    # check if the destination folder exist, then remove it
    if (os.path.isdir(dest_dir)):
        print (f"Found an already existing folder {dest_dir}. Will remove it!")
        logging.warning (f"Found an already existing folder {dest_dir}. Will remove it!")
        shutil.rmtree(dest_dir)
    shutil.copytree(ppsDir,dest_dir)
    print (f"Folder {ppsDir} copied to {dest_dir}")
    logging.info (f"Folder {ppsDir} copied to {dest_dir}")    
#    
logging.info (f"pnlw_nofilt_proc done on {obsid}, products saved in {subfolder}")
print ("All done")# will attache a canned response file for PN in LW, single events
