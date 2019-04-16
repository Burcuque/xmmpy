#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 09 12:00:00 2018

PN LW processing:  running epchain with the new testing CCF.

Adapted for spatial filtering for CalClosed data

Will only process PN Large Window mode.
 
The output will be in a foder called cti49

Using epchain

@author: ivaltchanov
"""
import os
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

os.environ['SAS_CCFPATH'] = '/xdata/ccf/pub'
#%%
#
# get the arguments
parser = argparse.ArgumentParser(description='XMM SAS full processing chain for PN LW')
parser.add_argument('obsid', type=str,
                    help='The OBSID to process')
parser.add_argument('expo', type=str,default="S003",
                    help='The exposure to process')
parser.add_argument('--wdir', type=str,default=os.getcwd(),
                    help='The working directory (root folder)')
#
args = parser.parse_args()
#
#%%
#
obsid = args.obsid
wdir = os.path.abspath(args.wdir)
odfdir = os.path.join(wdir,obsid)
nexpo = args.expo
#
if (not os.path.isdir(odfdir)):
    print (f"The ODF folder {odfdir} does not exist! Cannot continue.")
    raise FileNotFoundError
#
ppsDir = os.path.join(odfdir,"cti49")
#
if (not os.path.isdir(ppsDir)):
    print (f"PPS folder {ppsDir} does not exist. Will create it")
    os.mkdir(ppsDir)
#
logging.basicConfig(level=logging.DEBUG,
                    format='%(asctime)s %(levelname)s %(message)s',
                    filename=os.path.join(ppsDir,'pn_lw_cc_proc_cti49.log'),
                    filemode='w')
#
os.environ["SAS_ODF"] = odfdir

#%%
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
#
#%%
# === CIFBUILD ====
#
# adding the path where the testing CCF file is:
#
os.environ["SAS_CCFPATH"] = "/ccf/valid:/xdata/xcaldata/XMM/IVAN/ccfdev"
#
os.chdir(odfdir)
#
cif_file = "ccf_test.cif"
cifbuild = f"cifbuild calindexset={cif_file} withccfpath=yes ccfpath=$SAS_CCFPATH fullpath=yes -w 1"
#
status = run_command(cifbuild)
if (status != 0):
    raise Exception
#
ccffile = os.path.join(odfdir,cif_file)
if (not os.path.isfile(ccffile)):
    print (f"No ccf.cif file found in ODF folder {odfdir}.")
    raise FileNotFoundError
#
os.environ["SAS_CCF"] = ccffile
logging.info(f"Set SAS_CCF to {ccffile}.")
#
# Have to run ODFINGEST as I moved the files
#
print ("*** ODFINGEST")
odfcomm = f"odfingest withodfdir=yes odfdir=\"{odfdir}\""
status = run_command(odfcomm)
if (status != 0):
    raise Exception
#
logging.info(f"Set SAS_ODF to {odfdir}.")
#%%
# ======== EPCHAIN =============
#
os.chdir(ppsDir)
#
task = "epchain"
xtask = f"{task} odf={odfdir} withctilongterm=Y odfaccess=all backgroundtres=N "
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
for evlist in evfind:
    if (nexpo in evlist):
        haveEvlist = True
        break
if not haveEvlist:
    logging.error(f"Event list for {obsid} and exposure {nexpo} not found.")
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
#
print ("*** Generating spectra")
#
respfile = f"/home/ivaltchanov/IVAN/PN_LW/epn_e3_lw20_sY9_v17.0.rmf"
#
evsel = "(PATTERN==0) && (PAT_SEQ==0) && #XMMEA_EP"
#
# now loop over all 12 CCDs in EPN
#
for i in range(12):
    ccd = f"{i+1:02}"
    print (f"Processing CCD # {ccd}")
    #
    mnk_mask =  f"/xdata/xcaldata/XMM/PN/CTI/mask/refmask_{ccd}_mnk_lmap.fits"
    #
    sel1 = f"mask({mnk_mask},1,1,RAWX,RAWY) && (CCDNR == {i+1})"
    spec_file = f"pn_{nexpo}_{ccd}_all_spec5.fits"
    command = "evselect " + \
        f"table={evlist} energycolumn='PI' withspectrumset=yes" + \
        f" expression='{evsel} && {sel1}'" + \
        " withspecranges=yes specchannelmin=0 specchannelmax=20479" + \
        f" spectrumset={spec_file} spectralbinsize=5"
    status = run_command(command)
    if (status != 0):
        raise Exception
    #
    # adding some keywords needed for XSPEC
    #
    hdu = fits.open(spec_file)
    hdu[1].header['BACKFILE'] = "NONE"
    hdu[1].header['ANCRFILE'] = "NONE"
    hdu[1].header['RESPFILE'] = respfile
    hdu.writeto(spec_file,overwrite=True)
#    #
#    # now loop over RAWY in 20 pixels
#    #
#    for j in np.arange(100,200,20):
#        sel2 = f"(RAWY > {j} && RAWY <= {j+20})"
#        spec_file = f"pn_{nexpo}_{ccd}_{j}_spec5.fits"
#        command = "evselect " + \
#            f"table={evlist} energycolumn='PI' withspectrumset=yes" + \
#            f" expression='{evsel} && {sel1} && {sel2}'" + \
#            " withspecranges=yes specchannelmin=0 specchannelmax=20479" + \
#            f" spectrumset={spec_file} spectralbinsize=5"
#        status = run_command(command)
#        if (status != 0):
#            raise Exception
    # For CCD 4 we add the boresight with RAWY in [180,200]
    if (ccd == "04"):
        print (f"Extracting the boresight spectrum from CCD {ccd}")
        sel2 = f"(RAWY > 180 && RAWY <= 200)"
        spec_file = f"pn_{nexpo}_{ccd}_180_spec5.fits"
        command = "evselect " + \
            f"table={evlist} energycolumn='PI' withspectrumset=yes" + \
            f" expression='{evsel} && {sel1} && {sel2}'" + \
            " withspecranges=yes specchannelmin=0 specchannelmax=20479" + \
            f" spectrumset={spec_file} spectralbinsize=5"
        status = run_command(command)
        if (status != 0):
            raise Exception
        #
        # adding some keywords needed for XSPEC
        #
        hdu = fits.open(spec_file)
        hdu[1].header['BACKFILE'] = "NONE"
        hdu[1].header['ANCRFILE'] = "NONE"
        hdu[1].header['RESPFILE'] = respfile
        hdu.writeto(spec_file,overwrite=True)
#
logging.info ("pn_cc_cti49_proc done")
#
print ("All done")# will attache a canned response file for PN in LW, single events
