#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 09 12:00:00 2018

Step 03 of PN LW processing:  extract spectra from a spatial filter with mask 
 from Michael Smith

Output per CCD

Using the even list with no testing CTI in folder <subfolder>

Will only process PN Large Window mode.

No GTI filtering, will use the full PIEVL event list
 
Using evselect

@author: ivaltchanov
"""
import os
import glob
import subprocess
import sys
import logging
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
parser = argparse.ArgumentParser(description='XMM SAS step 03 for PN LW, making spectra from full CCDs')
parser.add_argument('obsid', type=str,
                    help='The OBSID to process')
parser.add_argument('subfolder', type=str,
                    help='The subfolder where to save the products')
#parser.add_argument('expo_name', type=str,
#                    help='The OBSID exposure to process')
parser.add_argument('--wdir', type=str,default=os.getcwd(),
                    help='The working directry (root folder)')
#
args = parser.parse_args()
#
#%%
#
home = os.path.expanduser('~')
obsid = args.obsid
wdir = args.wdir
#nexpo = args.expo_name
odfdir = os.path.join(wdir,obsid)
subfolder = args.subfolder
#
if (not os.path.isdir(odfdir)):
    print (f"The ODF folder {odfdir} does not exist! Cannot continue")
    raise FileNotFoundError
#
ppsDir = os.path.join(odfdir,subfolder)
#
if (not os.path.isdir(ppsDir)):
    print (f"PPS folder {ppsDir} does not exist. This means there are no event lists produced")
    raise FileNotFoundError
#
logging.basicConfig(level=logging.DEBUG,
                    format='%(asctime)s %(levelname)s %(message)s',
                    filename=os.path.join(ppsDir,f'pn_lw_step03_ccd_{subfolder}.log'),
                    filemode='w')
#
os.environ["SAS_ODF"] = odfdir
#os.environ['SAS_CCFPATH'] = '/xdata/ccf/pub'
os.environ["SAS_CCFPATH"] = f"/ccf/valid:{home}/XMM/CAL/ccfprod"

ccffile = os.path.join(odfdir,f"ccf_test_{subfolder}.cif")

if (not os.path.isfile(ccffile)):
    print (f"No ccf.cif file found in ODF folder {odfdir}")
    logging.error(f"No CCF file ccf.cif found in {odfdir}")
    raise FileNotFoundError
#
os.environ["SAS_CCF"] = ccffile
logging.info(f"Set SAS_ODF to {odfdir}")
logging.info(f"Set SAS_CCF to {ccffile}")
#
os.chdir(ppsDir)

#%%
# Select the cleaned event lists
#
evlists = glob.glob(f"{ppsDir}/*PIEVLI*")
nev = len(evlists)
haveEvlist = False
if (nev == 1):
    evlist = os.path.basename(evlists[0])
elif (nev > 1):
    print (f"More than one exposure is available(n={nev}), will pick up the scheduled one")
    logging.warning (f"More than one exposure is available(n={nev}), will pick up the scheduled one")
    for iev in evlists:
        if ('PNS' in iev):
            evlist = os.path.basename(iev)
            haveEvlist = True
            break
        # and unscheduled exposure
        if ((not haveEvlist) and ('PNU' in iev)):
            evlist = os.path.basename(iev)
else:
    print (f"EPIC-PN event list *PIEVLI not available. Have you run pn_lw_step02.py ?")
    logging.error(f"EPIC-PN event list *PIEVLI not available. Have you run pn_lw_step02.py ?")
    raise FileNotFoundError
#
#%%
#
print ("*** Generating spectra")
#
# will attache a canned response file for PN in LW, single events
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
    spec_file = f"pn_{ccd}_all_spec5.fits"
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
#
logging.info (f"pn_lw_step03_ccd_{subfolder} done")
#
