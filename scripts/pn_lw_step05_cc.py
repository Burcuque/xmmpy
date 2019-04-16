#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 09 12:00:00 2018

Step 05 of PN LW processing:  extract spetra from a spatial filter frmo Michael Smith

Output per CCD

Wil only process PN Large Window mode.

Will use the PN GTI and region files from the standard procesing
 
The output will be in a foder called PN_LW

Using evselect

@author: ivaltchanov
"""
import os
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
                print("Execution was terminated by signal", -retcode, file=sys.stderr)
            logging.warning("Execution was terminated by signal: {} \n {}".format(-retcode,result.stdout.decode()))
        else:
            if (verbose):
                print("Execution returned", retcode, file=sys.stderr)
            logging.info("Execution returned {}, \n {}".format(retcode,result.stdout.decode()))
    except OSError as e:
        print("Execution failed:", e, file=sys.stderr)
        logging.error("Execution of evselect failed: {}".format(e))
    return retcode
#
#%%
#
# get the arguments
parser = argparse.ArgumentParser(description='XMM SAS step 05 for PN LW, making spectra')
parser.add_argument('obsid', type=str,
                    help='The OBSID to process')
parser.add_argument('expo_name', type=str,
                    help='The OBSID exposure to process')
#
args = parser.parse_args()
#
#%%
#
obsid = args.obsid
nexpo = args.expo_name
odfdir = os.path.join(os.getcwd(),obsid)
#
if (not os.path.isdir(odfdir)):
    print (f"The ODF folder {odfdir} does not exist! Cannot continue")
    raise FileNotFoundError
#
ppsDir = os.path.join(odfdir,"PN_LW")
#
if (not os.path.isdir(ppsDir)):
    print (f"PPS folder {ppsDir} does not exist. This means there are no event lists produced")
    print (f"Have you run pn_sw_step02.py ?")
    raise FileNotFoundError
#
logging.basicConfig(level=logging.DEBUG,
                    format='%(asctime)s %(levelname)s %(message)s',
                    filename=os.path.join(ppsDir,'pn_lw_step05.log'),
                    filemode='w')
#
os.environ["SAS_ODF"] = odfdir
os.environ['SAS_CCFPATH'] = '/xdata/ccf/pub'
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
#
os.chdir(ppsDir)

#%%
# Select the cleaned event lists
#
evlist = os.path.join(ppsDir,f'pn_evlist_clean_{nexpo}.fits')
if (not os.path.isfile(evlist)):
    print(f"GTI filtered event list {evlist} not available. Have you run pn_lw_step03.py ?")
    logging.error(f"GTI filtered event list {evlist} not available.")
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
logging.info ("pn_lw_step05 done")
#
