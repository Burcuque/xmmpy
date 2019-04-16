#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 09 12:00:00 2018

Step 05 of PN LW processing:  extract spetra from a source region and background


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
cwd = os.getcwd()
odfdir = os.path.join(cwd,obsid)
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
print ("*** Generating spectra")
#
reg_file = f"{cwd}/regions/{obsid}_pn_regions.reg"
if (not os.path.isfile(reg_file)):
    print (f"No {reg_file} file found in current folder")
    print ("   The file have to have two regions with source and background")
    logging.error (f"No {reg_file} file found in current folder")
    raise FileNotFoundError
#
# now parse the regions file
with open(reg_file,'r') as regs:
    lines = regs.readlines()
for line in lines:
    qq = line.split()
    if (('circle' in qq[0]) and ('source' in line)):
        src_reg = qq[0]
    if (('circle' in qq[0]) and not ('source' in line)):
        bkg_reg = qq[0]
#
# PN specific filtering
spec_chan_max = 20479
expr1 = "(FLAG==0) && (PATTERN<=4)"
#
src_spec = "pn_src_spectrum.fits"
command = f"evselect table={evlist} withspectrumset=yes spectrumset={src_spec}" +  \
f" energycolumn=PI spectralbinsize=5 withspecranges=yes specchannelmin=0 specchannelmax={spec_chan_max}" +  \
f" expression='{expr1} && ((X,Y) IN {src_reg})'"

print (f"*** Extracting spectrum in the source region {src_reg}")
status = run_command(command,verbose=True)
if (status != 0):
    raise Exception
#
command = f"backscale spectrumset={src_spec} badpixlocation={evlist}"
print ("*** Backscale the source spectrum")
status = run_command(command,verbose=True)
if (status != 0):
    raise Exception
#
# extract the background spectrum
#        
#
bg_spec = "pn_bkg_spectrum.fits"
command = f"evselect table={evlist} withspectrumset=yes spectrumset={bg_spec}" +  \
f" energycolumn=PI spectralbinsize=5 withspecranges=yes specchannelmin=0 specchannelmax={spec_chan_max}" +  \
f" expression='{expr1} && ((X,Y) IN {bkg_reg})'"
print ("*** Extracting spectrum in the background region {}".format(bkg_reg))
status = run_command(command,verbose=True)
if (status != 0):
    raise Exception
#
# now the backscale for the background
#
command = f"backscale spectrumset={bg_spec} badpixlocation={evlist}"
print ("*** Backscale the background spectrum")
status = run_command(command,verbose=True)
if (status != 0):
    raise Exception
#
# now generate the RMF (can take time)
#
command = f"rmfgen spectrumset={src_spec} rmfset=pn.rmf"
print ("*** RMF generation for PN")
status = run_command(command,verbose=True)
if (status != 0):
    raise Exception
# now generate the ARF
command = f"arfgen spectrumset={src_spec} arfset=pn.arf withrmfset=yes rmfset=pn.rmf" + \
    f" badpixlocation={evlist} detmaptype=psf"
print ("*** ARF generation for PN")
status = run_command(command,verbose=True)
if (status != 0):
    raise Exception
#
#%%
#
# Use specgroup to add background, ARF and RMF files in headers. Without any grouping first
#
command = f"specgroup spectrumset={src_spec} addfilenames=yes rmfset=pn.rmf" + \
    f" arfset=pn.arf backgndset={bg_spec} groupedset=pn_spectrum_grp0.fits"
print ("*** Link bkg, ARF and RMF filenames to the PN spectrum")
status = run_command(command,verbose=True)
if (status != 0):
    raise Exception
#
logging.info ("pn_lw_step05 done")
#
print ("All done")# will attache a canned response file for PN in LW, single events
