#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 09 12:00:00 2018

Full workflow for EPIC-PN processing an XMM observation, filtering only single events, i.e.
PATTERN == 0

Will use an already generated clean event list

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
parser = argparse.ArgumentParser(description='XMM SAS workflow for EPN in SW with new CCF')
parser.add_argument('obsid', type=str,
                    help='The OBSID to process')
parser.add_argument('subfolder', type=str,
                    help='The subfolder to process: can be "cti49", "n_cti" or "xmmsas_20180620_1732-17.0.0"')
parser.add_argument('--odfdir', type=str, default=os.getcwd(),
                    help='Where the ODF files are')
#
args = parser.parse_args()
#
#%%
#
obsid = args.obsid
subf = args.subfolder
odfdir = os.path.join(os.path.abspath(args.odfdir),obsid)
#
if (not os.path.isdir(odfdir)):
    print ("The ODF folder {} does not exist! Cannot continue".format(odfdir))
    raise FileNotFoundError
#
ppsDir = os.path.join(odfdir,subf)
#
if (not os.path.isdir(ppsDir)):
    print (f"PPS folder {ppsDir} does not exist. Cannot continue")
    raise FileNotFoundError
#
logging.basicConfig(level=logging.DEBUG,
                    format='%(asctime)s %(levelname)s %(message)s',
                    filename=os.path.join(ppsDir,f'pn_proc_singles_{subf}.log'),
                    filemode='w')
#
os.environ["SAS_ODF"] = odfdir
os.environ["SAS_CCFPATH"] = "/ccf/valid:/xdata/xcaldata/XMM/IVAN/ccfdev"
cif_file = "ccf_test.cif"
ccffile = os.path.join(odfdir,cif_file)
if (not os.path.isfile(ccffile)):
    print (f"No ccf.cif file found in ODF folder {odfdir}. Cannot continue.")
    raise FileNotFoundError
#
os.environ["SAS_CCF"] = ccffile
logging.info("Set SAS_ODF to {}".format(odfdir))
logging.info("Set SAS_CCF to {}".format(ccffile))
os.chdir(ppsDir)
#
#%%
#
# find the filtered event list, which will be used to generate spectra with 
# additional filtering on pattern and flag
#
filt_evlist = f"pn_evlist_clean_{subf}.fits"
if ('xmmsas' in subf):
    filt_evlist = f"pn_evlist_clean.fits"

if (not os.path.isfile(filt_evlist)):
    logging.error("Cannot find cleaned event list")
    raise FileNotFoundError
print ("*** Generating spectra")
spec_chan_max = 20479
expr1 = "(FLAG==0) && (PATTERN == 0)"
#
# read the region file it must have a name source_region_[inst].txt
# and the background region must have a name bg_region_[inst].txt
#
src_reg_file = os.path.join(ppsDir,'../xmmsas_20180620_1732-17.0.0/regions','src_region_pn.reg')
bkg_reg_file = os.path.join(ppsDir,'../xmmsas_20180620_1732-17.0.0/regions','bkg_region_pn.reg')
#
if (not (os.path.isfile(src_reg_file) and os.path.isfile(bkg_reg_file))):
    print ("Missing region file for source or background")
    print ("Please create two files with names {} (source region) and {} (background region) ".format(src_reg_file,bkg_reg_file))
    print ("  each file should contain one line with the following content:")
    print ("    CIRCLE(DETX_CEN,DETY_CEN,DET_RAD)")
    print ("   where DETX_CEN, DETY_CEN is the centre in detector coordinates and DET_RAD is the radius.")
    raise FileNotFoundError
#
# extract the source spectrum from the region
with open(src_reg_file) as regfile:
    reg_line = regfile.readline()
src_reg = reg_line.strip().split()[0]
#
src_spec = "pn_src_spectrum_singles_cti49.fits"
command = f"evselect table={filt_evlist} withspectrumset=yes spectrumset={src_spec}" +  \
    f" energycolumn=PI spectralbinsize=5 withspecranges=yes specchannelmin=0 specchannelmax={spec_chan_max}" +  \
    f" expression='{expr1} && ((X,Y) IN {src_reg})'"
print ("*** Extracting spectrum in the source region")
status = run_command(command)
if (status != 0):
    raise Exception
#
# now calculate the backscale for the source region
#%%
command = f"backscale spectrumset={src_spec} badpixlocation={filt_evlist}"
print ("*** Backscale the source spectrum")
status = run_command(command)
if (status != 0):
    raise Exception
#%%
# extract the background spectrum
#        
with open(bkg_reg_file) as regfile:
    reg_line = regfile.readline()
bg_reg = reg_line.strip().split()[0]
#
bg_spec = "pn_bkg_spectrum_singles_cti49.fits"
command = f"evselect table={filt_evlist} withspectrumset=yes spectrumset={bg_spec}" + \
    f" energycolumn=PI spectralbinsize=5 withspecranges=yes specchannelmin=0 specchannelmax={spec_chan_max}" +  \
    f" expression='{expr1} && ((X,Y) IN {bg_reg})'"
print ("*** Extracting spectrum in the background region")
status = run_command(command)
if (status != 0):
    raise Exception
#
# now the backscale for the background
#
command = f"backscale spectrumset={bg_spec} badpixlocation={filt_evlist}"
print ("*** Backscale the background spectrum")
status = run_command(command)
if (status != 0):
    raise Exception
#%%
# now generate the RMF (can take time)
command = f"rmfgen spectrumset={src_spec} rmfset=pn_singles_cti49.rmf"
print ("*** RMF generation")
status = run_command(command)
if (status != 0):
    raise Exception
#%%
# now generate the ARF
command = f"arfgen spectrumset={src_spec} arfset=pn_singles_cti49.arf" + \
    " withrmfset=yes rmfset=pn_singles_cti49.rmf" + \
    " badpixlocation={filt_evlist} detmaptype=psf"
print ("*** ARF generation")
status = run_command(command)
if (status != 0):
    raise Exception
#
#%%
#
# Use specgroup to add background, ARF and RMF files in headers. Without any grouping first
#
command = f"specgroup spectrumset={src_spec} addfilenames=yes rmfset=pn_singles_cti49.rmf" + \
    f" arfset=pn_singles_cti49.arf backgndset={bg_spec} groupedset=pn_spectrum_singles_cti49_grp0.fits"
print ("*** Link bkg, ARF and RMF filenames to the spectrum")
status = run_command(command)
if (status != 0):
    raise Exception

print ("All done")