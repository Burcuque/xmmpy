#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 09 12:00:00 2018

Step 05 of PN SW processing:  extract source and background spectra

Wil only process PN SmallWindow mode.

Will use the PN GTI and region files from the standard procesing
 
The output will be in a foder called PN_SW

Using epchain

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
parser = argparse.ArgumentParser(description='XMM SAS step 05 for PN SW, extract spectra')
parser.add_argument('obsid', type=str,
                    help='The OBSID to process')
#
args = parser.parse_args()
#
#%%
#
obsid = args.obsid
odfdir = os.path.join(os.getcwd(),obsid)
#
if (not os.path.isdir(odfdir)):
    print (f"The ODF folder {odfdir} does not exist! Cannot continue")
    raise FileNotFoundError
#
ppsDir = os.path.join(odfdir,"PN_SW")
#
if (not os.path.isdir(ppsDir)):
    print (f"PPS folder {ppsDir} does not exist. This means there are no event lists produced")
    print (f"Have you run pn_sw_step02.py ?")
    raise FileNotFoundError
#
logging.basicConfig(level=logging.DEBUG,
                    format='%(asctime)s %(levelname)s %(message)s',
                    filename=os.path.join(ppsDir,'pn_sw_step05.log'),
                    filemode='w')
#
os.environ["SAS_ODF"] = odfdir
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
evlist = os.path.join(ppsDir,'pn_evlist_clean.fits')
if (not os.path.isfile(evlist)):
    print(f"GTI filtered event list {evlist} not available. Have you run pn_sw_step03.py ?")
    logging.error(f"GTI filtered event list {evlist} not available.")
    raise FileNotFoundError
#
#%%
#
# check if there is a region file with source and background
#
reg_file = "pn_sw_regions.reg"
if (not os.path.isfile(reg_file)):
    print (f"No {reg_file} file found in folder {ppsDir}")
    print ("   The file have to have two regions with source and background")
    logging.error(f"No regions file pn_sw_regions.reg found in {ppsDir}")
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
command = "evselect table={} withspectrumset=yes spectrumset={}".format(evlist,src_spec) +  \
" energycolumn=PI spectralbinsize=5 withspecranges=yes specchannelmin=0 specchannelmax={}".format(spec_chan_max) +  \
" expression='{} && ((X,Y) IN {})'".format(expr1,src_reg)

print ("*** Extracting spectrum in the source region {}".format(src_reg))
try:
    result = subprocess.run(command, shell=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
    retcode=result.returncode
    if retcode != 0:
        print("evselect was terminated by signal", -retcode, file=sys.stderr)
        logging.error("evselect was terminated by signal: {} \n {}".format(-retcode,result.stdout.decode()))
        raise Exception
    else:
        print("evselect returned", retcode, file=sys.stderr)
        logging.info("evselect returned {}, \n {}".format(retcode,result.stdout.decode()))
except OSError as e:
    print("Execution of evselect failed:", e, file=sys.stderr)
    logging.error("Execution of evselect failed: {}".format(e))

#
command = "backscale spectrumset={} badpixlocation={}".format(src_spec,evlist)
print ("*** Backscale the source spectrum")
try:
    result = subprocess.run(command, shell=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
    retcode=result.returncode
    if retcode != 0:
        print("backscale was terminated by signal", -retcode, file=sys.stderr)
        logging.warning("backscale was terminated by signal: {}".format(-retcode))
        raise Exception
    else:
        print("backscale returned", retcode, file=sys.stderr)
        logging.info("backscale returned {}, \n {}".format(retcode,result.stdout.decode()))
except OSError as e:
    print("Execution of backscale failed:", e, file=sys.stderr)
    logging.error("Execution of backscale failed: {}".format(e))
#%%
# extract the background spectrum
#        
#
bg_spec = "pn_bkg_spectrum.fits"
command = "evselect table={} withspectrumset=yes spectrumset={}".format(evlist,bg_spec) +  \
" energycolumn=PI spectralbinsize=5 withspecranges=yes specchannelmin=0 specchannelmax={}".format(spec_chan_max) +  \
" expression='{} && ((X,Y) IN {})'".format(expr1,bkg_reg)
print ("*** Extracting spectrum in the background region {}".format(bkg_reg))
try:
    result = subprocess.run(command, shell=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
    retcode=result.returncode
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
# now the backscale for the background
#
command = "backscale spectrumset={} badpixlocation={}".format(bg_spec,evlist)
print ("*** Backscale the background spectrum")
try:
    result = subprocess.run(command, shell=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
    retcode=result.returncode
    if retcode < 0:
        print("backscale was terminated by signal", -retcode, file=sys.stderr)
        logging.warning("backscale was terminated by signal: {}".format(-retcode))
    else:
        print("backscale returned", retcode, file=sys.stderr)
        logging.info("backscale returned {}, \n {}".format(retcode,result.stdout.decode()))
except OSError as e:
    print("Execution of backscale failed:", e, file=sys.stderr)
    logging.error("Execution of backscale failed: {}".format(e))
#%%
# now generate the RMF (can take time)
command = "rmfgen spectrumset={} rmfset=pn.rmf".format(src_spec)
print ("*** RMF generation for PN")
try:
    result = subprocess.run(command, shell=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
    retcode=result.returncode
    if retcode < 0:
        print("rmfgen was terminated by signal", -retcode, file=sys.stderr)
        logging.warning("rmfgen was terminated by signal: {} \n {}".format(-retcode,result.stdout.decode()))
    else:
        print("rmfgen returned", retcode, file=sys.stderr)
        logging.info("rmfgen returned {}, \n {}".format(retcode,result.stdout.decode()))
except OSError as e:
    print("Execution of rmfgen failed:", e, file=sys.stderr)
    logging.error("Execution of rmfgen failed: {}".format(e))
#%%
# now generate the ARF
command = "arfgen spectrumset={} arfset=pn.arf withrmfset=yes rmfset=pn.rmf".format(src_spec) + \
    " badpixlocation={} detmaptype=psf".format(evlist)
print ("*** ARF generation for PN")
try:
    result = subprocess.run(command, shell=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
    retcode=result.returncode
    if retcode < 0:
        print("arfgen was terminated by signal", -retcode, file=sys.stderr)
        logging.warning("arfgen was terminated by signal: {}".format(-retcode))
    else:
        print("arfgen returned", retcode, file=sys.stderr)
        logging.info("arfgen returned {}, \n {}".format(retcode,result.stdout.decode()))
except OSError as e:
    print("Execution of arfgen failed:", e, file=sys.stderr)
    logging.error("Execution of arfgen failed: {}".format(e))
#
#%%
#
# Use specgroup to add background, ARF and RMF files in headers. Without any grouping first
#
command = "specgroup spectrumset={} addfilenames=yes rmfset=pn.rmf".format(src_spec) + \
    " arfset=pn.arf backgndset={} groupedset=pn_spectrum_grp0.fits".format(bg_spec)
print ("*** Link bkg, ARF and RMF filenames to the PN spectrum")

try:
    result = subprocess.run(command, shell=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
    retcode=result.returncode
    if retcode < 0:
        print("specgroup was terminated by signal", -retcode, file=sys.stderr)
        logging.warning("specgroup was terminated by signal: {}".format(-retcode))
    else:
        print("specgroup returned", retcode, file=sys.stderr)
        logging.info("specgroup returned {}, \n {}".format(retcode,result.stdout.decode()))
except OSError as e:
    print("Execution of specgroup failed:", e, file=sys.stderr)
    logging.error("Execution of specgroup failed: {}".format(e))
#%%
# now group/bin the spectrum
#
min_cts = 25
command = "specgroup spectrumset={} mincounts={} oversample=3 rmfset=pn.rmf".format(src_spec,min_cts) + \
    " arfset=pn.arf backgndset={} groupedset=pn_spectrum_grp{}.fits".format(bg_spec,min_cts)
print ("*** Group the PN spectrum in bins of at least {} counts".format(min_cts))
try:
    result = subprocess.run(command, shell=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
    retcode=result.returncode
    if retcode < 0:
        print("specgroup was terminated by signal", -retcode, file=sys.stderr)
        logging.warning("specgroup was terminated by signal: {}".format(-retcode))
    else:
        print("specgroup returned", retcode, file=sys.stderr)
        logging.info("specgroup returned {}, \n {}".format(retcode,result.stdout.decode()))
except OSError as e:
    print("Execution of specgroup failed:", e, file=sys.stderr)
    logging.error("Execution of specgroup failed: {}".format(e))

print ("All done")