#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 14 16:58:36 2018

Create a spectrum after finding the source and background regions

@author: ivaltchanov
"""

import os
import subprocess
import sys

import logging
import argparse
#
# global folders
#
os.environ['SAS_CCFPATH'] = '/xdata/ccf/pub'
#%%
#
parser = argparse.ArgumentParser(description='XMMSAS workflow step 4 - spectra')
parser.add_argument('obsid', type=str,
                    help='The OBSID to process')
parser.add_argument('inst', type=str,
                    help='The instrument to process, can be "mos1", "mos2" or "pn"')
parser.add_argument('--odfdir', type=str, default=os.getcwd(),
                    help='Where the ODF files are')
args = parser.parse_args()
#
obsid = args.obsid
inst  = args.inst
odfdir = os.path.join(os.path.abspath(args.odfdir),obsid)
#
#%%
if (not os.path.isdir(odfdir)):
    print ("The ODF folder {} does not exist! Cannot continue".format(odfdir))
    raise FileNotFoundError
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
#
if (not os.path.isdir(ppsDir)):
    print ("PPS folder {} does not exist. This means step #01 was not run?".format(ppsDir))
    raise FileNotFoundError
#
logging.basicConfig(level=logging.DEBUG,
                    format='%(asctime)s %(levelname)s %(message)s',
                    filename=os.path.join(ppsDir,'proc_{}_step05.log'.format(inst)),
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
#%%
evlist = '{}_evlist_clean.fits'.format(inst)
if (not os.path.isfile(evlist)):
    logging.error("Cannot find cleaned event lists")
print ("*** Generating spectra for {}".format(inst))
if ('pn' in inst):
    spec_chan_max = 20479
    expr1 = "(FLAG==0) && (PATTERN<=4)"
else:
    spec_chan_max = 11999
    expr1 = "#XMMEA_EM && (PATTERN<=12)"
#
# read the region file it must have a name source_region_[inst].txt
# and the background region must have a name bg_region_[inst].txt
#
src_reg_file = os.path.join(ppsDir,'regions','src_region_{}.reg'.format(inst))
bkg_reg_file = os.path.join(ppsDir,'regions','bkg_region_{}.reg'.format(inst))
#
if (not (os.path.isfile(src_reg_file) and os.path.isfile(bkg_reg_file))):
    print ("Missing region file for source or background")
    print ("Please create two files with names {} (source region) and {} (background region) ".format(src_reg_file,bkg_reg_file))
    print ("  each file should contain one line with the following content:")
    print ("    CIRCLE(DETX_CEN,DETY_CEN,DET_RAD)")
    print ("   where DETX_CEN, DETY_CEN is the centre in detector coordinates and DET_RAD is the radius.")
    raise FileNotFoundError
#
#%%
# extract the source spectrum from the region
with open(src_reg_file) as regfile:
    reg_line = regfile.readline()
src_reg = reg_line.strip().split()[0]
#
src_spec = "{}_src_spectrum.fits".format(inst)
command = "evselect table={} withspectrumset=yes spectrumset={}".format(evlist,src_spec) +  \
" energycolumn=PI spectralbinsize=5 withspecranges=yes specchannelmin=0 specchannelmax={}".format(spec_chan_max) +  \
" expression='{} && ((X,Y) IN {})'".format(expr1,src_reg)
print ("*** Extracting spectrum in the source region")
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
# now calculate the backscale for the source region
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
with open(bkg_reg_file) as regfile:
    reg_line = regfile.readline()
bg_reg = reg_line.strip().split()[0]
#
bg_spec = "{}_bkg_spectrum.fits".format(inst)
command = "evselect table={} withspectrumset=yes spectrumset={}".format(evlist,bg_spec) +  \
" energycolumn=PI spectralbinsize=5 withspecranges=yes specchannelmin=0 specchannelmax={}".format(spec_chan_max) +  \
" expression='{} && ((X,Y) IN {})'".format(expr1,bg_reg)
print ("*** Extracting spectrum in the background region")
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
command = "rmfgen spectrumset={} rmfset={}.rmf".format(src_spec,inst)
print ("*** RMF generation for {}".format(inst))
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
command = "arfgen spectrumset={} arfset={}.arf withrmfset=yes rmfset={}.rmf".format(src_spec,inst,inst) + \
    " badpixlocation={} detmaptype=psf".format(evlist)
print ("*** ARF generation for {}".format(inst))
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
command = "specgroup spectrumset={} addfilenames=yes rmfset={}.rmf".format(src_spec,inst) + \
    " arfset={}.arf backgndset={} groupedset={}_spectrum_grp0.fits".format(inst,bg_spec,inst)
print ("*** Link bkg, ARF and RMF filenames to the {} spectrum".format(inst))
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
command = "specgroup spectrumset={} mincounts={} oversample=3 rmfset={}.rmf".format(src_spec,min_cts,inst) + \
    " arfset={}.arf backgndset={} groupedset={}_spectrum_grp{}.fits".format(inst,bg_spec,inst,min_cts)
print ("*** Group the {} spectrum in bins of at least {} counts".format(inst,min_cts))
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