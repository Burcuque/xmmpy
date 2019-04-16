#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 09 12:00:00 2018

Full workflow for EPIC processing an XMM observation with no CTI correction applied.

Wil only process PN SmallWindow mode.

Will use backgroundtres=N in epchain to speed up the processing

Will use the merged GTI and region files from the standard procesing

Will use the new EPN_CTI_0049_test.CCF for the long-term CTI correction

The output will be in a foder called cti49

@author: ivaltchanov
"""
import os
import subprocess
import sys
import logging
import glob
from astropy.io import fits

import argparse

os.environ['SAS_CCFPATH'] = '/xdata/ccf/pub'
#%%
#
# get the arguments
parser = argparse.ArgumentParser(description='XMM SAS workflow for EPN in SW with new CCF')
parser.add_argument('obsid', type=str,
                    help='The OBSID to process')
parser.add_argument('--newRun', default=False, action='store_true', 
                    help='New run or skip steps which were already done.')
parser.add_argument('--odfdir', type=str, default=os.getcwd(),
                    help='Where the ODF files are')
#
args = parser.parse_args()
#
#%%
#
obsid = args.obsid
newRun = args.newRun
odfdir = os.path.join(os.path.abspath(args.odfdir),obsid)
#skipDone = args.skip
#obsid = "0771000201"
#odfdir = f'/xdata/xcaldata/XMM/IVAN/PN_SW/ngc5548/{obsid}'
#skipDone = True
#
if (not os.path.isdir(odfdir)):
    print ("The ODF folder {} does not exist! Cannot continue".format(odfdir))
    raise FileNotFoundError
#
ppsDir = os.path.join(odfdir,"cti49")
#
if (not os.path.isdir(ppsDir)):
    print ("PPS folder {} does not exist. Will create it".format(ppsDir))
    newRun = True
    os.mkdir(ppsDir)
#
task = 'epchain'
#
logging.basicConfig(level=logging.DEBUG,
                    format='%(asctime)s %(levelname)s %(message)s',
                    filename=os.path.join(ppsDir,'pn_proc_cti49.log'.format(task)),
                    filemode='w')
#
#
os.environ["SAS_ODF"] = odfdir
#
os.environ["SAS_CCFPATH"] = "/ccf/valid:/xdata/xcaldata/XMM/IVAN/PN_SW/ccfdev"
#
os.chdir(odfdir)
#
cif_file = "ccf_test.cif"

cifbuild = f"cifbuild calindexset={cif_file} withccfpath=yes ccfpath=$SAS_CCFPATH fullpath=yes -w 1"

try:
    result = subprocess.run(cifbuild, shell=True,stderr=subprocess.STDOUT,stdout=subprocess.PIPE)
    retcode = result.returncode
    if retcode < 0:
        print("cifbuild was terminated by signal {}".format(-retcode))
        logging.warning("cifbuild was terminated by signal {}".format(-retcode))
    else:
        print("cifbuild returned {} \n {}".format(retcode,result.stdout.decode()))
        logging.info("cifbuild returned {} \n {}".format(retcode,result.stdout.decode()))
except OSError as e:
    print("Execution of cifbuild failed: {}".format(e))
    logging.error("Execution of cifbuild failed: {}".format(e))
#
#
ccffile = os.path.join(odfdir,cif_file)
if (not os.path.isfile(ccffile)):
    print ("No ccf.cif file found in ODF folder {}. This means step #01 was not run?".format(odfdir))
    raise FileNotFoundError
#
os.environ["SAS_CCF"] = ccffile
logging.info("Set SAS_ODF to {}".format(odfdir))
logging.info("Set SAS_CCF to {}".format(ccffile))
#%%
# Have to run ODFINGEST as I moved the files
#
print ("*** ODFINGEST")
os.chdir(odfdir)

sumFile = glob.glob("{}/*.SAS".format(odfdir))
if (len(sumFile) != 1 or newRun):
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
else:
            logging.info("Skipping odfingest as already done")
            print("Skipping odfingest as already done")
    
#%%
# parse the ODF Summary SAS file to only select observations with SmallWindow mode for PN
#
pnSW = False
with open(sumFile[0],'r') as sas:
    lines = sas.readlines()
for qline in lines:
    if ("MODE = PrimeSmallWindow" in qline):
        pnSW = True
        break
#
if (not pnSW):
    print (f"Not processing {obsid} as it is not with PN Small Window Mode")
    raise Exception
else:
    print (f"Processing {obsid} as it is with PN Small Window Mode")
#
#
#%%
os.chdir(ppsDir)
#
evfind = glob.glob(f"{ppsDir}/*PIEVLI*")
#
if ((len(evfind) < 1) or newRun):
    xtask = task + f" odf={odfdir} odfaccess=all backgroundtres=N"
    print ("*** RUNNING %s"%xtask)
    try:
        #result = subprocess.run(xtask, shell=True)
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
else:
            logging.info("epchain skipped as already done")
            print("epchain skipped as already done")    
#
print ("*** all done for {} and {}".format(obsid,task))
logging.info ("*** all done for {} and {}".format(obsid,task))
#
evfind = glob.glob(f"{ppsDir}/*PIEVLI*")
#
# now check if more than one event list, select the one with larger (non-zero)
# exposure time
#
if (len(evfind) >= 1):
    lt_max = 0.0
    for xev in evfind:
        hdu = fits.open(xev)
        lt = hdu[1].header['LIVETIME']
        if (lt >= lt_max):
            evlist = xev
            lt_max = lt
        #
    haveEvlist = True
else:
    logging.error("Event list not produced")
    raise FileNotFoundError

#%%
#
# Jump straight to step04 GTI filtering of the standard workflow
# will reuse the merged GTI and region files from the previous run
#
#
# now filter the event lists with the merged GTI
# will copy the merged GTI from the standard processing
gti_file = "{}/../xmmsas_20180620_1732-17.0.0/gti_pn.fits".format(ppsDir)
#gti_file = "{}/../xmmsas_20180620_1732-17.0.0/gti_merged.fits".format(ppsDir)
if (not os.path.isfile(gti_file)):
    logging.error ("GTI file {} not found".format(gti_file))
    raise FileNotFoundError
#
filt_evlist = f"{ppsDir}/pn_evlist_clean_cti49.fits"

expr = "#XMMEA_EP && gti({},TIME) && (PI>150)".format(gti_file)
inst = 'pn'
#
ev_command = 'evselect table={} withfilteredset=Y filteredset={}_evlist_clean_cti49.fits'.format(os.path.basename(evlist),inst) +  \
   ' destruct=Y keepfilteroutput=T expression=\'{}\''.format(expr)
#
print ("Filtering {}".format(inst))
try:
    result = subprocess.run(ev_command, shell=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
    retcode=result.returncode
    if retcode < 0:
        print("evselect was terminated by signal", -retcode, file=sys.stderr)
        logging.warning("evselect was terminated by signal {}".format(-retcode))
    else:
        print("evselect returned", retcode, file=sys.stderr)
        logging.info("evselect returned {}, \n {}".format(retcode,result.stdout.decode()))
except OSError as e:
    print("Execution of evselect failed:", e, file=sys.stderr)
    logging.error("Execution of evselect failed: {}".format(e))
    #
#%%
if (not os.path.isfile(filt_evlist)):
    logging.error("Cannot find cleaned event list")
    raise FileNotFoundError
print ("*** Generating spectra")
spec_chan_max = 20479
expr1 = "(FLAG==0) && (PATTERN<=4)"
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
src_spec = "pn_src_spectrum_cti49.fits"
command = "evselect table={} withspectrumset=yes spectrumset={}".format(filt_evlist,src_spec) +  \
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
#%%
command = "backscale spectrumset={} badpixlocation={}".format(src_spec,filt_evlist)
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
bg_spec = "pn_bkg_spectrum_cti49.fits"
command = "evselect table={} withspectrumset=yes spectrumset={}".format(filt_evlist,bg_spec) +  \
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
command = "backscale spectrumset={} badpixlocation={}".format(bg_spec,filt_evlist)
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
command = "rmfgen spectrumset={} rmfset=pn_cti49.rmf".format(src_spec,inst)
print ("*** RMF generation")
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
command = "arfgen spectrumset={} arfset=pn_cti49.arf withrmfset=yes rmfset=pn_cti49.rmf".format(src_spec) + \
    " badpixlocation={} detmaptype=psf".format(evlist)
print ("*** ARF generation")
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
command = "specgroup spectrumset={} addfilenames=yes rmfset=pn_cti49.rmf".format(src_spec) + \
    " arfset=pn_cti49.arf backgndset={} groupedset=pn_spectrum_cti49_grp0.fits".format(bg_spec)
print ("*** Link bkg, ARF and RMF filenames to the spectrum")
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
command = "specgroup spectrumset={} mincounts={} oversample=3 rmfset=pn_cti49.rmf".format(src_spec,min_cts) + \
    " arfset=pn_cti49.arf backgndset={} groupedset=pn_spectrum_cti49_grp{}.fits".format(bg_spec,min_cts)
print ("*** Group the spectrum in bins of at least {} counts".format(min_cts))
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