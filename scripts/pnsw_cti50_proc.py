#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 09 12:00:00 2018

Full workflow for EPIC-PN processing an XMM observation.

Wil only process PN SmallWindow mode.

Will use backgroundtres=N in epchain to speed up the processing

Will use the merged GTI and region files from the standard procesing

Will use the new EPN_CTI_0049.CCF for the long-term CTI correction

The output will be in a foder called <subfolder>

@author: ivaltchanov
"""
import os
import subprocess
import sys
import logging
import glob

import argparse

home = os.path.expanduser('~')

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
parser = argparse.ArgumentParser(description='XMM SAS workflow for EPN in SW with new CCF')
parser.add_argument('obsid', type=str,
                    help='The OBSID to process')
parser.add_argument('subfolder', type=str,
                    help='The subfolder where to save the products')
#
args = parser.parse_args()
#
#%%
#
obsid = args.obsid
cwd = os.getcwd()
odfdir = os.path.join(f'{home}/IVAN/PN_SW/sources',cwd,obsid)
subfolder = args.subfolder
#
if (not os.path.isdir(odfdir)):
    print ("The ODF folder {} does not exist! Cannot continue".format(odfdir))
    raise FileNotFoundError
#
ppsDir = os.path.join(odfdir,subfolder)
#
if (not os.path.isdir(ppsDir)):
    print (f"PPS folder {ppsDir} does not exist. Will create it")
    newRun = True
    os.mkdir(ppsDir)
#
logging.basicConfig(level=logging.DEBUG,
                    format='%(asctime)s %(levelname)s %(message)s',
                    filename=os.path.join(ppsDir,f'pn_proc_{subfolder}.log'),
                    filemode='w')
#
#
os.environ["SAS_ODF"] = odfdir
#
# add the testing fiolder at the end of the CCFPATH
#
os.environ["SAS_CCFPATH"] = f"/ccf/valid:{home}/XMM/CAL/ccfprod"
#
os.chdir(odfdir)
#
cif_file = f"ccf_test_{subfolder}.cif"
#
cifbuild = f"cifbuild calindexset={cif_file} withccfpath=yes ccfpath=$SAS_CCFPATH fullpath=yes -w 1"
#
status = run_command(cifbuild)
if (status != 0):
    print (f"Execution of {cifbuild} failed.")
    raise Exception
#
ccffile = os.path.join(odfdir,cif_file)
if (not os.path.isfile(ccffile)):
    print (f"No CIF file found in ODF folder {odfdir}. Run ccfbuild")
    raise FileNotFoundError
#
os.environ["SAS_CCF"] = ccffile
logging.info(f"Set SAS_ODF to {odfdir}")
logging.info(f"Set SAS_CCF to {ccffile}")
#%%
# Have to run ODFINGEST as I moved the files
#
print ("*** ODFINGEST")
os.chdir(odfdir)
#
odf_comm = f"odfingest withodfdir=yes odfdir=\"{odfdir}\""
status = run_command(odf_comm)
if (status != 0):
    print (f"Execution of {cifbuild} failed.")
    raise Exception
#%%
os.chdir(ppsDir)
#
ep_comm = f"epchain odf={odfdir} odfaccess=all backgroundtres=N"
print (f"*** RUNNING {ep_comm}")
status = run_command(ep_comm)
if (status != 0):
    print (f"Execution of {ep_comm} failed.")
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
# Jump straight to step04 GTI filtering of the standard workflow
# will reuse the merged GTI and region files from the previous run
#
#
# now filter the event lists with the merged GTI
# will copy the merged GTI from the standard processing
gti_file = f"../xmmsas_20180620_1732-17.0.0/gti_pn.fits"
if (not os.path.isfile(gti_file)):
    logging.error (f"GTI file {gti_file} not found")
    raise FileNotFoundError
#
filt_evlist = f"pn_evlist_clean_{subfolder}.fits"

expr = f"#XMMEA_EP && gti({gti_file},TIME) && (PI>150)"
#
ev_comm = f'evselect table={evlist} withfilteredset=Y filteredset=pn_evlist_clean_{subfolder}.fits' +  \
   f' destruct=Y keepfilteroutput=T expression=\'{expr}\''
#
print ("Filtering pn")
status = run_command(ev_comm)
if (status != 0):
    print (f"Execution of {ev_comm} failed.")
    raise Exception
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
src_reg_file = '../xmmsas_20180620_1732-17.0.0/regions/src_region_pn.reg'
bkg_reg_file = '../xmmsas_20180620_1732-17.0.0/regions/bkg_region_pn.reg'
#
if (not (os.path.isfile(src_reg_file) and os.path.isfile(bkg_reg_file))):
    print ("Missing region file for source or background")
    print (f"Please create two files with names {src_reg_file} (source region) and {bkg_reg_file} (background region) ")
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
src_spec = f"pn_src_spectrum_{subfolder}.fits"
comm = f"evselect table={filt_evlist} withspectrumset=yes spectrumset={src_spec}" +  \
    f" energycolumn=PI spectralbinsize=5 withspecranges=yes specchannelmin=0 specchannelmax={spec_chan_max}" +  \
    f" expression='{expr1} && ((X,Y) IN {src_reg})'"
print ("*** Extracting spectrum in the source region")
status = run_command(comm)
if (status != 0):
    print (f"Execution of {comm} failed.")
    raise Exception
#
# now calculate the backscale for the source region
#%%
comm = f"backscale spectrumset={src_spec} badpixlocation={filt_evlist}"
print ("*** Backscale the source spectrum")
status = run_command(comm)
if (status != 0):
    print (f"Execution of {comm} failed.")
    raise Exception
#%%
# extract the background spectrum
#        
with open(bkg_reg_file) as regfile:
    reg_line = regfile.readline()
bg_reg = reg_line.strip().split()[0]
#
bg_spec = f"pn_bkg_spectrum_{subfolder}.fits"
comm = f"evselect table={filt_evlist} withspectrumset=yes spectrumset={bg_spec}" +  \
    f" energycolumn=PI spectralbinsize=5 withspecranges=yes specchannelmin=0 specchannelmax={spec_chan_max}" +  \
    f" expression='{expr1} && ((X,Y) IN {bg_reg})'"
print ("*** Extracting spectrum in the background region")
status = run_command(comm)
if (status != 0):
    print (f"Execution of {comm} failed.")
    raise Exception
#
# now the backscale for the background
#
comm = f"backscale spectrumset={bg_spec} badpixlocation={filt_evlist}"
print ("*** Backscale the background spectrum")
status = run_command(comm)
if (status != 0):
    print (f"Execution of {comm} failed.")
    raise Exception
#%%
# now generate the RMF (can take time)
rmfset = f"pn_{subfolder}.rmf"
comm = f"rmfgen spectrumset={src_spec} rmfset={rmfset}"
print ("*** RMF generation")
status = run_command(comm)
if (status != 0):
    print (f"Execution of {comm} failed.")
    raise Exception
#%%
# now generate the ARF
arfset = f"pn_{subfolder}.arf"
comm = f"arfgen spectrumset={src_spec} arfset={arfset} withrmfset=yes rmfset={rmfset}" + \
    f" badpixlocation={evlist} detmaptype=psf"
print ("*** ARF generation")
status = run_command(comm)
if (status != 0):
    print (f"Execution of {comm} failed.")
    raise Exception
#
#%%
#
# Use specgroup to add background, ARF and RMF files in headers. Without any grouping first
#
comm = f"specgroup spectrumset={src_spec} addfilenames=yes rmfset={rmfset}" + \
    f" arfset={arfset} backgndset={bg_spec} groupedset=pn_spectrum_{subfolder}_grp0.fits"
print ("*** Link bkg, ARF and RMF filenames to the spectrum")
status = run_command(comm)
if (status != 0):
    print (f"Execution of {comm} failed.")
    raise Exception

print ("All done")