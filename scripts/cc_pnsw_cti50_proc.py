#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 09 12:00:00 2018

Full workflow for EPIC processing an XMM observation with test CCF for a new CTI correction.
Will pick up the new CCFs from the "production" ccf folder in HOME/XMM/CAL/ccfprod folder

Will only process PN SmallWindow mode in CALCLOSED
 
The output will be in a foder called <subfolder> (input parameter)

Will use the new EPN_CTI_0049.CCF and EPN_CTI_0050.CCF

@author: ivaltchanov
"""
import os
import subprocess
import sys
import logging
import glob

import argparse

from astropy.io import fits
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
parser = argparse.ArgumentParser(description='XMM SAS workflow for PN SW in CALCLOSED')
parser.add_argument('obsid', type=str,
                    help='The OBSID to process')
parser.add_argument('subfolder', type=str,
                    help='The subfolder where to save the processed products')
#
args = parser.parse_args()
#
#%%
#
obsid = args.obsid
subfolder = args.subfolder
#
odfdir = f'/xdata/xcaldata/XMM/IVAN/PN_SW/CalClosed/{obsid}'
#
if (not os.path.isdir(odfdir)):
    print ("The ODF folder {} does not exist! Cannot continue".format(odfdir))
    raise FileNotFoundError
os.environ["SAS_ODF"] = odfdir
#
ppsDir = os.path.join(odfdir,subfolder)
#
if (not os.path.isdir(ppsDir)):
    print ("PPS folder {} does not exist. Will create it".format(ppsDir))
    os.mkdir(ppsDir)
#
# Now rerun cifbuild with the new CCf file
#
logging.basicConfig(level=logging.DEBUG,
                    format='%(asctime)s %(levelname)s %(message)s',
                    filename=os.path.join(ppsDir,f"pn_proc_{subfolder}.log"),
                    filemode='w')

#os.environ["SAS_CCFPATH"] = "/ccf/valid:/xdata/xcaldata/XMM/IVAN/PN_SW/ccfdev"
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
#
#
#%%
# Have to run ODFINGEST as I moved the files
#
print ("*** ODFINGEST")
os.chdir(odfdir)
odf_comm = f"odfingest withodfdir=yes odfdir=\"{odfdir}\""
status = run_command(odf_comm)
if (status != 0):
    print (f"Execution of {cifbuild} failed.")
    raise Exception
#
#
#%%
os.chdir(ppsDir)
#
ep_comm = f"epchain odf={odfdir} withphagaincolumn=yes propagatecolumns=all backgroundtres=N"
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
# For CALCLOSED no need to filter for GTI, directly extract a spectrum
#
sbin = 5 # spectral bin size
specFile = f'spectrum_bin{sbin}.fits'

evfile = glob.glob(f"{ppsDir}/*PIEVLI*")
if (len(evfile) == 0):
    print ("No event list file was created. Check the processing, cannot continue!")
    raise FileNotFoundError
#
ev_comm = "evselect " + \
        f"table={evfile[0]} energycolumn='PI' withspectrumset=yes" + \
        " expression='(((DETX,DETY) in BOX(517,1430,2500,2250,0)) && (FLAG==0) " + \
        " && (PATTERN==0) && (PAT_SEQ==0))' " + \
        " xcolumn=DETX ycolumn=DETY" +  \
        " withspecranges=yes specchannelmin=0 specchannelmax=20479" + \
        f" spectrumset={specFile} spectralbinsize={sbin}"
status = run_command(ev_comm)
if (status != 0):
    print (f"Execution of {ev_comm} failed.")
    raise Exception

#%%
# Now preparing for the XSPEC
#
hdu = fits.open(specFile)
hdu[1].header['BACKFILE'] = 'NONE'
hdu.writeto(specFile,overwrite=True)
hdu.close()
#%%
#
# Will use already available canned RMF and ARF
#
rmfset = f'{home}/IVAN/PN_SW/epn_e3_sw20_sY9_v17.0.rmf'
arfset = '{home}/IVAN/PN_SW/cc_ancillary.ds'
outfile = 'pn_spectrum_grp0.fits'
#
comm = f"specgroup spectrumset={specFile} addfilenames=yes rmfset={rmfset} " + \
    f"arfset={arfset} groupedset={outfile}"
print ("*** Link bkg, ARF and RMF filenames to the {obsid} spectrum")
status = run_command(comm)
if (status != 0):
    print (f"Execution of {comm} failed.")
    raise Exception

print ("All done")    