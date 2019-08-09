#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 09 12:00:00 2018

PN LW processing:  running epchain with testing Long-term CTI (LTCTI) CCF or with no LTCTI
    using source and background regions for AGNs

For no LTCTI set subfolder input to 'nocti', to test the CTI file set a different name, e.g. cti51x

Will only process PN LargeWindow mode.
 
The output will be in a foder called <subfolder>

Using epchain

Modification history:
    v1 (this): will use the ODF from /xdata/xcaldata/XMM/IVAN but will create the products in $TMPDIR

@author: ivaltchanov
"""
import os
import shutil
import subprocess
import sys
import logging
import glob

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

#%%
#
# get the arguments
parser = argparse.ArgumentParser(description='XMM SAS full processing chain for PN LW')
parser.add_argument('obsid', type=str,
                    help='The OBSID to process')
parser.add_argument('expo', type=str,default="S003",
                    help='The exposure to process')
parser.add_argument('subfolder', type=str,default="cti51x",
                    help='The subfolder where to save the products')
parser.add_argument('--wdir', type=str,default=os.getcwd(),
                    help='The working directory (root folder)')
parser.add_argument('-c','--clean',  action='store_true',
                    help='If set will clean the subfolder before processing')

#
args = parser.parse_args()
#
#%%
#
# check if we are on a grid node, then use the /scratch area, encoded in TMPDIR env variable
#
copy_to_node = False
if ('TMPDIR' in os.environ.keys()):
    copy_to_node = True
    tmpdir = os.environ["TMPDIR"]

obsid = args.obsid
wdir = os.path.abspath(args.wdir)
odfdir = os.path.join(wdir,obsid)
nexpo = args.expo
subf = args.subfolder
clean = args.clean
#
#
if (not os.path.isdir(odfdir)):
    print (f"The ODF folder {odfdir} does not exist! Cannot continue.")
    raise FileNotFoundError
#
logging.basicConfig(level=logging.DEBUG,
                    format='%(asctime)s %(levelname)s %(message)s',
                    filename=os.path.join(odfdir,f'pnlw_{subf}_agn_proc.log'),
                    filemode='w')
#
ppsDir = os.path.join(odfdir,subf)
#
if (copy_to_node):
    topdir = os.path.join(tmpdir,obsid)
    if (not os.path.isdir(topdir)):
        print (f"Creating grid node top folder {topdir}")
        logging.info (f"Creating grid node top folder {topdir}")
        os.mkdir(topdir)
    ppsDir = os.path.join(topdir,subf)
    print (f"Will save PPS products to the grid node folder {ppsDir}")
    logging.info (f"Will save PPS products to the grid node folder {ppsDir}")

#
# Now clean the folder if asked to
#
if (os.path.isdir(ppsDir) and clean):
    print (f"Found an already existing folder {ppsDir}, removing the files as clean is {clean}.")
    logging.info (f"Found an already existing folder {ppsDir}, removing the files as clean is {clean}.")
    shutil.rmtree(ppsDir)
#
if (not os.path.isdir(ppsDir)):
    print (f"PPS folder {ppsDir} does not exist. Will create it")
    logging.info (f"PPS folder {ppsDir} does not exist. Will create it")
    os.mkdir(ppsDir)
#
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
if ('nocti' in subf):
    os.environ["SAS_CCFPATH"] = "/ccf/valid"
else:
    # testing ccf in ccfdev
    os.environ["SAS_CCFPATH"] = "/ccf/valid:/xdata/xcaldata/XMM/IVAN/ccfdev"
#
os.chdir(odfdir)
#
cif_file = f"{ppsDir}/ccf_{subf}.cif"
cifbuild = f"cifbuild calindexset={cif_file} withccfpath=yes ccfpath=$SAS_CCFPATH fullpath=yes -w 1"
#
status = run_command(cifbuild)
if (status != 0):
    raise Exception
#
if (not os.path.isfile(cif_file)):
    print (f"No {cif_file} file found in ODF folder {ppsDir}.")
    raise FileNotFoundError
#
os.environ["SAS_CCF"] = cif_file
logging.info(f"Set SAS_ODF to {odfdir}.")
logging.info(f"Set SAS_CCF to {cif_file}.")
#%%
# ======== EPCHAIN =============
#
os.chdir(ppsDir)
#
task = "epchain"
if ('nocti' in subf):
    xtask = f"{task} odf={odfdir} withctilongterm=N odfaccess=all backgroundtres=N "
else:
    xtask = f"{task} odf={odfdir} withctilongterm=Y odfaccess=all backgroundtres=N "
#
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
# now filter the event lists with the merged GTI
print (f"Filtering the calibrated event lists with the GTI for expo {nexpo}")
gti_file = f"{odfdir}/PN_LW/gti_pn_{nexpo}.fits"
if (not os.path.isfile(gti_file)):
    logging.error ("GTI file {} not found".format(gti_file))
    raise FileNotFoundError
#
# check the duration of the GTI
#
filt_evlist = f"{ppsDir}/pn_evlist_clean_{nexpo}.fits"

expr = f"#XMMEA_EP && gti({gti_file},TIME) && (PI>150)"
#
ev_command = f'evselect table={evlist} withfilteredset=Y filteredset={filt_evlist}' +  \
   f' destruct=Y keepfilteroutput=T expression=\'{expr}\''
status = run_command(ev_command)
if status != 0:
    raise Exception
#
#%%
#
# Now generate the spectra, will reuse the regions from previous run
#
evlist = os.path.join(ppsDir,f'pn_evlist_clean_{nexpo}.fits')
if (not os.path.isfile(evlist)):
    print(f"GTI filtered event list {evlist} not available. Have you run pn_lw_step03.py ?")
    logging.error(f"GTI filtered event list {evlist} not available.")
    raise FileNotFoundError
#
print ("*** Generating spectra")
#
reg_file = f"{wdir}/regions/{obsid}_pn_regions.reg"
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
# now copy back the files from the node
#
if (copy_to_node):
    dest_dir = os.path.join(odfdir,subf)
    # check if the destination folder exist, then remove it
    if (os.path.isdir(dest_dir)):
        print (f"Found an already existing folder {dest_dir}. Will remove it!")
        logging.warning (f"Found an already existing folder {dest_dir}. Will remove it!")
        shutil.rmtree(dest_dir)
    shutil.copytree(ppsDir,dest_dir)
    print (f"Folder {ppsDir} copied to {dest_dir}")
    logging.info (f"Folder {ppsDir} copied to {dest_dir}")
#
logging.info (f"pnlw_{subf}_agn_proc done")
#
print ("All done")# will attache a canned response file for PN in LW, single events
