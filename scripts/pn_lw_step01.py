#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 09 12:00:00 2018

Step 01 of the workflow for processing of an XMM observation

1. download from NXSA
2. create a subfolder with name = obsid
4. untar the files in obsid folder
5. sets the ODF_DIR
6. runs cifbuild
7. sets the ccf.cif 
8. run odfingest

version 1

@author: ivaltchanov
"""
import os
import tarfile
import subprocess
import sys
import requests
import logging

import argparse

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
os.environ['SAS_CCFPATH'] = '/xdata/ccf/pub'
#%%
#
# get the arguments
parser = argparse.ArgumentParser(description='XMM SAS workflow step 01, CIFBUILD and ODFINGEST')
parser.add_argument('obsid', type=str,
                    help='The OBSID to process')
parser.add_argument('--wdir', type=str, default=os.getcwd(),
                    help='Where the ODF will be stored')
#
args = parser.parse_args()
#
#%%
#
# first download the ODF from XMM archive AIO interface with URL
#
obsid = args.obsid
url = f"http://nxsa.esac.esa.int/nxsa-sl/servlet/data-action-aio?obsno={obsid}&level=ODF"
#
# create output folder for ODF and CCF
#
wdir = os.path.abspath(args.wdir)
if (not os.path.isdir(wdir)):
    print (f"The input working folder {wdir} does not exist. Create it and run the script from within.")
    raise FileNotFoundError
#
outputDir = os.path.join(wdir,obsid)
if (not os.path.isdir(outputDir)):
    os.mkdir(outputDir)
#
logging.basicConfig(level=logging.DEBUG,
                    format='%(asctime)s %(levelname)s %(message)s',
                    filename=f'{outputDir}/pn_lw_step01.log',
                    filemode='w')
#
# crate output folder for PPS and logging
#
tarFile = f"{outputDir}/{obsid}_ODF.tar"
#
if (not os.path.isfile(tarFile)):
    print (f"No ODF tar file for OBSID {obsid}, downloading from NXSA")
    logging.info (f"No ODF tar file for OBSID {obsid}, downloading from NXSA")
    with requests.get(url) as r:
        r.raise_for_status() # ensure we notice bad responses
        with open(tarFile,"wb") as tmp:
            tmp.write(r.content)
else:
    print (f"ODF tar file {tarFile} found, will skip downloading it again")
    logging.info (f"ODF tar file {tarFile} found, will skip downloading it again")
#
# now untar downloaded tar file and untar the TAR file inside
#
# will do it member by member to avoid extracting files that already exist
#
foundTar = False
os.chdir(outputDir)

with tarfile.open(tarFile,'r') as tar:
    for member in tar.getmembers():
        # check member is alredy in folder
        if (not os.path.isfile(f"{outputDir}/{member.name}")):
            print (f"Extracting {member.name}")
            f=tar.extract(member,path=outputDir)
        else:
            print (f"File {member.name} already extracted, skipping")
        if ("TAR" in member.name):
            foundTar = True
            print (f"Extracting the content of {member.name}")
            tarcomm = f"tar xf {member.name}"
            status = run_command(tarcomm)
            if (status != 0):
                print (f"Failed to extract {member.name}. Cannot continue!")
                logging.error (f"Failed to extract {member.name}. Cannot continue!")
                raise Exception
            else:
                # remove the TAR file, we already extracted its content
                os.remove(f"{outputDir}/{member.name}")
if (not foundTar):
    raise FileNotFoundError
    print(f"TAR file for {obsid} not found")
    logging.error(f"TAR file for {obsid} not found")
#
#%%
#
# now the SAS processing, must initialise HEADAS and XMMSAS before running the code!
#
odfdir = outputDir
os.environ["SAS_ODF"] = odfdir
ccffile = f"{odfdir}/ccf.cif"
logging.info(f"Set SAS_ODF to {odfdir}")
#
#%%
# cifbuild
print ("*** CIFBUILD")
comm = "cifbuild"
status = run_command(comm)
if (status != 0):
    print("Execution of cifbuild failed. Cannot continue.")
    logging.error("Execution of cifbuild failed. Cannot continue.")
    raise Exception
#
# add the CCF file to environment
#
os.environ["SAS_CCF"] = ccffile
logging.info(f"Set SAS_CCF to {ccffile}")

#%%
print ("*** ODFINGEST")
comm = "odfingest"
status = run_command(comm)
if (status != 0):
    print("Execution of odfingest failed. Cannot continue.")
    logging.error("Execution of odfingest failed. Cannot continue.")
    raise Exception

print (f"*** all done for {obsid} and step #01")
logging.info (f"*** all done for {obsid} and step #01")
#
