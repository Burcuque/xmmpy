#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 09 12:00:00 2018

Step 00 of the workflow for processing of an XMM observation

1. download from NXSA
2. create a subfolder with name = obsid
3. create a sub-subfolder with name = PN_SW 
4. untar the files in obsid folder

@author: ivaltchanov
"""
import os
import tarfile
import subprocess
import sys
import requests
import logging

import argparse

#
os.environ['SAS_CCFPATH'] = '/xdata/ccf/pub'
#%%
#
# get the arguments
parser = argparse.ArgumentParser(description='XMM SAS workflow step 01, CIFBUILD and ODFINGEST')
parser.add_argument('obsid', type=str,
                    help='The OBSID to process')
parser.add_argument('--wdir', type=str, default=os.getcwd(),
                    help='Where the ODF and the PPS will be stored')
#
args = parser.parse_args()
#
#%%
#
# first download the ODF from XMM archive AIO interface with URL
#
obsid = args.obsid
url = "http://nxsa.esac.esa.int/nxsa-sl/servlet/data-action-aio?obsno=%s&level=ODF"%obsid
#
# create output folder for ODF and CCF
#
wdir = os.path.abspath(args.wdir)
if (not os.path.isdir(wdir)):
    print ("The input working folder {} does not exist. Create it and run the script from within.".format(wdir))
    raise FileNotFoundError
#
outputDir = os.path.join(wdir,obsid)
if (not os.path.isdir(outputDir)):
    os.mkdir(outputDir)
#
logging.basicConfig(level=logging.DEBUG,
                    format='%(asctime)s %(levelname)s %(message)s',
                    filename='{}/pn_sw_step00.log'.format(outputDir),
                    filemode='w')
#
# crate output folder for PPS and logging
#
ppsDir = os.path.join(outputDir,"PN_SW")
if (not os.path.isdir(ppsDir)):
    os.mkdir(ppsDir)
#
tarFile = "{}/{}_ODF.tar".format(outputDir,obsid)
#
if (not os.path.isfile(tarFile)):
    print ("No ODF tar file for OBSID {}, downloading from NXSA".format(obsid))
    logging.info ("No ODF tar file for OBSID {}, downloading from NXSA".format(obsid))
    with requests.get(url) as r:
        r.raise_for_status() # ensure we notice bad responses
        with open(tarFile,"wb") as tmp:
            tmp.write(r.content)
else:
    print ("ODF tar file {} found, will skip downloading it again".format(tarFile))
    logging.info ("ODF tar file {} found, will skip downloading it again".format(tarFile))
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
        if (not os.path.isfile("{}/{}".format(outputDir,member.name))):
            print ("Extracting {}".format(member.name))
            f=tar.extract(member,path=outputDir)
        else:
            print ("File {} already extracted, skipping".format(member.name))
        if ("TAR" in member.name):
            foundTar = True
            print ("Extracting the content of {}".format(member.name))
            try:
                tarcomm = "tar xf {}".format(member.name)
                result = subprocess.run(tarcomm, shell=True,stderr=subprocess.STDOUT,stdout=subprocess.PIPE)
                retcode=result.returncode
                if retcode < 0:
                    print("tar was terminated by signal {} \n {}".format(-retcode,result.stdout.decode()))
                    logging.warning("tar was terminated by signal {}".format(-retcode))
                else:
                    print("tar returned {} \n {}".format(retcode,result.stdout.decode()))
                    logging.info("tar returned {} \n {}".format(retcode,result.stdout.decode()))
            except OSError as e:
                print("Execution of tar failed: {}".format(e))
                logging.error("Execution of tar failed: {}".format(e))
            # remove the TAR file, we already extracted its content
            os.remove("{}/{}".format(outputDir,member.name))
if (not foundTar):
    raise FileNotFoundError
    logging.error("TAR file not found")
    
print ("*** all done for {} and step #00".format(obsid))
logging.info ("*** all done for {} and step #00".format(obsid))
#
