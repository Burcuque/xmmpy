#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 09 12:00:00 2018

Step 01 of the workflow for processing of an XMM observation

1. download from NXSA
2. create a subfolder with name = obsid
3. create a sub-subfolder with name = sasversion 
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
                    filename='{}/proc_step01.log'.format(outputDir),
                    filemode='w')
#
# crate output folder for PPS and logging
#
# get the SAS version, will be used for the output folder 
try:
    result = subprocess.run(["sasversion -v"], shell=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
    retcode = result.returncode
    if retcode < 0:
        print("sasversion was terminated by signal {} \n {}".format(-retcode,result.stdout.decode()))
        logging.warning("sasversion was terminated by signal {} \n {}".format(-retcode,result.stdout.decode()))
    else:
        logging.info("sasversion returned {} \n {}".format(retcode,result.stdout.decode()))
except OSError as e:
    print("Execution of sasversion failed:", e, file=sys.stderr)
    logging.error("Execution of sasversion failed: {}".format(e))
#
sasversion = result.stdout.decode().split('[')[1].split(']')[0]
#
ppsDir = os.path.join(outputDir,sasversion)
if (not os.path.isdir(ppsDir)):
    os.mkdir(ppsDir)
#
tarFile = "{}/{}_ODF.tar".format(outputDir,obsid)
#
if (not os.path.isfile(tarFile)):
    print ("No ODF tar file for OBSID {}, downloading from NXSA".format(obsid))
    logging.info ("No ODF tar file for OBSID {}, downloading from NXSA".format(obsid))
    with requests.post(url, auth=('ivaltch', 'XXXX')) as r:
    #with requests.get(url) as r:
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
    
#
#%%
#
# now the SAS processing, must initialise HEADAS and XMMSAS before running the code!
#
odfdir = outputDir
os.chdir(odfdir)
os.environ["SAS_ODF"] = odfdir
ccffile = "{}/ccf.cif".format(odfdir)
logging.info("Set SAS_ODF to {}".format(odfdir))
#
#%%
# cifbuild
print ("*** CIFBUILD")
try:
    result = subprocess.run("cifbuild", shell=True,stderr=subprocess.STDOUT,stdout=subprocess.PIPE)
    retcode=result.returncode
    if retcode < 0:
        print("cifbuild was terminated by signal {} \n {}".format(-retcode))
        logging.warning("cifbuild was terminated by signal {}".format(-retcode))
    else:
        print("cifbuild returned {}".format(retcode))
        logging.info("cifbuild returned {} \n {}".format(retcode,result.stdout.decode()))
except OSError as e:
    print("Execution of cifbuild failed: {}".format(e))
    logging.error("Execution of cifbuild failed: {}".format(e))
#
# add the CCF file to environment
#
os.environ["SAS_CCF"] = ccffile
logging.info("Set SAS_CCF to {}".format(ccffile))

#%%
print ("*** ODFINGEST")
try:
    result = subprocess.run("odfingest", shell=True,stderr=subprocess.STDOUT,stdout=subprocess.PIPE)
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
print ("*** all done for {} and step #01".format(obsid))
logging.info ("*** all done for {} and step #01".format(obsid))
#
