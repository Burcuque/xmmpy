#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 09 12:00:00 2018

EPIC workflow for processing an XMM observation

@author: ivaltchanov
"""
import os
import tarfile
import subprocess
import sys
import requests
#import ipyparallel as ipp

#import argparse

#import numpy as np 

#from astropy.io import fits
#
# global folders
#
home = os.path.expanduser('~')
wdir = home + '/XMM/CAL'
#outputDir = wdir + '/CLS'
#dataDir = '/xdata/xcaldata/XMM'
#xcalDir = '/xdata/xcaldata/CAL_CLOSED/PN_CAL_CLOSED'
# set the SAS version XMMSAS 17.0
sasversion = 'xmmsas_20180504_1731'
os.environ['SAS_CCFPATH'] = '/xdata/ccf/pub'
#%%
#
# get the arguments
#parser = argparse.ArgumentParser(description='XMM SAS pipeline')
#parser.add_argument('obsid', type=str,default='0111240101',
#                    help='The OBSID to process')
#parser.add_argument('cifDo', type=bool, default=False,
#                    help='Create a new CCF file')
#parser.add_argument('odfDo', type=bool, default=False,
#                    help='Do the ODF ingest')
#parser.add_argument('task', type=str,default='emproc',
#                    help='Which proc to run: emproc or epproc')
#
#args = parser.parse_args()
#
#%%
#
# first download the ODF from XMM archive AIO interface with URL
#
# Circinus galaxy
#obsid = "0111240101"
obsid = "0656580601"
#obsid = args.obsid
url = "http://nxsa.esac.esa.int/nxsa-sl/servlet/data-action-aio?obsno=%s&level=ODF"%obsid
#
# create output folder
#
outputDir = "{}/{}".format(wdir,obsid)
if (not os.path.isdir(outputDir)):
    os.mkdir(outputDir)
#
tarFile = "{}/{}_ODF.tar".format(outputDir,obsid)
#
if (not os.path.isfile(tarFile)):
    print ("No ODF tar file for OBSID {}, downloading from NXSA".format(obsid))
    with requests.get(url) as r:
        r.raise_for_status() # ensure we notice bad responses
        with open(tarFile,"wb") as tmp:
            tmp.write(r.content)
else:
    print ("ODF tar file {} found, will skip downloading it again".format(tarFile))
#
# now untar downloaded tar file and untar the TAR file inside
#
# will do it member by member to avoid extracting files that already exist
#
foundTar = False
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
            with tarfile.open("{}/{}".format(outputDir,member.name)) as tx:
                for xmember in tx.getmembers():
                    # check member is alredy in folder
                    if (not os.path.isfile("{}/{}".format(outputDir,xmember.name))):
                        print ("Extracting {}".format(xmember.name))
                        f=tar.extract(xmember,path=outputDir)
                    else:
                        print ("File {} already extracted, skipping".format(xmember.name))
            #
            # remove the TAR file, we already extracted its content
            os.remove("{}/{}".format(outputDir,member.name))
if (not foundTar):
    raise FileNotFoundError
    
#
#%%
#
# now the SAS processing, must initialise HEADAS and XMMSAS before running the code!
#
odfdir = outputDir
os.chdir(odfdir)
os.environ["SAS_ODF"] = odfdir
ccffile = "{}/ccf.cif".format(odfdir)
#
#%%
# cifbuild
#cifDo = args.cifDo
cifDo = False
# now check if CCF file exists if cifDo = False
if ((not cifDo) and (not os.path.isfile(ccffile))):
    cifDo = True
    print ("Will process cifbuild although you said to skip it <= there is no CCF file in folder")
#
if (cifDo):
    print ("*** CIFBUILD")
    try:
        result = subprocess.run("cifbuild", shell=True)
        retcode=result.returncode
        if retcode < 0:
            print("cifbuild was terminated by signal", -retcode, file=sys.stderr)
        else:
            print("cifbuild returned", retcode, file=sys.stderr)
    except OSError as e:
        print("Execution of cifbuild failed:", e, file=sys.stderr)
    # good run of cifbuild
    if (retcode == 0):
        if (not os.path.isfile(ccffile)):
            print ("cifbuild executed OK but no ccf.cif file is found in %s"%odfdir)    
#
# add the CCF file to environment
#
os.environ["SAS_CCF"] = ccffile
#%%
# now ODF ingest
#odfDo = args.odfDo
odfDo = False
if (odfDo):
    print ("*** ODFINGEST")
    try:
        result = subprocess.run("odfingest", shell=True)
        retcode = result.returncode
        if retcode < 0:
            print("odfingest was terminated by signal", -retcode, file=sys.stderr)
        else:
            print("odfingest returned", retcode, file=sys.stderr)
    except OSError as e:
        print("Execution of odfingest failed:", e, file=sys.stderr)
#%%
def doProcess(task="emproc",folder=""):
    #
    # The emproc/epproc processing
    #
    os.chdir(folder)
    logfile = open("%s/%s_processing.log"%(folder,task),'wb')    
    try:
        result = subprocess.run(task, shell=False,stderr=subprocess.STDOUT,stdout=subprocess.PIPE)
        retcode = result.returncode
        logfile.write(result.stdout)
        if retcode < 0:
            print("%s was terminated by signal"%task, -retcode, file=sys.stderr)
        else:
            print("%s returned"%task, retcode, file=sys.stderr)
    except OSError as e:
        print("Execution of %s failed:"%task, e, file=sys.stderr)
    logfile.close()
    return retcode
#
#%%
# create a folder to save the processed files
procDir = "%s/%s"%(odfdir,sasversion)
if (not os.path.isdir(procDir)):
    print ("Creating folder {}".format(procDir))
    os.mkdir(procDir)    
#c = ipp.Client()
#%%
#task = 'emproc'
task = 'emproc'
print ("*** RUNNING %s"%task)
doProcess(task=task,folder=procDir)
print ("*** all done for {} and {}".format(obsid,task))
#%%
task = 'epproc'
print ("*** RUNNING %s"%task)
doProcess(task=task,folder=procDir)
print ("*** all done for {} and {}".format(obsid,task))
#
