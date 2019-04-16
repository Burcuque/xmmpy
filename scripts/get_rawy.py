#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 14 12:22:44 2018

Get he RAWY coordinates for source regions

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
parser = argparse.ArgumentParser(description='XMM SAS coordinate conversion')
parser.add_argument('obsid', type=str,
                    help='The OBSID to process')
parser.add_argument('--odfdir', type=str, default=os.getcwd(),
                    help='Where the ODF files are')
#
args = parser.parse_args()
#
#%%
#
obsid = args.obsid
odfdir = os.path.join(os.path.abspath(args.odfdir),obsid)

if (not os.path.isdir(odfdir)):
    print ("The ODF folder {} does not exist! Cannot continue".format(odfdir))
    raise FileNotFoundError
#
ppsDir = os.path.join(odfdir,"xmmsas_20180620_1732-17.0.0")
#
if (not os.path.isdir(ppsDir)):
    print ("PPS folder {} does not exist. Will create it".format(ppsDir))
    raise FileNotFoundError
#
task = 'ecoordconv'
#
os.environ["SAS_ODF"] = odfdir
#
os.environ["SAS_CCF"] = os.path.join(odfdir,"ccf.cif")
#
#%%
# read the region file for the source
#
#
os.chdir(ppsDir)
#
src_reg_file = os.path.join(ppsDir,'regions/src_region_pn.reg')
#
if (not os.path.isfile(src_reg_file)):
    print ("Missing region file for source")
    raise FileNotFoundError
#
# extract the source spectrum from the region
with open(src_reg_file) as regfile:
    reg_line = regfile.readline()
src_reg = reg_line.strip().split()[0]
#
# now extract X and Y
#
qqq = src_reg.split(',')
detx = qqq[0].replace('circle(','')
dety = qqq[1]

#command = f'{task} imageset=pn_image_2000_10000.fits withcoords=no srcexp=\'(DETX,DETY) in {src_reg}\' coordtype=DET'
command = f'{task} imageset=pn_image_2000_10000.fits x={detx} y={dety} coordtype=POS'

try:
    result = subprocess.run(command, shell=True,stderr=subprocess.STDOUT,stdout=subprocess.PIPE)
    retcode = result.returncode
    if retcode < 0:
        print("{} was terminated by signal {}".format(task,-retcode))
    else:
        print("{} returned {} \n {}".format(task,retcode,result.stdout.decode()))
except OSError as e:
    print("Execution of {} failed: {}".format(task,e))
#
print ("*** RESULT ***")
for line in result.stdout.decode().split('\n'):
    if ('RAWY' in line):
        print (line)
        rawy = line.split()[3]
        print (f"rawy is {rawy}")