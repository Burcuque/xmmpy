#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 18 16:38:44 2018

call ds9 and plot the images with the source and background regions

@author: ivaltchanov
"""

import os
import subprocess
import sys
#import glob

import argparse

#import numpy as np 

#from astropy.table import Table
#from astropy.io import fits
#from astropy.stats import median_absolute_deviation as mad
#from astropy.stats import mad_std

#import matplotlib as mpl
#mpl.use('Agg')

#import matplotlib.pylab as plt

import logging
#
# global folders
#
os.environ['SAS_CCFPATH'] = '/xdata/ccf/pub'
#%%
#
# get the arguments
parser = argparse.ArgumentParser(description='XMM SAS workflow step 05 check regions')
parser.add_argument('obsid', type=str,
                    help='The OBSID to process')
parser.add_argument('inst', type=str, default='pn',
                    help='Instrument')
parser.add_argument('--wdir', type=str, default=os.getcwd(),
                    help='Where the ODF files are')
#
args = parser.parse_args()
#
#%%
#
wdir = os.path.abspath(args.wdir)
obsid = args.obsid
inst = args.inst

#wdir = "/home/ivaltchanov/XMM/CAL/PN_SW/WR140"
#inst = 'mos2'
#obsid = "0555470701"
#
if (not os.path.isdir(wdir)):
    print ("The ODF folder {} does not exist! Cannot continue".format(wdir))
    raise FileNotFoundError
#    
#
ppsDir = os.path.join(wdir,obsid,'xmmsas_20180620_1732-17.0.0')
if (not os.path.isdir(ppsDir)):
    print ("PPS folder {} does not exist. This means step #01 was not run?".format(ppsDir))
    raise FileNotFoundError

os.chdir(ppsDir)

image_file = "{}_image_2000_10000.fits".format(inst)
if (not os.path.isfile(image_file)):
    print ("No image file found {}".format(image_file))
    raise FileNotFoundError
    

ds9comm = "/home/ivaltchanov/bin/ds9 "

ds9comm += image_file
ds9comm += " -regions load \"regions/*_{}.reg\" -zoom 3 -scale log -scale limits 0 50".format(inst)
ds9comm += " -saveimage {}/images/{}_{}_regions.png -quit".format(wdir,obsid,inst)

try:
    result = subprocess.run([ds9comm], shell=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
    retcode = result.returncode
    if retcode < 0:
        print("ds9 was terminated by signal {} \n {}".format(-retcode,result.stdout.decode()))
    #else:
    #    print("ds9 returned {} \n {}".format(retcode,result.stdout.decode()))
except OSError as e:
    print("Execution of ds9 failed:", e, file=sys.stderr)
