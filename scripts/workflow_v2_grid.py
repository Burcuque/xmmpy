#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 15 11:28:29 2018

workflow for PN analysis of CalClosed data
v2: process with epchain to create calibrated event lists
v2_grid: special version to run on the ESAC grid
v3: select events to create a spectrum
v4: simultaneous fit for the 3 lines in the spectrum

@author: ivaltchanov
"""
import os
import glob
from astropy.io import fits
from astropy.table import Table
import tarfile
import subprocess
import sys
import argparse

from lmfit.models import GaussianModel, PolynomialModel

import matplotlib.pyplot as plt

#
# global folders
#
home = os.path.expanduser('~')
wdir = home + '/XMM/CAL'
outputDir = wdir + '/calclosed'
# here I will save the ODF
dataDir = '/xdata/xcaldata/epic_ex/calclosed'
# this place is where I check available ODFs
calDir = '/xdata/xcaldata/CAL_CLOSED/PN_CAL_CLOSED'

scripts = home + '/Dropbox/Work/XMM/scripts'
defsFile = scripts + '/workflow_defs.py'
exec(open(defsFile).read())
#
#t = Table.read('{}/PN_CalClosed_SW.csv'.format(wdir),format='ascii.csv')
t = Table.read('{}/PN_CalClosed_FF.csv'.format(wdir),format='ascii.csv')
obsid_list = t['observation_id']
#
parser = argparse.ArgumentParser()
parser.add_argument("istart", type=int, help="Start index")
parser.add_argument("iend", type=int, help="End index")
#
args = parser.parse_args()
istart = args.istart
iend = args.iend
#
#%%
#
# step 2, loop over obs_id and run epchain
#
newRun = False
nobs = len(obsid_list)
iend = max(iend,nobs)
nsel = iend - istart + 1
#
for i in range(istart,iend):
    iobs = obsid_list[i]
    print ("Doing OBS_ID %i/%i, %i"%(i+1,iend,nobs))
    tt = epchain_for_obsid(iobs, newRun=False)
print ("all done")
#
