#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 20 10:21:48 2018

@author: ivaltchanov
"""
from astropy.table import Table
import os
import shutil

# read the table with obs_id
t = Table.read('/xdata/xcaldata/XMM/IVAN/PN_calclosed/PN_CalClosed_SW.csv',format='ascii.csv')
obsid_list = t['observation_id']
#
# now copy the ODF folders to outDir
#
outDir = '/xdata/xcaldata/XMM/IVAN/PN_calclosed'
oldDir = '/xdata/xcaldata/CAL_CLOSED/PN_CAL_CLOSED'
#
for iobs in obsid_list:
    inDir = f"{oldDir}/{iobs}"
    xobsid = "%010i"%int(iobs)
    xdir = f"{outDir}/{xobsid}"
    print (f"Copying files from  {inDir} to {xdir}")
    if (not os.path.isdir(xdir)):
        os.mkdir(xdir)
    #
    names = os.listdir(inDir)
    for name in names:
        src = os.path.join(inDir,name)
        shutil.copy(src,xdir)
    #
#