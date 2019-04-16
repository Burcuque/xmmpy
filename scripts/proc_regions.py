#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 24 15:02:16 2018

process the region files, to source and background
the input file should have two lines starting with circle:
    one with text={source}
    and another one with no text={} field

@author: ivaltchanov
"""

import os
from shutil import copyfile
import argparse

home = os.path.expanduser('~')

parser = argparse.ArgumentParser(description='Process region files')
parser.add_argument('obsid', type=str,
                    help='The OBSID to process')
parser.add_argument('--wdir', type=str, default=os.getcwd(),
                    help='The working directory')
#
args = parser.parse_args()
#
#%%
#
wdir = os.path.abspath(args.wdir)
obsid = args.obsid
#wdir = home + "/XMM/CAL/PN_SW/ngc4151"
#obsid = "0112310101"

regdir = "{}/{}/xmmsas_20180620_1732-17.0.0/regions".format(wdir,obsid)

for j in ["mos1","mos2","pn"]:
    fx = '{}/src_region_{}.reg'.format(regdir,j)
    if (not os.path.isfile(fx)):
        continue
    temp = '{}/temp.reg'.format(regdir)
    copyfile(fx,temp)
    haveSourceRegion = False
    bg = '{}/bkg_region_{}.reg'.format(regdir,j)
    with open(temp) as tmp:
        lines = tmp.readlines()
    if (len(lines) <= 2):
        print ("{} already processed".format(j))
        continue
    for q in lines:
        if ('circle' in q):
            if ('text' in q):
                #print ("Found source line")
                with open(fx,'w') as f:
                    print (q,file=f)
                haveSourceRegion = True
            else:
                #print ("Found background line")
                with open(bg,'w') as f:
                    print (q,file=f)
        #
    if (not haveSourceRegion):
        print ("No source region file in file {}".format(fx))
        raise FileNotFoundError
    os.remove(temp)
#
