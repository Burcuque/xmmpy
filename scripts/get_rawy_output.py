#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 14 12:22:44 2018

Get he RAWY coordinates for source regions

@author: ivaltchanov
"""

import os
import subprocess
#from astropy.io import fits
from astropy.table import Table, Column

def get_rawy(iobs,odfdir,ppsdir):
    #
    # Extracts the RAWY value of s source region for a OBS_ID
    #
    # uses XMM-SAS ecoordconv task which needs ccf.cif and image file 
    # 
    #
    if (not os.path.isdir(odfdir)):
        print ("The ODF folder {} does not exist! Cannot continue".format(odfdir))
        #raise FileNotFoundError
        return None
    #
    #ppsDir = os.path.join(odfdir,"xmmsas_20180620_1732-17.0.0")
    #
    if (not os.path.isdir(ppsdir)):
        print ("PPS folder {} does not exist. Will create it".format(ppsdir))
        #raise FileNotFoundError
        return None
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
    os.chdir(ppsdir)
    #
    src_reg_file = os.path.join(ppsdir,'regions/src_region_pn.reg')
    #
    if (not os.path.isfile(src_reg_file)):
        print (f"Missing PN region file for source: {src_reg_file}")
        #raise FileNotFoundError
        return None
    img_file = "pn_image_2000_10000.fits"
    if (not os.path.isfile(img_file)):
        print (f"Missing PN image file for source: {img_file}")
        #raise FileNotFoundError
        return None
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
    command = f'{task} imageset={img_file} x={detx} y={dety} coordtype=POS'
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
    rawy = None
    for line in result.stdout.decode().split('\n'):
        if ('RAWY' in line):
            print (line)
            rawy = line.split()[3]
            return rawy
#%%
#
t = Table.read('output_xspec.csv')
obsids = t['obsid']
wdir = os.getcwd()
#
xrawy = []
for obsid in obsids:
    odfDir = os.path.join(os.path.abspath(wdir),"{:010}".format(obsid))
    ppsDir = os.path.join(odfDir,"xmmsas_20180620_1732-17.0.0")    
    qrawy = get_rawy(obsid,odfDir,ppsDir)
    xrawy.append(qrawy)
#
tabfile = f'{wdir}/output_xspec_rawy.csv'
t.add_column(Column(xrawy),name='rawy')
t.write(tabfile,overwrite=True)