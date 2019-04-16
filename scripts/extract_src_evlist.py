#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 23 10:50:09 2018

Extract the event list in the source area

At the end of each extraction it will calculate the average RAWY

@author: ivaltchanov
"""

import os
import sys
import subprocess
import numpy as np

from astropy.table import Table, Column
from astropy.io import fits

#import matplotlib.pylab as plt
#from matplotlib.patches import Patch
targets = ["ngc4151","ngc3227",'ngc3783',"ngc5506","ngc5548", "ngc4593","ngc2992"]
#
wdir = "/xdata/xcaldata/XMM/IVAN/PN_SW"
os.chdir(wdir)

fout = open(f"{wdir}/rawy_stats.csv","w")
print ("target,obsid,rev,time,rawy",file=fout)
fout2 = open(f"{wdir}/rawy_target_stats.csv","w")
print ("target,rawy,rawy_err",file=fout2)

#%%
for target in targets:
    # the celaned event list will be the input
    tabfile = f'{wdir}/{target}/output_xspec_cti49.csv'
    if (os.path.isfile(tabfile)):
        t = Table.read(tabfile)
        nt = len(t)
    else:
        print (f"{tabfile} not found.")
        raise FileNotFoundError        
    # get the FITS header keyword SRCPOSY
    #
    obsids = t['obsid']
    qrawy = []
    for j,obsid in enumerate(obsids):
        revj = t['rev'].data[j]
        timej = t['delta_time'].data[j]
        ppsDir = f'{target}/{obsid:010}'
        qDir = f'{wdir}/{target}/{obsid:010}/cti49'
        os.chdir(qDir)
        # now all files are relative to qDir
        evlist_c = "pn_evlist_clean_cti49.fits"
        regfile = "../xmmsas_20180620_1732-17.0.0/regions/src_region_pn.reg"
        if (os.path.isfile(evlist_c)):
            if os.path.isfile(regfile):
                with open(regfile) as regfile:
                    reg_line = regfile.readline()
                src_reg = reg_line.strip().split()[0]
                #
                src_evlist = "pn_src_evlist_cti49.fits"
                command = "evselect table=pn_evlist_clean_cti49.fits withfilteredset=yes filteredset={}".format(src_evlist) +  \
                " expression='(FLAG==0) && (PATTERN<=4) && ((X,Y) IN {})'".format(src_reg)
                print ("*** Extracting events in the source region")
                print (command)
                try:
                    result = subprocess.run(command, shell=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
                    retcode=result.returncode
                    if retcode != 0:
                        print("evselect was terminated by signal", -retcode, file=sys.stderr)
                        raise Exception
                    else:
                        print("evselect returned", retcode, file=sys.stderr)
                except OSError as e:
                    print("Execution of evselect failed:", e, file=sys.stderr)
            else:
                print (f"{regfile} not found.")
                raise FileNotFoundError
        else:
            print (f"{evlist_c} not found.")
            raise FileNotFoundError
        hdu = fits.open(src_evlist)
        energy = hdu['EVENTS'].data['PI']
        ix = np.where((energy >= 6000)*(energy <= 7000))[0]
        print ("Will calculate RAWY of {} selected events (out of {})".format(len(ix),len(energy)))
        rawy = np.mean(hdu['EVENTS'].data['RAWY'][ix])
        qrawy.append(rawy)
        #print ("Mean RAWY for target {} and OBS_ID {}: {}".format(target,obsid,rawy))
        print (f"{target},{obsid:010},{revj},{timej:.5f},{rawy:.2f}")
        print (f"{target},{obsid:010},{revj},{timej:.5f},{rawy:.2f}",file=fout)
        #print (f"{target},{obsid},{revj},{timej},{rawy}",file=fout)
    #
    print ("*** Average RAWY for all observations of {}: {} +/- {}".format(target,np.mean(qrawy),np.std(qrawy)))
    print ("*** Average RAWY for all observations of {}: {} +/- {}".format(target,np.mean(qrawy),np.std(qrawy)),file=fout2)
    
#
fout.close()
fout2.close()