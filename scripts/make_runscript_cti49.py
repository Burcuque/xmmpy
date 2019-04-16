#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 15 15:10:55 2018

Generate a run script for the grid for the selected AGN targets to process on grid

@author: ivaltchanov
"""
import os
from astropy.table import Table

wdir = '/xdata/xcaldata/XMM/IVAN/PN_SW/sources'
scripts = '/xdata/xcaldata/XMM/IVAN/scripts'
python = "/home/ivaltchanov/miniconda3/bin/python"
#
t = Table.read(f"{wdir}/selected_targets_nocti.csv",comment="\s*#")
targets = t['target']
obsids = t['obsid']

fout1 = open("/xdata/xcaldata/XMM/IVAN/scripts/grid_run_selected_t8.sh","w")
fout2 = open("/xdata/xcaldata/XMM/IVAN/scripts/cti49_fit_selected_t8.sh","w")

print ("#!/bin/sh",file=fout1)
print ("#",file=fout1)
print ("#!/bin/sh",file=fout2)
print ("#",file=fout2)

for target,iobs in zip(targets,obsids):
    xobs = "{:010}".format(iobs)
    xpath = os.path.join(wdir,target)
    print (f"cd {xpath}",file=fout1)
    odfidr = os.path.join(xpath,xobs)
    print (f"qsub -N {target}_{iobs} -v obsid=\"{xobs}\" /xdata/xcaldata/XMM/IVAN/scripts/cti49_proc_on_grid.sh",file=fout1)
    #
    print (f"cd {xpath}",file=fout2)
    print (f"{python} -u {scripts}/plot_and_fit_spectra_cti49.py \"{xobs}\" \"pn\" --output {target}_output_xspec_cti49.csv",file=fout2)
#
fout1.close()
fout2.close()


    