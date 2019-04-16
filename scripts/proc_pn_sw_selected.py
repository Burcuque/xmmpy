#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 14 15:35:12 2018

combine all no_cti results in one table, the target name will be the first column

@author: ivaltchanov
"""
import os
import numpy as np
from astropy.table import Table, Column, vstack

targets = ["ngc3783","ngc4151","ngc4593","ngc2992","ngc3227","ngc5548"]
#

wdir = "/home/ivaltchanov/IVAN/PN_SW"
#
stack = False
#
for target in targets:
    tfile = f"{wdir}/{target}/output_xspec_nocti.csv"
    if (not os.path.isfile(tfile)):
        print (f"No CTI table for {target} not found")
        continue
    t = Table.read(tfile)
    t['obsid'].format = "010d"
    t['delta_time'].format = "5.2f"
    t['chi2r'].format = "6.3f"
    t['rawy'].format = "5.1f"
    name_col = np.chararray(len(t),itemsize=7)
    name_col[:] = target
    t.add_column(Column(name_col),name='target',index=0)
    if (not stack):
        tout = t
        stack = True
    else:
        tout = vstack([tout,t])
#
tout.write(f"{wdir}/all_targets_nocti.csv",overwrite=True)