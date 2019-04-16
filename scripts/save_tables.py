#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Dump the csv tables with results to Latex

@author: ivaltchanov
"""
import os

import numpy as np

from astropy.table import Table,Column

targets = ["ngc2992","ngc3227","ngc3783","ngc4151","ngc4593",'ngc5548','ngc7213']
redshift = {'ngc4151': 0.003262, 'ngc3227': 0.00386, 'mrk1048': 0.0427, 'ngc3783': 0.009755,\
            'ngc4593': 0.008344, 'ngc5506': 0.00589, 'mcg-5-23-16': 0.008226, 'ngc3516': 0.008816,\
            'ngc5548': 0.01627, 'ngc2992':  0.007296, 'ngc1566': 0.005036, 'iras09149': 0.057150,\
            "iras05078": 0.017879, 'ngc7213': 0.005869}
feK = 6.399 # the weitghed mean of Kalpha_2 at 6.3908 (intensity 50) and Kalpha_1 at 6.40308 (intensity 100)

wdir = "/xdata/xcaldata/XMM/IVAN/PN_SW/sources"
os.chdir(wdir)

#%%
#
# version of CSV tables to use
#
ver = "t28"
#
for target in targets:
    print (f"Process table for {target}")
    out_tab = f'{target}/{target}_output_xspec_cti49_{ver}.csv'
    if (not os.path.isfile(out_tab)):
        print (f"No CTIv49 results found for target {target}: {out_tab}")
        continue
    t = Table.read(out_tab,data_start=0,names=("obsid","rev","delta_time","submode",\
                                  "full","xfilt","inst","ontime",\
                                  "lineE","lineE_low","lineE_up","cstat","chi2r","dof"))
    nt = len(t)
    #
    lineX =  feK/(1.0 + redshift[target]) # redshifted line position
    line = t['lineE'].data
    lineErrLo = t['lineE_low'].data
    lineErrUp = t['lineE_up'].data
    diff = (line - lineX)*1000.0 # in eV
    diff_err = 1000*(lineErrLo + lineErrUp)/2.0 # in eV
    t.add_column(Column(diff),index=11,name='deltaE')
    t.add_column(Column(diff_err),index=12,name='deltaE_err')
    t['deltaE'].format = '.1f'
    t['deltaE_err'].format = '.1f'
    t['delta_time'].format = '.4f'
    tx = t['obsid','rev','delta_time','ontime','lineE','lineE_low','lineE_up',\
           'deltaE','deltaE_err','cstat','chi2r','dof']
    tx.write(f'{wdir}/{target}_xspec_cti49_{ver}.tex',format='ascii.latex',overwrite=True)
#