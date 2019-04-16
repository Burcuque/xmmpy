#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 16 14:06:43 2019

@author: ivaltchanov
"""
from astropy.table import Table
import numpy as np

ccdir = '/xdata/xcaldata/XMM/IVAN/PN_SW/CalClosed'
#tcc = Table.read(f"{ccdir}/pn_sw_calclosed_fit_results_nocti_sorted.csv",comment="\s*#")
tcc = Table.read(f"{ccdir}/cc_output_xspec_no_cti.csv",comment="\s*#")
tcc.sort('time')
rtime = tcc['time']
ix = np.where(rtime >= 1.0)[0]
for j in ix:
    print ('"{:010}" '.format(tcc['obsid'].data[j]),end='')

