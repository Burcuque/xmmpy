#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  7 17:35:58 2018

@author: ivaltchanov
"""
import os
from astropy.table import Table
#
home = os.path.expanduser('~')
wdir = "{}/XMM/CAL/PN_SW".format(home)
#
# csv file from NXSA
#
# fname = 'pn_sw_agn.csv'
fname = 'NXSA-results-WR140.csv'

t = Table.read('{}/{}'.format(wdir,fname))
#
tt = t['TARGET_TYPE.DESCRIPTION']
nt = len(tt)
ix = []
for i in range(nt):
    if ('TYPE 2' in tt[i]):
        ix.append(i)
#
tx = t[ix]
#
print (tx["OBSERVATION.OBSERVATION_ID","OBSERVATION.REVOLUTION","OBSERVATION.TARGET","TARGET_TYPE.DESCRIPTION"])

tq = Table()
tq["obsid"] = tx["OBSERVATION.OBSERVATION_ID"]
tq["rev"] = tx["OBSERVATION.REVOLUTION"]
tq["target"] = tx["OBSERVATION.TARGET"]
tq["expo"] = tx["OBSERVATION.DURATION"]
