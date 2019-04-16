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
wdir = "{}/XMM/CAL/PN_SW/WR140".format(home)
#
# csv file from NXSA
#
# fname = 'pn_sw_agn.csv'
fname = 'NXSA-results-WR140.csv'

t = Table.read('{}/{}'.format(wdir,fname))
#%%
#
obsid = t['OBSERVATION.OBSERVATION_ID']
#
