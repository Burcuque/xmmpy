#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 11 11:09:13 2019

@author: ivaltchanov
"""

from astropy.table import Table

def read_pn_obstable():
    #
    # load the observations from EPICMON account
    #
    
    t = Table.read('/home/epicmon/monitoring/mjss/BadPix/db/pn_all.dat',delimiter='|',\
                   format='ascii.no_header',data_start=1,\
                   names=('rev','ac_obs_id','ac_exp_id','exp_start_time','exp_end_time','perf_dur',\
                          'mode','filter','target','s_obs_id','s_exp_id','ro','rx','rs'))
    return t
#
