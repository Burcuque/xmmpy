#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 11 15:39:47 2018

TAP query the XMM archive (XSA)

@author: ivaltchanov
"""
import io
import os
import time

import requests
import numpy as np
from astropy.io.votable import parse_single_table
from astropy.table import Table

today = time.strftime("%Y-%m-%d %H:%M")

home = os.path.expanduser('~')
outDir = home + '/IVAN/PN_SW/CalClosed'
#

tap_url = 'http://nxsa.esac.esa.int:8080/tap-server/tap'

#tap_url = 'http://archives.esac.esa.int/hsa/whsa-tap-server/tap'
#
# instrument_mode_oid=66 for small window
#query = "select t0.*,t1.* from v_exposure t0, v_public_observations t1 " + \
#    "where t0.observation_oid=t1.observation_oid and t0.filter='CalClosed' and t0.instrument_mode_oid=66"

# instrument_mode_oid=65 for large window
query = "select t0.*,t1.* from v_exposure t0, v_public_observations t1 " + \
    "where t0.observation_oid=t1.observation_oid and " +\
    "(t0.filter='CalClosed' or t0.filter='UNDEFINED') and " + \
    "t0.instrument_mode_oid=65"
#
tap_string = "{}/sync?REQUEST=doQuery&LANG=ADQL&QUERY={}".format(tap_url,query)
fail = False
#
# columns to use
useCol = ['observation_id','exposure_id','instrument_mode_oid','revolution','start_utc','end_utc',\
                 'duration','is_scientific','epic_bkg_count_rate','filter','position_angle',\
                 'ra','ra_nom','dec','lii','bii']
with requests.get(tap_string) as response:
    data = response.content
    if ('ERROR' in str(data)):
        print ("Check your query as the ADQL returned an error")
        print (tap_string)
        fail = True
if (not fail):
    votable = parse_single_table(io.BytesIO(data), pedantic=False)
    tlist = votable.to_table()
    #
    tsel = Table()
    for icol in useCol:
        if (tlist[icol].dtype == np.dtype('O')):
            tsel[icol] = tlist[icol].astype(str)
        else:
            tsel[icol] = tlist[icol]
    tsel.sort('revolution')
    tsel.write(f"{outDir}/PN_CalClosed_LW_{today}.csv",overwrite=True)
    #tsel.write(f"{outDir}/PN_CalClosed_SW_{today}.csv",overwrite=True)
    #tsel.write("{}/PN_CalClosed_FF.csv".format(outDir),overwrite=True)
