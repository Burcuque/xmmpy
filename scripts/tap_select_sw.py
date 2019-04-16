#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 11 15:39:47 2018

TAP query the XMM archive (XSA)

@author: ivaltchanov
"""
import io, os
import requests
from astropy.io.votable import parse_single_table
from astropy.table import Table
import numpy as np

home = os.path.expanduser('~')
#%%
tap_url = 'http://nxsa.esac.esa.int:8080/tap-server/tap'
#tap_url = 'http://archives.esac.esa.int/hsa/whsa-tap-server/tap'
query = "SELECT * from tables;"
#
tap_string = "{}/sync?REQUEST=doQuery&LANG=ADQL&QUERY={}".format(tap_url,query)
fail = False
with requests.get(tap_string) as response:
    data = response.content
    if ('ERROR' in str(data)):
        print ("Check your query as the ADQL returned an error")
        print (tap_string)
        fail = True
if (not fail):
    votable = parse_single_table(io.BytesIO(data), pedantic=False)
    tlist = Table(votable.to_table())
    print (tlist.info)
    print (tlist['table_name'])
#%%
#qtable = "v_public_observations"
qtable = "v_exposure"
q2 = "SELECT top 1 * from {};".format(qtable)
tap_string = "{}/sync?REQUEST=doQuery&LANG=ADQL&QUERY={}".format(tap_url,q2)
fail = False
with requests.get(tap_string) as response:
    data = response.content
    if ('ERROR' in str(data)):
        print ("Check your query as the ADQL returned an error")
        print (tap_string)
        fail = True
if (not fail):
    votable = parse_single_table(io.BytesIO(data), pedantic=False)
    ttt = Table(votable.to_table())
    print (ttt.info)
#%%
#qtable = "v_public_observations"
qtable = "v_instrument_mode"
q2 = "SELECT top 100 * from {};".format(qtable)
tap_string = "{}/sync?REQUEST=doQuery&LANG=ADQL&QUERY={}".format(tap_url,q2)
fail = False
with requests.get(tap_string) as response:
    data = response.content
    if ('ERROR' in str(data)):
        print ("Check your query as the ADQL returned an error")
        print (tap_string)
        fail = True
if (not fail):
    votable = parse_single_table(io.BytesIO(data), pedantic=False)
    tinst = Table(votable.to_table())
    print (tinst)
#%%
# now select all SW observations
q3 = "select * from v_exposure as t1, v_all_observations as t2" + \
    " where t1.instrument_mode_oid = 66 and " + \
    " t1.filter != 'CalClosed' and t1.observation_oid = t2.observation_oid"
tap_string = "{}/sync?REQUEST=doQuery&LANG=ADQL&QUERY={}".format(tap_url,q3)
fail = False
with requests.get(tap_string) as response:
    data = response.content
    if ('ERROR' in str(data)):
        print ("Check your query as the ADQL returned an error")
        print (tap_string)
        fail = True
if (not fail):
    votable = parse_single_table(io.BytesIO(data), pedantic=False)
    ttt = Table(votable.to_table())
    print (ttt.info)
#
#%%
pn_sw_obsids, indeces = np.unique(ttt['observation_id'],return_index=True)
tx = ttt["observation_id","revolution","citext","filter","duration"][indeces]

for i in np.arange(len(tx)):
    tx["observation_id"][i] = tx["observation_id"][i].decode()
    tx["citext"][i] = tx["citext"][i].decode()
    tx["filter"][i] = tx["filter"][i].decode()
    
#tx["revolution"].astype(int)
#tx["citext"].astype(str)

tx.write("{}/XMM/CAL/PN_SW_all.csv".format(home),overwrite=True)
