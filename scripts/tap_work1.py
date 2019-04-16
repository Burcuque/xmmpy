#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 11 15:39:47 2018

TAP query the XMM archive (XSA)

@author: ivaltchanov
"""
import io
import requests
from astropy.io.votable import parse_single_table
from astropy.table import Table

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
# now select PN SW CalClosed
qtable = "v_exposure"
qfilter = 'CalClosed'
q3 = "SELECT * from {} where filter=\'{}\' and instrument_mode_oid=66;".format(qtable,qfilter)
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
