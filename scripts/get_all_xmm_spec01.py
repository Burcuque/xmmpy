#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May  6 17:31:41 2019

Using the TAP interface, loop over all XMM public observations, 
download the PPS in a temporary file, extract the spectrum from the brightest source 
and save it to a FITS file

Wil only use PN observations of modes between 60 and 69.

@author: ivaltchanov
"""
import os
import tarfile
import requests
import pandas as pd
from io import StringIO
#import numpy as np

import matplotlib.pylab as plt
import matplotlib as mpl
mpl.use('Agg')

from astropy.io import fits
#from astropy.stats import histogram

#%%
def get_all_xmm_obsids(min_expo=20000,verbose=True):
    #
    # min_expo os the minimum exposure time (duration) to consider
    #
    tap_url = "http://nxsa.esac.esa.int/tap-server/tap/sync?REQUEST=doQuery&LANG=ADQL&FORMAT=csv&"
    query = f"QUERY=select t2.observation_id,t2.ra,t2.dec,t2.revolution,t2.duration,t3.mode_friendly_name,t1.filter" + \
    " from v_exposure as t1, v_public_observations as t2, v_instrument_mode as t3" + \
    " WHERE " + \
    f"(t1.duration >= {min_expo}) and " + \
    "(t1.observation_oid = t2.observation_oid) and " + \
    "(t1.instrument_mode_oid BETWEEN 60 and 69) and " + \
    "(t1.instrument_mode_oid = t3.instrument_mode_oid)"
    #"((t1.filter ILIKE \"Thin1\") or (t1.filter ILIKE \"Thick\") or (t1.filter ILIKE \"Medium\"))"
    xreq = tap_url + query
    with requests.get(xreq) as r:
        jx = r.content.decode()
    #
    if ('ERROR' in jx):
        print ("ADQL query returns error:",jx)
        raise Exception
    buff = StringIO(jx)
    df = pd.read_csv(buff, sep=",")
    df.rename(columns={'observation_id': 'obs_id','mode_friendly_name': 'obs_mode','revolution': 'rev'},inplace=True)
    #
    # now select those with 'MEDIUM', 'THICK', 'THIN2', 'THIN1'
    #
    df['filter'] = df['filter'].str.upper()
    #
    out_df = df.query("filter in [\"THIN1\",\"THIN2\",\"THICK\",\"MEDIUM\"]")
    if (verbose):
        print (f"Found {len(out_df)} XMM-Newton EPIC-pn public observations in the archive")
    return out_df
#%%
def get_xmm_spec(obs_id,src_index=1, tmpdir=".",expo_limit = 10.0,\
                 verbose=False, outputFile="output.fits"):
    #
    # get the PPS spectrum product for source with src_index from an XMM observation obs_id
    #
    # special version to reuse the same tar filename 
    # will skip spectra of less than expo_limit (default 10 ks)
    #
    ftype = "SRSPEC"
    extn = "FTZ"
    inst = "PN"
    url = 'http://nxsa.esac.esa.int/nxsa-sl/servlet'
    #req  = f"{url}/data-action-aio?obsno={iobs}&name={ftype}&level=PPS"
    req  = f"{url}/data-action-aio?obsno={obs_id}&extension={extn}&name={ftype}&instname={inst}&level=PPS"
    #
    
    tarFile = f"{tmpdir}/pps_spec.tar"
    if (verbose):
        print (f"Downloading tar file for {obs_id} from NXSA... wait")
    with requests.get(req) as r:
        if (b'No results' in r.content):
            print (f"No PPS spectral products found for {obs_id}")
            return None
        else:
            with open(tarFile,"wb") as tmp:
                tmp.write(r.content)
    if (verbose):
        print (f"{inst} spectra in all bands saved to {tarFile}")
    # 
    # Extract the spectrum
    #
    hdu = None
    try:
        with tarfile.open(tarFile,'r') as tar:
            src = f'SRSPEC{src_index:04}'
            for member in tar.getmembers():
                if (src in member.name):
                    print (f"Extracting {member.name}")
                    f=tar.extract(member,path=tmpdir)
                    hdu = fits.open(f"{tmpdir}/{member.name}")
                    texp = hdu[1].header['EXPOSURE']/1000.0
                    if (texp >= expo_limit):
                        hdu.writeto(outputFile,overwrite=True)
                        print (f"Will save PN spectrum 01 for {obs_id} as the total exposure after GTI is {texp:.1f} ks (>= {expo_limit} ks)")                        
                    else:
                        print (f"Will skip PN spectrum 01 for {obs_id} as the total exposure after GTI is {texp:.1f} ks (< {expo_limit} ks)")
                    break
    except:
        hdu = fits.open(tarFile)
        texp = hdu[1].header['EXPOSURE']/1000.0
        if (texp >= expo_limit):
            hdu.writeto(outputFile,overwrite=True)
            #print (f"Will save PN spectrum 01 for {obs_id} as the total exposure after GTI is {texp:.1f} ks (>= {expo_limit} ks)")                        
        else:
            print (f"Will skip PN spectrum 01 for {obs_id} as the total exposure after GTI is {texp:.1f} ks (< {expo_limit} ks)")
            return None
    #
    if hdu == None:
        print (f"No PN spectral product found for source {src}")
    return hdu

#%%
def plot_xmm_spec(spec_hdu, pngName = None):
    #
    #
    #
    plot_title = f"OBS_ID: {spec_hdu[0].header['EXP_ID']}, " + \
        f"Target: {spec_hdu[0].header['OBJECT']}, \n" + \
        f"expo: {spec_hdu[1].header['EXPOSURE']/1000.0:.1f} ks"
    spec = tt['SPECTRUM']
    binsize = spec.header['SPECDELT']
    x = spec.data['CHANNEL']*binsize/1000.0
    y = spec.data['COUNTS']
    xlims  = (0.2, 12.0)
    #ysel = y[np.where((x <= xlims[1]) & (x >= xlims[0]))[0]]
    #a = histogram(ysel,bins='blocks')
    #
    #
    fig = plt.figure(figsize=(15,5))
    ax = fig.subplots()
    #ax.semilogx(x,y)
    #ax.step(a[1][1:],a[0],color='red')
    ax.step(x,y,color='red')
    ax.set_xlim(xlims)
    ax.set_xlabel("Energy (keV)")
    ax.set_ylabel("Counts")
    ax.set_title(plot_title)
    ax.grid()
    if (pngName != None):
        plt.savefig(pngName,dpi=100)
        plt.close()
    
    
#%%
home = os.path.expanduser('~')

tmpdir = "/lhome/ivaltchanov/XMM-data"
if (not os.path.isdir(tmpdir)):
    print ("Please set the tmpdir to an aleady existing folder. It is needed.")
    raise FileNotFoundError
#
output_dir = "/lhome/ivaltchanov/XMM-data/spectra01"
if (not os.path.isdir(output_dir)):
    print ("Please set the output_dir to an aleady existing folder. It is needed.")
    raise FileNotFoundError

#%%
obs_tab = get_all_xmm_obsids()
#%%
# limit to spectra with at least 10 ks after GTI
#
fout = open(f'{tmpdir}/xmm_spectra_for_ml.csv','w')
#
print ("exp_id,rev,object,submode,filter,duration,exposure,filename",file=fout)
nt = obs_tab.shape[0]

for i,iobs in enumerate(obs_tab["obs_id"]):
    sobs = f"{iobs:010}"
    print (f"Doing {sobs} ({i+1}/{nt})")
    output_spec_name = f"{output_dir}/{sobs}_spec01.fits"
    pngName = f"{output_dir}/{sobs}_spec01.png"
    tt = get_xmm_spec(sobs,verbose=False,expo_limit=10.0,tmpdir=tmpdir,outputFile=output_spec_name)
    if (tt != None):
        print (f"{tt[0].header['EXP_ID']}, " + \
                  f"{tt[0].header['REVOLUT']}," + \
                  f"{tt[0].header['OBJECT']}," + \
                  f"{tt[0].header['SUBMODE']}," + \
                  f"{tt[0].header['FILTER']}," + \
                  f"{tt[0].header['DURATION']/1000.0:.1f}," + \
                  f"{tt[1].header['EXPOSURE']/1000.0:.1f}," + \
                  f"{os.path.basename(output_spec_name)}",file=fout)
        plot_xmm_spec(tt,pngName=pngName)
#
print ("all done")
fout.close()
#
#%%
#



