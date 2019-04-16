#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  7 12:26:41 2018

@author: ivaltchanov
"""

import os
import tarfile
import requests
import shutil

import numpy as np

from astropy.io import fits
from astropy.convolution import convolve, Box1DKernel
from astropy.table import Table


#import scipy.special
#from scipy import stats

import matplotlib.pylab as plt

home = os.path.expanduser('~')
wdir = "{}/XMM/CAL/PN_SW".format(home)

#
#%%
tin = Table.read("{}/XMM/CAL/PN_SW_all.csv".format(home))
#
obsid = tin['observation_id']
nn = len(obsid)
#%%
ftype = 'SRSPEC0001'
url = 'http://nxsa.esac.esa.int/nxsa-sl/servlet'
#
os.chdir(wdir)
#
for i,qobs in enumerate(obsid):
#for i,iobs in enumerate(['0125110101','0147190101','0693781501','0721600201','0721600301','0721600401','0721600501']):
    iobs = '{:010d}'.format(qobs)
    req  = "{}/data-action-aio?obsno={}&name={}&level=PPS".format(url,iobs,ftype)
    #req  = "{}/data-action-aio?obsno={}&name={}&level=PPS".format(url,iobs,ftype)
    tarFile = "pps_{:}.tar".format(iobs)
    print ("Doing {} ({}/{}), ".format(iobs,i+1,nn),end="")
    with requests.get(req) as r:
        #r.raise_for_status() # ensure we notice bad responses
        if (b'No results' in r.content):
            print ("Skipping {}: no {} file found".format(iobs,ftype))
            continue
        with open(tarFile,"wb") as tmp:
            tmp.write(r.content)
    pn_found = False
    try:
        with tarfile.open(tarFile,'r') as tar:
            for member in tar.getmembers():
                # check member is alredy in folder
                if ('PN' in member.name and 'FTZ' in member.name):
                    pn_found = True
                    print ("Extracting {}".format(member.name))
                    f=tar.extract(member,path=wdir)
                    specFile = member.name
    except:
        print ("Cannot open the tarfile")
    #
    os.remove(tarFile)
    #
    if (pn_found):
        hdu = fits.open(specFile)
        target = hdu[0].header['OBJECT']
        rev = hdu[0].header['REVOLUT']
        xspec = hdu['SPECTRUM']
        chan = xspec.header['SPECDELT']
        xx = xspec.data['CHANNEL']*chan/1000.0
        yy = xspec.data['COUNTS']
        yy_smooth = convolve(yy, Box1DKernel(11))
        #
        m1 = xx > 4.0
        m2 = xx < 10.0
        ix = np.where(m1*m2)[0]
        ymax = np.max(yy_smooth[ix])
        m1 = xx > 5.0
        m2 = xx < 10.0
        ix = np.where(m1*m2)[0]
        countsx = np.sum(yy[ix])
        if (countsx > 5000):
            #
            #bin_means, bin_edges, binnumber = stats.binned_statistic(xx, yy,statistic='mean', bins=25)
            #bin_width = (bin_edges[1] - bin_edges[0])
            #bin_centers = bin_edges[1:] - bin_width/2
            #
            fig = plt.figure(figsize=(10,8))
            ax = fig.subplots()
            ax.loglog(xx,yy)
            ax.loglog(xx,yy_smooth,color='red')
            #plt.hlines(bin_means, bin_edges[:-1], bin_edges[1:], colors='g', lw=2,
            #           label='binned statistic of data')
            #plt.plot((binnumber - 0.5) * bin_width, yy, 'g.', alpha=0.5)
            ax.set_xlim((0.2,12.0))
            ax.set_xlabel("Energy (keV)")
            ax.set_ylabel("Counts")
            ax.set_title("PN spectrum for {} ({}_{})".format(target,rev,iobs))
            #
            #ylim = ax.get_ylim()
            #
            ay = plt.axes([0.2, 0.2, .2, .2])
            ay.plot(xx,yy_smooth)
            ay.set_xlim((4.0,10.0))
            ay.set_ylim((0,ymax*1.1))
            plt.savefig("{}/images/{}_{}_{}_spec0001.png".format(wdir,rev,iobs,target),dpi=100)
            plt.close()
        else:
            print ("Skip as counts in [5,10] keV = {} < 5000".format(countsx))
    else:
        print ("OBSID {} has no PN spectrum".format(iobs))
    #
    # remove the folders
    rootDir = member.name.split('/')[0]
    if (os.path.isdir(rootDir)):
        shutil.rmtree(rootDir,ignore_errors=True)
    #
#
