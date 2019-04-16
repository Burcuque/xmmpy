#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov  9 10:36:31 2018

Modify CCF file "by hand "

@author: ivaltchanov
"""

from astropy.io import fits
import numpy as np
#
#%%
# the ascii dumps should follow the format as in Michael Smith's folder:
# ~msmith/CCFdev/ccfdev/packages/epn-cti/src/
# 


wdir = "/home/ivaltchanov/IVAN/ccfdev"

ver = "0049"
test = "test53"
#ver = "0048"
#ccf_file = f'EPN_CTI_{ver}.CCF_{test}'
ccf_file = f'EPN_CTI_{ver}.CCF'
hdu = fits.open(f"{wdir}/{ccf_file}")
#%%
# first, dump the LONG_TERM_CTI extension
#
extname = 'LONG_TERM_CTI'
ltc = hdu[extname]
nt = len(ltc.data)
fout = open(f'{wdir}/long_term_cti_v{ver}_{test}.dat','w')

for j in np.arange(nt):
    print ("{:3} : {:3} :   {:6.4f} : ".format(ltc.data['MODE_ID'][j],ltc.data['CCD_ID'][j],\
           ltc.data['ENERGY'][j]),end="",file=fout)
    tcoeff = ltc.data['T_COEFF']
    nrow,ncol = tcoeff.shape
    for l in np.arange(ncol):
        print ("  {:3.6e}".format(tcoeff[j,l]),end="",file=fout)
    print (" :   1",file=fout)
    #
#
    
fout.close()
#%%
# second, dump the 'COMB_EVT_OFFSET' extension
#
extname = 'COMB_EVT_OFFSET'
ltc = hdu[extname]
nt = len(ltc.data)
fout = open(f'{wdir}/comb_evt_offset_v{ver}_{test}.dat','w')

for j in np.arange(nt):
    print ("{:3} : {:3} : {:3} :   {:6.4f} : {:6.1f}".format(ltc.data['MODE_ID'][j],ltc.data['CCD_ID'][j],\
           ltc.data['PATTERN'][j],ltc.data['ENERGY'][j], ltc.data['SHIFT'][j]),file=fout)
    #
#
    
fout.close()
#%%
#
# Now the system testgin CCFs, dump them to ASCII and then use diff
#
ccfdir = "/ccf/valid"
#ccfdir = "/ccf/priv/dt"
wdir = "/home/ivaltchanov/IVAN/ccfdev"
vers = "0048"

ccf_file = f'EPN_CTI_{vers}.CCF'
hdu = fits.open(f"{ccfdir}/{ccf_file}")

extname = 'LONG_TERM_CTI'
ltc = hdu[extname]
nt = len(ltc.data)
fout = open(f'{wdir}/long_term_cti_v{vers}_testing.dat','w')

for j in np.arange(nt):
    print ("{:3} : {:3} :   {:6.4f} : ".format(ltc.data['MODE_ID'][j],ltc.data['CCD_ID'][j],\
           ltc.data['ENERGY'][j]),end="",file=fout)
    tcoeff = ltc.data['T_COEFF']
    nrow,ncol = tcoeff.shape
    for l in np.arange(ncol):
        print ("  {:3.6e}".format(tcoeff[j,l]),end="",file=fout)
    print (" :   1",file=fout)
    #
#    
fout.close()
# second, dump the 'COMB_EVT_OFFSET' extension
#
extname = 'COMB_EVT_OFFSET'
ltc = hdu[extname]
nt = len(ltc.data)
fout = open(f'{wdir}/comb_evt_offset_v{vers}_testing.dat','w')

for j in np.arange(nt):
    print ("{:3} : {:3} : {:3} :   {:6.4f} : {:6.1f}".format(ltc.data['MODE_ID'][j],ltc.data['CCD_ID'][j],\
           ltc.data['PATTERN'][j],ltc.data['ENERGY'][j], ltc.data['SHIFT'][j]),file=fout)
    #
#
hdu.close()
fout.close()
