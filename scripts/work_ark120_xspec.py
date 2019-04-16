#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  2 16:43:10 2018

Work on Ark 120 and the Fe K_alpha line at 6.4 keV
Comparing MOS1&2 and PN

@author: ivaltchanov
"""
import os
import glob
from astropy.table import Table
import matplotlib.pyplot as plt
import numpy as np
#
wdir = '/xdata/xcaldata/XMM'
obsList = ['0147190101','0693781501','0721600201','0721600301','0721600401','0721600501']
#
orbit = []
emos11 = []
emos11up = []
emos11down = []
emos21 = []
emos21up = []
emos21down = []
epn1 = []
epn1up = []
epn1down = []
emos12 = []
emos12up = []
emos12down = []
emos22 = []
emos22up = []
emos22down = []
epn2 = []
epn2up = []
epn2down = []
#
for iobs in obsList:
    #
    # processing #1
    #
    xdir1 = '{}/{}/xmmsas_16.1_MOSCTI_PSF'.format(wdir,iobs)
    #xdir1 = '{}/{}/xmmsas_20170719_1539-16.1.0_release'.format(wdir,iobs)
    
    # now find the XSPEC fit results 
    #
    #MOS1
    #
    m1par = glob.glob('{}/*_EMOS1.parameter'.format(xdir1))
    if (len(m1par) != 1):
        print ("Cannot find MOS1 results in {}".format(xdir1))
        break
    else:
        orbit.append(float(os.path.basename(m1par[0]).split('_')[1]))
        t = Table.read(m1par[0],format='ascii')
        emos11.append(t['col4'].data[5])
        emos11down.append(t['col5'].data[5])
        emos11up.append(t['col6'].data[5])
    #
    # MOS2
    #
    m2par = glob.glob('{}/*_EMOS2.parameter'.format(xdir1))
    if (len(m2par) != 1):
        print ("Cannot find MOS2 results in {}".format(xdir1))
        break
    else:
        t = Table.read(m2par[0],format='ascii')
        emos21.append(t['col4'].data[5])
        emos21down.append(t['col5'].data[5])
        emos21up.append(t['col6'].data[5])
    #
    # PN
    #
    pnpar = glob.glob('{}/*_EPN.parameter'.format(xdir1))
    if (len(m1par) != 1):
        print ("Cannot find PN results in {}".format(xdir1))
        break
    else:
        t = Table.read(pnpar[0],format='ascii')
        epn1.append(t['col4'].data[5])
        epn1down.append(t['col5'].data[5])
        epn1up.append(t['col6'].data[5])
    #
    # processing #2
    #
    vers = '20180504_beta_detmap'
    xdir2 = '{}/{}/xmmsas_{}'.format(wdir,iobs,vers)
    #xdir2 = '{}/{}/xmmsas_20180403_alpha'.format(wdir,iobs)
    # now find the XSPEC fit results 
    #
    #MOS1
    #
    m1par = glob.glob('{}/*_EMOS1.parameter'.format(xdir2))
    if (len(m1par) != 1):
        print ("Cannot find MOS1 results in {}".format(xdir2))
        break
    else:
        t = Table.read(m1par[0],format='ascii')
        emos12.append(t['col4'].data[5])
        emos12down.append(t['col5'].data[5])
        emos12up.append(t['col6'].data[5])
    #
    # MOS2
    #
    m2par = glob.glob('{}/*_EMOS2.parameter'.format(xdir2))
    if (len(m2par) != 1):
        print ("Cannot find MOS2 results in {}".format(xdir2))
        break
    else:
        t = Table.read(m2par[0],format='ascii')
        emos22.append(t['col4'].data[5])
        emos22down.append(t['col5'].data[5])
        emos22up.append(t['col6'].data[5])
    #
    # PN
    #
    pnpar = glob.glob('{}/*_EPN.parameter'.format(xdir2))
    if (len(m1par) != 1):
        print ("Cannot find PN results in {}".format(xdir2))
        break
    else:
        t = Table.read(pnpar[0],format='ascii')
        epn2.append(t['col4'].data[5])
        epn2down.append(t['col5'].data[5])
        epn2up.append(t['col6'].data[5])
#%%
# now plotting
#
idx = np.arange(len(orbit))
fig = plt.figure(figsize=(10,8))
ax1 = fig.add_subplot(111)
#ax1.errorbar(orbit,emos1,yerr=[np.abs(emos1down),emos1up],fmt='-o',label='MOS1')
#ax1.errorbar(np.add(orbit,0.1),emos2,yerr=[np.abs(emos2down),emos2up],fmt='-o',label='MOS2')
#ax1.errorbar(np.add(orbit,0.2),epn,yerr=[np.abs(epndown),epnup],fmt='-o',label='PN')
#ax1.set_xlabel("Orbit")
ax1.errorbar(idx,emos11,yerr=[np.abs(emos11down),emos11up],fmt='o',ms=10,color='blue',label='MOS1 CTI PSF')
ax1.errorbar(idx+0.1,emos21,yerr=[np.abs(emos21down),emos21up],fmt='o',ms=10,color='red',label='MOS2 CTI PSF')
ax1.errorbar(idx+0.2,epn1,yerr=[np.abs(epn1down),epn1up],fmt='o',ms=10,color='green',label='PN CTI PSF')
#ax1.errorbar(idx,emos11,yerr=[np.abs(emos11down),emos11up],fmt='o',ms=10,mfc='blue',label='MOS1 CTI PSF')
#ax1.errorbar(idx+0.1,emos21,yerr=[np.abs(emos21down),emos21up],fmt='o',ms=10,mfc='red',label='MOS2 CTI PSF')
#ax1.errorbar(idx+0.2,epn1,yerr=[np.abs(epn1down),epn1up],fmt='o',ms=10,mfc='green',label='PN CTI PSF')
#
#ax1.errorbar(idx,emos12,yerr=[np.abs(emos12down),emos12up],marker='s',ms=10,mfc='cyan',label='MOS1 %s'%vers)
#ax1.errorbar(idx+0.1,emos22,yerr=[np.abs(emos22down),emos22up],marker='s',ms=10,mfc='pink',label='MOS2 %s'%vers)
#ax1.errorbar(idx+0.2,epn2,yerr=[np.abs(epn2down),epn2up],marker='s',ms=10,mfc='lime',label='PN %s'%vers)
ax1.errorbar(idx,emos12,yerr=[np.abs(emos12down),emos12up],marker='s',ms=10,color='cyan',label='MOS1 %s'%vers)
ax1.errorbar(idx+0.1,emos22,yerr=[np.abs(emos22down),emos22up],marker='s',ms=10,color='pink',label='MOS2 %s'%vers)
ax1.errorbar(idx+0.2,epn2,yerr=[np.abs(epn2down),epn2up],marker='s',ms=10,color='lime',label='PN %s'%vers)
#
ax1.set_xlabel("Obs index")
for ik in idx:
    ax1.text(ik,6.3,'{}'.format(int(orbit[ik])))
ax1.set_ylabel(r'Fe K$_\alpha$ line centre (keV)')
ax1.grid(True)
ax1.legend(loc=4)
plt.title(r'Ark 120 results on Fe K$_\alpha$ line')
#plt.title("{} (n={})".format(fitfile,nt))
#plt.show()
plt.savefig('Ark120_xspec_FeK_comparison_moscti_psf.png',dpi=300)
plt.close()
    

    