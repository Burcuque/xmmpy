#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 28 14:35:17 2019

@author: ivaltchanov
"""
import numpy as np

from astropy.table import Table
import matplotlib.pylab as plt
import seaborn as sns
sns.set(style="white")

plt.rc('text', usetex=False)
plt.rc('font', family='serif')

#%%
wdir = '/xdata/xcaldata/XMM/IVAN/PN_LW'

obslist = {"0727360201": "S003", "0203720201": "S020"}
colx1 = {"0727360201": "black", "0203720201": "red"}
colx2 = {"0727360201": "black", "0203720201": "red"}
colx3 = {"0727360201": "black", "0203720201": "red"}

fig, axs = plt.subplots(3,1,sharex=True,figsize=(15,10))

for iobs in obslist.keys():
    #
    #fout = open(f"{wdir}/pn_sw_calclosed_fit_results_nocti.csv","a")
    #t = Table.read(f"{wdir}/pn_lw_calclosed_fit_results_{iobs}_{obslist[iobs]}.csv")
    #tx = Table.read(f"{wdir}/{iobs}_{obslist[iobs]}_xspec_single_results.csv")
    tx = Table.read(f"{wdir}/CalClosed/{iobs}_{obslist[iobs]}_xspec_results_cti49.csv")
    #
    alka = 1.486
    mnka = 5.8988
    cuka = 8.04
    #
    #d_al = 1000*(t['al_cen'].data - alka)
    #ix = np.where(t['al_cen_err'].data == "None")[0]
    #al_cen_err = t['al_cen_err'].data.copy()
    #al_cen_err[ix] = 0.0
    #d_al_err = 1000*al_cen_err.astype(float)
    #
    #d_mn = 1000*(t['mn1_cen'].data - mnka)
    #ix = np.where(t['mn1_cen_err'].data == "None")[0]
    #mn_cen_err = t['mn1_cen_err'].data.copy()
    #mn_cen_err[ix] = 0.0
    #d_mn_err = 1000*mn_cen_err.astype(float)
    ##
    #d_cu = 1000*(t['cu_cen'].data - cuka)
    #ix = np.where(t['cu_cen_err'].data == "None")[0]
    ##cu_cen_err = t['cu_cen_err'].data.copy()
    #cu_cen_err[ix] = 0.0
    #d_cu_err = 1000*cu_cen_err.astype(float)   
    #
    # now the XSPEC results
    d2_al = 1000*(tx['line1'].data - alka)
    d2_al_err = 1000*(tx['line1err'].data)
    d2_mn = 1000*(tx['line2'].data - mnka)
    d2_mn_err = 1000*(tx['line2err'].data)
    d2_cu = 1000*(tx['line4'].data - cuka)
    d2_cu_err = 1000*(tx['line4err'].data)
    # the boresight has ccdnr == 4.2
    #
    iq = np.where(t['ccdnr'].data == 4.2)[0]
    iq2 = np.where(tx['ccd'].data == 4.2)[0]
    #axs[0].errorbar(t['ccdnr'].data,d_al,yerr=(d_al_err,d_al_err),fmt='o',color=colx1[iobs],\
    #   label=rf"{iobs}: Al K$\alpha$ (1.49 keV)")
    #axs[0].errorbar(t['ccdnr'].data[iq],d_al[iq],yerr=(d_al_err[iq],d_al_err[iq]),fmt='s',color=colx1[iobs],\
    #   label="",markersize=15,fillstyle='none')
    #
    axs[0].errorbar(tx['ccd'].data,d2_al,yerr=(d2_al_err,d2_al_err),fmt='o',\
       markersize=7,fillstyle='none',color=colx1[iobs],\
       label=rf"{iobs}: XSPEC Al K$\alpha$ (1.49 keV)")
    axs[0].errorbar(tx['ccd'].data[iq2],d2_al[iq2],yerr=(d2_al_err[iq2],d2_al_err[iq2]),fmt='s',color=colx1[iobs],\
       label="",markersize=15,fillstyle='none')
    axs[0].set_ylim((-60.0,60.0))
    axs[0].xaxis.set_major_locator(plt.MaxNLocator(13))
    axs[0].grid(True)
    axs[0].legend()
    axs[0].axhline(0.0,color='k',ls='dashed')
    #axs[0].set_title(f"PN Large Window OBS_ID {iobs}")
    #axs[1].errorbar(t['ccdnr'].data,d_mn,yerr=(d_mn_err,d_mn_err),fmt='o',color=colx2[iobs],\
    #   label=fr"{iobs}: Mn K$\alpha$ (5.9 keV)")
    #axs[1].errorbar(t['ccdnr'].data[iq],d_mn[iq],yerr=(d_mn_err[iq],d_mn_err[iq]),fmt='s',color=colx2[iobs],\
    #   label="",markersize=15,fillstyle='none')
    axs[1].errorbar(tx['ccd'].data,d2_mn,yerr=(d2_mn_err,d2_mn_err),fmt='o',\
       markersize=7,fillstyle='none',color=colx2[iobs],\
       label=fr"{iobs}: XSPEC Mn K$\alpha$ (5.9 keV)")
    axs[1].errorbar(tx['ccd'].data[iq2],d2_mn[iq2],yerr=(d2_mn_err[iq2],d2_mn_err[iq2]),fmt='s',color=colx2[iobs],\
       label="",markersize=15,fillstyle='none')
    axs[1].set_ylim((-60.0,60.0))
    axs[1].xaxis.set_major_locator(plt.MaxNLocator(13))
    axs[1].grid(True)
    axs[1].legend()
    axs[1].axhline(0.0,color='k',ls='dashed')
    #axs[2].errorbar(t['ccdnr'].data,d_cu,yerr=(d_cu_err,d_cu_err),fmt='o',color=colx3[iobs],\
    #   label=fr"{iobs}: Cu K$\alpha$ (8.0 keV)")
    #axs[2].errorbar(t['ccdnr'].data[iq],d_cu[iq],yerr=(d_cu_err[iq],d_cu_err[iq]),fmt='s',color=colx3[iobs],\
    #   label="",markersize=15,fillstyle='none')
    axs[2].errorbar(tx['ccd'].data,d2_cu,yerr=(d2_cu_err,d2_cu_err),fmt='o',\
       markersize=7,fillstyle='none',color=colx3[iobs],\
       label=fr"{iobs}: XSPECT Cu K$\alpha$ (8.0 keV)")
    axs[2].errorbar(tx['ccd'].data[iq2],d2_cu[iq2],yerr=(d2_cu_err[iq2],d2_cu_err[iq2]),fmt='s',color=colx3[iobs],\
       label="",markersize=15,fillstyle='none')
    axs[2].set_ylim((-60.0,60.0))
    axs[2].xaxis.set_major_locator(plt.MaxNLocator(13))
    axs[2].grid(True)
    axs[2].legend()
    axs[2].axhline(0.0,color='k',ls='dashed')
    axs[2].set_xlabel("CCD number")
#
axs[0].set_title("PN LW Analysis, EPN_CTI_0048.CCF")
plt.subplots_adjust(wspace=0, hspace=0)
#plt.savefig(f"{wdir}/pn_lw_xspec_fit_results_cti48.png")
plt.show()
plt.close()
