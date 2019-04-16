#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 27 16:46:34 2018

Using the CalClosed results in SW mode for PN
comparing CTI v48 and v49 results 

plot results

@author: ivaltchanov
"""
import numpy as np

from astropy.table import Table
import matplotlib.pylab as plt
from matplotlib.lines import Line2D

#wdir = "/home/ivaltchanov/XMM/CAL/PN_SW"

wdir = "/xdata/xcaldata/XMM/IVAN/PN_SW/CalClosed"


v48 = Table.read(f'{wdir}/Jan2018/PN_SW_fit_results.csv')
n48 = len(v48)
v49 = Table.read(f'{wdir}/pn_sw_calclosed_fit_results_cti49z.csv')
n49 = len(v49)
# only double events selection
v49d = Table.read(f'{wdir}/pn_sw_calclosed_fit_results_cti49_doubles.csv')
n49d = len(v49d)

#%%
msize=10
fig = plt.figure(figsize=(12,8))
ax = fig.subplots()
#
# Current CTI v48 results
#
xrev = v48['revolution']
alc = 1.486 # keV
al_dif = 1000.0*(v48['al_cen']/200.0 - alc)
al_dif_err = 5.0*v48['al_cen_err']
#
mn1c = (0.162*5.888 + 0.6*5.899)/0.762 #eV probabilities from Wikipedia on Iron-55
mn1_dif = 1000.0*(v48['mn1_cen']/200.0 - mn1c)
mn1_dif_err = 5.0*v48['mn1_cen_err']
#
mn2c = 6.490 #eV probabilities from Wikipedia on Iron-55
mn2_dif = 1000.0*(v48['mn2_cen']/200.0 - mn2c)
mn2_dif_err = 5.0*v48['mn2_cen_err']
#
ax.errorbar(xrev,al_dif,yerr=(al_dif_err,al_dif_err),marker='h',markersize=15,color='cyan',
            fillstyle='none',linestyle='none')
ax.errorbar(xrev,mn1_dif,yerr=(mn1_dif_err,mn1_dif_err),marker='^',markersize=15,color='k',
            fillstyle='none',linestyle='none')
ax.errorbar(xrev,mn2_dif,yerr=(mn2_dif_err,mn2_dif_err),marker='*',markersize=15,color='y',
            fillstyle='none',linestyle='none')
#
# New CTI v49 results
#
xrev = v49['rev']
alc = 1486 # keV
al_dif = v49['al_cen'] - alc
al_dif_err = v49['al_cen_err']
#
mn1c = (0.162*5888 + 0.6*5899)/0.762 #eV probabilities from Wikipedia on Iron-55
mn1_dif = v49['mn1_cen'] - mn1c
mn1_dif_err = v49['mn1_cen_err']
#
mn2c = 6490 #eV probabilities from Wikipedia on Iron-55
mn2_dif = v49['mn2_cen'] - mn2c
mn2_dif_err = v49['mn2_cen_err']
#
# statistics
#
mean_mn1 = np.average(mn1_dif)
std_mn1 = np.std(mn1_dif)
# and all above revolution 500
irr = np.where(v49['rev'] >= 500)[0]
xmean_mn1 = np.average(mn1_dif[irr])
xstd_mn1 = np.std(mn1_dif[irr])
#
#
ax.errorbar(xrev,al_dif,yerr=(al_dif_err,al_dif_err),marker='h',markersize=15,color='cyan',
            fillstyle='full',linestyle='none')
ax.errorbar(xrev,mn1_dif,yerr=(mn1_dif_err,mn1_dif_err),marker='^',markersize=15,color='k',
            fillstyle='full',linestyle='none')
ax.errorbar(xrev,mn2_dif,yerr=(mn2_dif_err,mn2_dif_err),marker='*',markersize=15,color='y',
            fillstyle='full',linestyle='none')
#
#
# New CTI v49 results with double events
#
xrev = v49d['rev']
alc = 1486 # keV
al_dif = v49d['al_cen'] - alc
al_dif_err = v49d['al_cen_err']
#
mn1c = (0.162*5888 + 0.6*5899)/0.762 #eV probabilities from Wikipedia on Iron-55
mn1_dif = v49d['mn1_cen'] - mn1c
mn1_dif_err = v49d['mn1_cen_err']
#
mn2c = 6490 #eV probabilities from Wikipedia on Iron-55
mn2_dif = v49d['mn2_cen'] - mn2c
mn2_dif_err = v49d['mn2_cen_err']
#
# statistics
#
mean_mn1 = np.average(mn1_dif)
std_mn1 = np.std(mn1_dif)
# and all above revolution 500
irr = np.where(v49d['rev'] >= 500)[0]
xmean_mn1 = np.average(mn1_dif[irr])
xstd_mn1 = np.std(mn1_dif[irr])
#
#
#ax.errorbar(xrev,al_dif,yerr=(al_dif_err,al_dif_err),marker='h',markersize=15,color='blue',
#            fillstyle='full',linestyle='none')
ax.errorbar(xrev,mn1_dif,yerr=(mn1_dif_err,mn1_dif_err),marker='^',markersize=15,color='red',
            fillstyle='full',linestyle='none',alpha=0.3)
#
#
ax.text(450,65,r"Mean $\Delta E($Mn K$\alpha)$ = {:.2f} +/- {:.2f} eV ({:.2f} +/- {:.2f}, rev > 500)".format(mean_mn1,std_mn1,xmean_mn1,xstd_mn1),fontsize=14)
ax.set_ylim((-90.0,90.0))
ax.set_xlim((0,3300))
ax.set_title("PN Small Window in CalClosed")
ax.grid(True)
y_major_ticks = np.arange(-80, 80.0, 20)
ax.set_yticks(y_major_ticks)
ax.axhline(color='k')
ax.set_xlabel("Revolution")
ax.set_ylabel(r"$E_{obs} - E_{lab}$ (eV)")
#
cc_elements = [Line2D([0],[0],marker='h',color='w',markerfacecolor='cyan',markersize=15,label=r'Al K$\alpha$, {:.1f} keV'.format(alc)),\
                   Line2D([0],[0],marker='^',color='w',markerfacecolor='k',markersize=15,label=r'Mn K$\alpha$, {:.1f} keV'.format(mn1c)),\
                   Line2D([0],[0],marker='^',color='w',markerfacecolor='red',alpha=0.3,markersize=15,label=r'Mn K$\alpha$, {:.1f} keV doubles'.format(mn1c)),\
                   Line2D([0],[0],marker='*',color='w',markerfacecolor='y',markersize=15,label=r'Mn K$\beta$, {:.1f} keV'.format(mn2c))]

leg2 = ax.legend(handles=cc_elements,loc=3)

ax.add_artist(leg2)

plt.savefig(f'{wdir}/pn_fit_results_cti49z_with_doubles.png',dpi=100)
plt.show()
plt.close()
