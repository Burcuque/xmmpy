#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 27 16:46:34 2018

Using the Fe fit in XSPEC for AGNs and CalClosed results in SW mode for PN
plot results

Only plot the results for PN Small Window mode
@author: ivaltchanov
"""
import numpy as np

from astropy.table import Table
import matplotlib.pylab as plt
from matplotlib.lines import Line2D

#wdir = "/home/ivaltchanov/XMM/CAL/PN_SW"

wdir = "/xdata/xcaldata/XMM/IVAN/PN_SW"

target_mark = {'ngc4151': 'o', 'ngc3227': 'D', 'mrk1048': "s", 'ngc3783': 'v',\
            'ngc4593': '^', 'ngc5506': '<', 'ngc3516': '>', 'mcg-5-23-16': '*',\
            'ngc5548': 'p', 'ngc2992': 'H'}

tcc = Table.read('{}/../PN_calclosed/Jan2018/PN_SW_fit_results.csv'.format(wdir))
ncc = len(tcc)

#%%
t = Table.read('{}/Fe_diff_results.csv'.format(wdir))
nt = len(t)

msize=10
fig = plt.figure(figsize=(12,8))
ax = fig.subplots()

for q in np.arange(nt):
    colx = 'red'
    fstyle = 'none'
    labx = 'PN SW'
    if (t['inst'].data[q] != 'pn' or t['mode'].data[q] == 1):
        continue
    ax.errorbar([t['rev'].data[q]],[t['diff'].data[q]],yerr=([t['diffErr'].data[q]],[t['diffErr'].data[q]]),\
            fmt=target_mark[t['target'].data[q]],\
            markersize=msize,fillstyle=fstyle,color=colx)

ax.set_ylim((-90.0,90.0))
ax.set_xlim((0,3300))
ax.set_title("Energy scale analysis")
ax.grid(True)
y_major_ticks = np.arange(-80, 80.0, 20)
ax.set_yticks(y_major_ticks)
ax.axhline(color='k')
ax.set_xlabel("Revolution")
ax.set_ylabel("Line Centre Differnece (eV)")
#
# now the CalClosed results
#
xrev = tcc['revolution']
alc = 1.486 # keV
al_dif = 1000.0*(tcc['al_cen']/200.0 - alc)
al_dif_err = 5.0*tcc['al_cen_err']
#
mn1c = (0.162*5.888 + 0.6*5.899)/0.762 #eV probabilities from Wikipedia on Iron-55
mn1_dif = 1000.0*(tcc['mn1_cen']/200.0 - mn1c)
mn1_dif_err = 5.0*tcc['mn1_cen_err']
#
mn2c = 6.490 #eV probabilities from Wikipedia on Iron-55
mn2_dif = 1000.0*(tcc['mn2_cen']/200.0 - mn2c)
mn2_dif_err = 5.0*tcc['mn2_cen_err']
#
ax.errorbar(xrev,al_dif,yerr=(al_dif_err,al_dif_err),marker='h',markersize=15,color='cyan',
            fillstyle='full',linestyle='none')
ax.errorbar(xrev,mn1_dif,yerr=(mn1_dif_err,mn1_dif_err),marker='^',markersize=15,color='k',
            fillstyle='full',linestyle='none')
ax.errorbar(xrev,mn2_dif,yerr=(mn2_dif_err,mn2_dif_err),marker='*',markersize=15,color='y',
            fillstyle='full',linestyle='none')

legend_elements = [Line2D([0],[0],marker='o',color='w',markerfacecolor='r',markersize=10,fillstyle='none',label='PN')]
leg0 = ax.legend(handles=legend_elements,loc=2,numpoints=1)
#legend_elements = [Line2D([0],[0],marker='o',color='w',markerfacecolor='b',markersize=10,label='MOS1'),\
#                   Line2D([0],[0],marker='o',color='w',markerfacecolor='g',markersize=10,label='MOS2'),\
#                   Line2D([0],[0],marker='o',color='w',markerfacecolor='r',markersize=10,label='PN')]
#leg0 = ax.legend(handles=legend_elements,loc=2)
#
target_elements = [Line2D([0],[0],marker='o',color='w',markerfacecolor='b',markersize=10,label='NGC4151'),\
                   Line2D([0],[0],marker='D',color='w',markerfacecolor='b',markersize=10,label='NGC3227'),\
                   Line2D([0],[0],marker='s',color='w',markerfacecolor='b',markersize=10,label='MRK1048'),\
                   Line2D([0],[0],marker='v',color='w',markerfacecolor='b',markersize=10,label='NGC3783'),\
                   Line2D([0],[0],marker='<',color='w',markerfacecolor='b',markersize=10,label='NGC5506'),\
                   Line2D([0],[0],marker='>',color='w',markerfacecolor='b',markersize=10,label='NGC3516'),\
                   Line2D([0],[0],marker='*',color='w',markerfacecolor='b',markersize=10,label='MCG-5-23-16'),\
                   Line2D([0],[0],marker='p',color='w',markerfacecolor='b',markersize=10,label='NGC5548'),\
                   Line2D([0],[0],marker='H',color='w',markerfacecolor='b',markersize=10,label='NG2992'),\
                   Line2D([0],[0],marker='^',color='w',markerfacecolor='b',markersize=10,label='NGC4593')]
leg1 = ax.legend(handles=target_elements,loc=4,numpoints=1)

cc_elements = [Line2D([0],[0],marker='h',color='w',markerfacecolor='cyan',markersize=15,label=r'Al K$\alpha$, {:.3f} keV'.format(alc)),\
                   Line2D([0],[0],marker='^',color='w',markerfacecolor='k',markersize=15,label=r'Mn K$\alpha$, {:.3f} keV'.format(mn1c)),\
                   Line2D([0],[0],marker='*',color='w',markerfacecolor='y',markersize=15,label=r'Mn K$\beta$, {:.3f} keV'.format(mn2c))]

leg2 = ax.legend(handles=cc_elements,loc=3,numpoints=1)

ax.add_artist(leg0)
ax.add_artist(leg1)
ax.add_artist(leg2)

plt.savefig('{}/ccf48_results_pn_sw.png'.format(wdir),dpi=100)
plt.show()
plt.close()
