#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 27 16:46:34 2018

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
            'ngc5548': 'p'}

t = Table.read('{}/Fe_diff_results.csv'.format(wdir))
nt = len(t)

msize=10
fig = plt.figure(figsize=(12,8))
ax = fig.subplots()

for q in np.arange(nt):
    colx = 'blue'
    fstyle = 'none'
    if (t['inst'].data[q] == 'mos1'):
        colx = 'blue'
        labx = 'MOS1'
    if (t['inst'].data[q] == 'mos2'):
        colx = 'green'
        labx = 'MOS2'
    if (t['inst'].data[q] == 'pn'):
        colx = 'red'
        labx = 'PN'
    if (t['mode'].data[q] == 1):
        fstyle = 'full'
    ax.errorbar([t['rev'].data[q]],[t['diff'].data[q]],yerr=([t['diffErr'].data[q]],[t['diffErr'].data[q]]),\
            fmt=target_mark[t['target'].data[q]],\
            markersize=msize,fillstyle=fstyle,color=colx)

ax.set_ylim((-90.0,90.0))
ax.set_xlim((2000,3300))
ax.set_title(r"Fe K$_\alpha$ analysis")
ax.grid(True)
y_major_ticks = np.arange(-80, 80.0, 20)
ax.set_yticks(y_major_ticks)
ax.axhline(color='k')
ax.set_xlabel("Revolution")
ax.set_ylabel("Line Centre Differnece (eV)")
#
legend_elements = [Line2D([0],[0],marker='o',color='w',markerfacecolor='b',markersize=10,label='MOS1'),\
                   Line2D([0],[0],marker='o',color='w',markerfacecolor='g',markersize=10,label='MOS2'),\
                   Line2D([0],[0],marker='o',color='w',markerfacecolor='r',markersize=10,label='PN')]
leg0 = ax.legend(handles=legend_elements,loc=0)
#
target_elements = [Line2D([0],[0],marker='o',color='w',markerfacecolor='b',markersize=10,label='NGC4151'),\
                   Line2D([0],[0],marker='D',color='w',markerfacecolor='b',markersize=10,label='NGC3227'),\
                   Line2D([0],[0],marker='s',color='w',markerfacecolor='b',markersize=10,label='MRK1048'),\
                   Line2D([0],[0],marker='v',color='w',markerfacecolor='b',markersize=10,label='NGC3783'),\
                   Line2D([0],[0],marker='<',color='w',markerfacecolor='b',markersize=10,label='NGC5506'),\
                   Line2D([0],[0],marker='>',color='w',markerfacecolor='b',markersize=10,label='NGC3516'),\
                   Line2D([0],[0],marker='*',color='w',markerfacecolor='b',markersize=10,label='MCG-5-23-16'),\
                   Line2D([0],[0],marker='p',color='w',markerfacecolor='b',markersize=10,label='NGC5548'),\
                   Line2D([0],[0],marker='^',color='w',markerfacecolor='b',markersize=10,label='NGC4593')]
leg1 = ax.legend(handles=target_elements,loc=4)
ax.add_artist(leg0)
ax.add_artist(leg1)
plt.savefig('{}/Fe_fit_results_2000.png'.format(wdir),dpi=100)
plt.show()
plt.close()
