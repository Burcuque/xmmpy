#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 20 12:27:08 2018

Plot the line fit results, the difference is with respect to the 
redshifted Fe Kalpha line

@author: ivaltchanov
"""
import os

import numpy as np

from astropy.table import Table

import matplotlib.pylab as plt
#from matplotlib.patches import Patch
from matplotlib.lines import Line2D

plt.rc('text', usetex=False)
plt.rc('font', family='serif')
#
#target = "ngc4593"
#target = "ngc4151"
#target = "mrk1048"
#target = "ngc3227"
target = "ngc3783"
#target = "ngc5506"
#target = "mcg-5-23-16"
#target = "ngc5548"
#target = "ngc3516"
#target = "WR140"
#target = "ngc2992"
#
# redshifts
#
redshift = {'ngc4151': 0.003262, 'ngc3227': 0.00386, 'mrk1048': 0.0427, 'ngc3783': 0.009755,\
            'ngc4593': 0.008344, 'ngc5506': 0.00589, 'mcg-5-23-16': 0.008226, 'ngc3516': 0.008816,\
            'ngc5548': 0.01627, 'WR140': 0.0, 'ngc2992': 0.007296}
feK = 6.399 # the weitghed mean of Kalpha_2 at 6.3908 (intensity 50) and Kalpha_1 at 6.40308 (intensity 100)
lineX =  feK/(1.0 + redshift[target]) # redshifted line position
if (target == 'WR140'):
    lineX = 6.697 # see Sugawara et al. (2015)
#
wdir = "/xdata/xcaldata/XMM/IVAN/PN_SW/{}".format(target)
#wdir = "/home/ivaltchanov/XMM/CAL/PN_SW/{}".format(target)
#wdir = "/home/ivaltchanov/XMM/CAL/PN_SW/WR140"
os.chdir(wdir)

fout = open('{}/../Fe_diff_results.csv'.format(wdir),'a')

#%%
t = Table.read('output_xspec.csv')
nt = len(t)
#
rev = t['revol'].data
ff = t['ff'].data
inst = t['instrument'].data
line = t['line'].data
lineErrUp = t['lineErr_up'].data
lineErrLo = t['lineErr_lo'].data
rchi2 = t["chi2r"].data
ymax = lineX
#ymax = np.max(line + lineErrUp)
ymin = np.min(line - lineErrLo)
#ymin = ymax - 0.
#
#%%
msize=10
fig = plt.figure(figsize=(12,8))
ax = fig.subplots(3,sharex=True)

istart = 0
ratmin = -80.0
ratmax = 80.0
#
for j in np.unique(rev):
    ix = np.where(rev == j)[0]
    nx = len(ix)
    #if (nx < 2):
    #    continue
    for q in np.arange(nx):
        colx = 'blue'
        fstyle = 'none'
        if (inst[ix[q]] == 'mos1'):
            colx = 'blue'
            labx = 'MOS1'
            #ratio = line[ix[q]]/lineX
            #ratioErr = ratio*np.sqrt(np.power(lineXerr/lineX,2) + np.power(lineErrUp[ix[q]]/line[ix[q]],2))
            #ratioErr = 0.0
            diff = line[ix[q]] - lineX
            diffErr = lineErrUp[ix[q]]
        if (inst[ix[q]] == 'mos2'):
            colx = 'green'
            labx = 'MOS2'
            #ratio = line[ix[q]]/lineX
            #ratioErr = ratio*np.sqrt(np.power(lineXerr/lineX,2) + np.power(lineErrUp[ix[q]]/line[ix[q]],2))
            diff = line[ix[q]] - lineX
            diffErr = lineErrUp[ix[q]]
            #diffErr = np.sqrt(np.power(lineXerr,2) + np.power(lineErrUp[ix[q]],2))
        if (inst[ix[q]] == 'pn'):
            colx = 'red'
            labx = 'PN'
            #ratio = line[ix[q]]/lineX
            #ratioErr = ratio*np.sqrt(np.power(lineXerr/lineX,2) + np.power(lineErrUp[ix[q]]/line[ix[q]],2))
            diff = line[ix[q]] - lineX
            diffErr = lineErrUp[ix[q]]
            #diffErr = np.sqrt(np.power(lineXerr,2) + np.power(lineErrUp[ix[q]],2))
        if (ff[ix[q]] == 1):
            fstyle = 'full'
        #
        #if (istart == 0):
        #qq = ax[0].errorbar([q+istart],[line[ix[q]]],yerr=([lineErrLo[ix[q]]],[lineErrUp[ix[q]]]),fmt='o',\
        #       markersize=msize,fillstyle=fstyle,label=labx,color=colx)
        #else:
        ax[0].errorbar([q+istart],[line[ix[q]]],yerr=([lineErrLo[ix[q]]],[lineErrUp[ix[q]]]),fmt='o',\
                   markersize=msize,fillstyle=fstyle,color=colx)
        #ax[1].errorbar([q+istart],[ratio],yerr=([ratioErr],[ratioErr]),fmt='o',\
        #            markersize=msize,fillstyle=fstyle,color=colx)
        ax[1].errorbar([q+istart],[1000*diff],yerr=([1000*diffErr],[1000*diffErr]),fmt='o',\
                    markersize=msize,fillstyle=fstyle,color=colx)
        ax[2].plot([q+istart],[rchi2[ix[q]]],'o',\
                    markersize=msize,fillstyle=fstyle,color=colx)
        print ("{},{},{},{},{},{}".format(target,j,inst[ix[q]],ff[ix[q]],diff*1000,diffErr*1000), file=fout)
    ylim = ax[0].get_ylim()
    ax[0].text(istart,ymin,j)
    ax[0].fill([istart-1,istart+3,istart+3,istart-1], [ymax-0.08,ymax-0.08,ymax+0.08,ymax+0.08], 'grey', alpha=0.2, edgecolor='r')
    ax[1].text(istart,ratmin,j)
    ax[1].fill([istart-1,istart+3,istart+3,istart-1], [ratmin+5,ratmin+5,ratmax-5,ratmax-5], 'grey', alpha=0.2, edgecolor='r')
    ax[2].text(istart,0.75,j)
    ax[2].fill([istart-1,istart+3,istart+3,istart-1], [0.75,0.75,1.75,1.75], 'grey', alpha=0.2, edgecolor='r')
    istart += 10
#
ax[0].set_ylabel(r"Fe Line Energy (keV)")
ax[0].get_xaxis().set_visible(True)
ax[0].set_ylim((ymax-0.1,ymax+0.1))
ax[0].set_title("{} Analysis".format(target.capitalize()))
ax[0].grid(True)
#
# custom legend
#
legend_elements = [Line2D([0],[0],marker='o',color='w',markerfacecolor='b',markersize=10,label='MOS1'),\
                   Line2D([0],[0],marker='o',color='w',markerfacecolor='g',markersize=10,label='MOS2'),\
                   Line2D([0],[0],marker='o',color='w',markerfacecolor='r',markersize=10,label='PN')]
ax[0].legend(handles=legend_elements,loc=2)
#
# now middle panel, the difference
#
ax[1].get_xaxis().set_visible(True)
ax[1].set_ylim((ratmin,ratmax))
ax[1].set_ylabel(r"Difference to Fe K$_\alpha$ (eV)")
ax[1].grid(True)
ax[1].tick_params(labelbottom='off')
y_major_ticks = np.arange(-80, 80.0, 20)
ax[1].set_yticks(y_major_ticks)
ax[1].axhline(color='k')
#
# the bottom panel, reduced chi2
#
ax[2].set_xlabel("Revolution index")
ax[2].set_ylabel(r"Reduced $\chi^2$")
ax[2].set_ylim((0.7,2.0))
ax[2].grid(True)
ax[2].tick_params(labelbottom='off')   
fig.subplots_adjust(hspace=0)
plt.savefig('{}/{}_fe_fit_results.png'.format(wdir,target),dpi=100)
plt.show()
plt.close()
#
fout.close()
