#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 20 12:27:08 2018

Plot the line fit results, the difference is with respect to the 
redshifted Fe Kalpha line

Version with the current CTI CCF v48

CalClosed fit resutls are from XSPEC

@author: ivaltchanov
"""
import os

import numpy as np

from astropy.table import Table

import matplotlib.pylab as plt
#from matplotlib.patches import Patch
#from matplotlib.lines import Line2D

import seaborn as sns
sns.set(style="white")

plt.rc('text', usetex=False)
plt.rc('font', family='serif')
#
#target = "ngc4593"
#target = "ngc4151"
#target = "mrk1048"
#target = "ngc3227"
#target = "ngc3783"
#target = "ngc5506"
#target = "mcg-5-23-16"
#target = "ngc5548"
#target = "ngc3516"
#target = "WR140"testing run with
#targets = ["ngc4593","ngc5506","ngc5548", "ngc4151","ngc3516","ngc2992","ngc3227",'ngc3783','ngc1566']
targets = ["ngc2992","ngc3227","ngc3783","ngc4151","ngc4593",'ngc5548','ngc7213']
#
# redshifts
#
redshift = {'ngc4151': 0.003262, 'ngc3227': 0.00386, 'mrk1048': 0.0427, 'ngc3783': 0.009755,\
            'ngc4593': 0.008344, 'ngc5506': 0.00589, 'mcg-5-23-16': 0.008226, 'ngc3516': 0.008816,\
            'ngc5548': 0.01627, 'ngc2992':  0.007296, 'ngc1566': 0.005036, 'iras09149': 0.057150,\
            "iras05078": 0.017879, 'ngc7213': 0.005869}
feK = 6.399 # the weitghed mean of Kalpha_2 at 6.3908 (intensity 50) and Kalpha_1 at 6.40308 (intensity 100)

control = ["ngc1566","iras09149","ngc5506","iras05078","ngc7213"]

#if (target == 'WR140'):
#    lineX = 6.697 # see Sugawara et al. (2015)
#
wdir = "/xdata/xcaldata/XMM/IVAN/PN_SW/sources"
#wdir = "/home/ivaltchanov/XMM/CAL/PN_SW/{}".format(target)
#wdir = "/home/ivaltchanov/XMM/CAL/PN_SW/WR140"
os.chdir(wdir)

#%%
#
# version to use
#
all_diffs = []
xrtime = []
msize=10
fig = plt.figure(figsize=(12,12))
(ax,ax2) = fig.subplots(2,1)
for target in targets:
    print (f"Adding {target}")
    #out_tab = f'{target}/{target}_output_xspec_cti49_t11.csv'
    out_tab = f'{target}/{target}_output_xspec_cti48.csv'
    if (not os.path.isfile(out_tab)):
        print (f"No CTIv48 results found for target {target}: {out_tab}")
        continue
    t = Table.read(out_tab,data_start=0,names=("obsid","rev","delta_time","submode",\
                                  "full","xfilt","inst","ontime",\
                                  "lineE","lineE_low","lineE_up","cstat","chi2r","dof"))
    nt = len(t)
    #
    lineX =  feK/(1.0 + redshift[target]) # redshifted line position
    rev = t['rev'].data
    rtime = t['delta_time'].data
    xrtime = np.concatenate((xrtime,rtime))
    ff = t['full'].data
    inst = t['inst'].data
    line = t['lineE'].data
    lineErrUp = t['lineE_low'].data
    lineErrLo = t['lineE_up'].data
    rchi2 = t["chi2r"].data
    ymax = lineX
    #ymax = np.max(line + lineErrUp)
    ymin = np.min(line - lineErrLo)
    #ymin = ymax - 0.
    #
    #
    #ax.errorbar(rtime,line,yerr=([lineErrLo,lineErrUp]),fmt='o',\
    #           markersize=msize)
    #ax.axhline(lineX,color='k',ls='dashed')
    diff = (line - lineX)*1000.0
    all_diffs = np.concatenate((all_diffs,diff))
    errs = [lineErrLo*1000,lineErrUp*1000.0]
    if (target in control):
        symb='*'
    else:
        symb = 'o'
    ax.errorbar(rtime,diff,yerr=(errs),fmt=symb,\
               markersize=msize,label=target.capitalize())
    ax2.errorbar(rev,diff,yerr=(errs),fmt=symb,\
               markersize=msize,label=target.capitalize())
#
#
# mean and stdev of difference
#
#
ax.axhline(1.0,color='k',ls='dashed')
#ax.set_ylabel(r"Fe Line Energy (keV)")
ax.set_ylabel(r"$\Delta$ E (obs-lab) (eV)")
ax.set_xlabel(r"Years since 2000-01-01T00:00:00")
#ax.set_title("{} Analysis".format(target.capitalize()))
ax.grid(True)
#
ax2.axhline(1.0,color='k',ls='dashed')
#ax.set_ylabel(r"Fe Line Energy (keV)")
ax2.set_ylabel(r"$\Delta$ E (obs-lab) (eV)")
ax2.set_xlabel(r"Revolution")
ax2.grid(True)
#
# Now add the calclosed data, only Mn Kalpha double results
#
#t = Table.read("/xdata/xcaldata/XMM/IVAN/PN_SW/CalClosed/pn_sw_calclosed_fit_results_cti49z.csv")
#tcc = Table.read(f"/xdata/xcaldata/XMM/IVAN/PN_SW/CalClosed/pn_sw_calclosed_fit_results_cti49_{ver}.csv")
tcc = Table.read(f"/xdata/xcaldata/XMM/IVAN/PN_SW/CalClosed/cc_output_xspec_cti48_new.csv")
#tcc = Table.read(f"/xdata/xcaldata/XMM/IVAN/PN_SW/CalClosed/cc_output_xspec_cti49_t19.csv")
#tcc = Table.read(f"/xdata/xcaldata/XMM/IVAN/PN_SW/CalClosed/pn_sw_calclosed_fit_results_cti49_t9.csv")
nt = len(tcc)
#
rev = tcc['rev'].data
rtime = tcc['time'].data
#
# AlKa
line1_lab = 1.486 #keV probabilities from Wikipedia on Iron-55
line1 = tcc['al'].data
line1Err = tcc['al_err'].data
cc1_diff = (line1 - line1_lab)*1000.0
cc1_errs = [line1Err*1000,line1Err*1000]
ax.errorbar(rtime,cc1_diff,yerr=(cc1_errs),fmt='h',\
           markersize=msize,color='cyan',fillstyle='full', \
           linestyle='none',alpha=0.5, label=r'AlK$\alpha$')
ax2.errorbar(rev,cc1_diff,yerr=(cc1_errs),fmt='h',\
           markersize=msize,color='cyan',fillstyle='full', \
           linestyle='none',alpha=0.5, label=r'AlK$\alpha$')
# MnKa
line2_lab = (0.162*5.888 + 0.6*5.899)/0.762 #keV probabilities from Wikipedia on Iron-55
line2 = tcc['mn1'].data
line2Err = tcc['mn1_err'].data
cc2_diff = (line2 - line2_lab)*1000.0
cc2_errs = [line2Err*1000,line2Err*1000]
ax.errorbar(rtime,cc2_diff,yerr=(cc2_errs),fmt='^',\
           markersize=msize,color='k',fillstyle='full', \
           linestyle='none',alpha=0.5, label=r'MnK$\alpha$')
ax2.errorbar(rev,cc2_diff,yerr=(cc2_errs),fmt='^',\
           markersize=msize,color='k',fillstyle='full', \
           linestyle='none',alpha=0.5, label=r'MnK$\alpha$')
# MnKb
line3_lab = 6.490 #eV probabilities from Wikipedia on Iron-55
line3 = tcc['mn2'].data
line3Err = tcc['mn2_err'].data
cc3_diff = (line3 - line3_lab)*1000.0
cc3_errs = [line3Err*1000,line3Err*1000]
#ax.errorbar(rtime,cc3_diff,yerr=(cc3_errs),fmt='*',\
#            marker='*',markersize=msize,color='y',fillstyle='full',\
#            linestyle='none', label=r'MnK$\beta$')
#ax2.errorbar(rev,cc3_diff,yerr=(cc3_errs),fmt='*',\
#            marker='*',markersize=msize,color='y',fillstyle='full',\
#            linestyle='none', label=r'MnK$\beta$')

print (cc2_diff[np.where(rtime > 3.0)[0]])
#fig.subplots_adjust(hspace=0)
ax.set_ylim([-100.0,100.0])
ax2.set_ylim([-100.0,100.0])
ax2.set_xlim([0.0,4000.0])
ax.set_title(f"PN SmallWindow Energy Scale Analysis\n EPN_CTI_0048.CCF")
plt.legend(numpoints=1)
plt.savefig('{}/pn_sw_cti48_results.png'.format(wdir),dpi=100)
#plt.savefig(f'{wdir}/pn_sw_cti48_results.png',dpi=100)
plt.show()
plt.close()
#
#%%
# 
# plot as a function of revolution
#
#fout.close()
