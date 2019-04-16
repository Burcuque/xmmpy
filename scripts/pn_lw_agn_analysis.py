#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Plot the line fit results, the difference is with respect to the 
redshifted Fe Kalpha line

PN LW mode, AGNs

@author: ivaltchanov
"""
import os
import glob

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
#target = "WR140"
#targets = ["ngc4593","ngc5506","ngc5548", "ngc4151","ngc3516","ngc2992","ngc3227",'ngc3783','ngc1566']
targets = ["ngc3227","mkn1040","mkn915","ngc6814",'mkn883']
#
# redshifts
#
redshift = {'ngc4151': 0.003262, 'ngc3227': 0.00386, 'mrk1048': 0.0427, 'ngc3783': 0.009755,\
            'ngc4593': 0.008344, 'ngc5506': 0.00589, 'mcg-5-23-16': 0.008226, 'ngc3516': 0.008816,\
            'ngc5548': 0.01627, 'ngc2992':  0.007296, 'ngc1566': 0.005036, 'iras09149': 0.057150,\
            "iras05078": 0.017879, 'ngc7213': 0.005869,"mkn915": 0.024043,"mkn1040":0.016338,\
            "ngc6814": 0.005227, 'mkn883': 0.03787}
feK = 6.399 # the weitghed mean of Kalpha_2 at 6.3908 (intensity 50) and Kalpha_1 at 6.40308 (intensity 100)

#
wdir = "/xdata/xcaldata/XMM/IVAN/PN_LW/sources"
os.chdir(wdir)

#%%
#
all_diffs = []
xrtime = []
msize=10
fig = plt.figure(figsize=(12,12))
(ax,ax2) = fig.subplots(2,1)
for target in targets:
    print (f"Adding {target}")
    #out_tab = f'{target}/{target}_output_xspec_cti49_t11.csv'
    out_tab = f'{target}/output_xspec.csv'
    if (not os.path.isfile(out_tab)):
        print (f"No CTIv48 results found for target {target}: {out_tab}")
        continue
    t = Table.read(out_tab,data_start=0,names=("obsid","expo","rev","delta_time","submode",\
                                  "xfilt","inst","ontime","ccd",\
                                  "lineE","lineE_err","cstat","chi2r","dof"))
    nt = len(t)
    #
    lineX =  feK/(1.0 + redshift[target]) # redshifted line position
    rev = t['rev'].data
    rtime = t['delta_time'].data
    xrtime = np.concatenate((xrtime,rtime))
    inst = t['inst'].data
    line = t['lineE'].data
    lineErr = t['lineE_err'].data
    rchi2 = t["chi2r"].data
    diff = (line - lineX)*1000.0
    errs = [lineErr*1000,lineErr*1000.0]
    symb = 'o'
    ax.errorbar(rtime,diff,yerr=(errs),fmt=symb,\
               markersize=msize,label=target.capitalize())
    ax2.errorbar(rev,diff,yerr=(errs),fmt=symb,\
               markersize=msize,label=target.capitalize())
#
ax.axhline(0.0,color='k',ls='dashed')
#ax.set_ylabel(r"Fe Line Energy (keV)")
ax.set_ylabel(r"$\Delta$ E (obs-lab) (eV)")
ax.set_xlabel(r"Years since 2000-01-01T00:00:00")
#ax.set_title("{} Analysis".format(target.capitalize()))
ax.grid(True)
#
# the calclosed data
#
ccdir = "/xdata/xcaldata/XMM/IVAN/PN_LW/CalClosed"
cc_files = glob.glob(f"{ccdir}/*_xspec_single_results.csv")
mnka_lab = 5.8988
cclab = r"Mn K$\alpha$"
for jj in cc_files:
    tcc = Table.read(jj)
    rtime_cc = tcc['delta_time'].data
    rev_cc = tcc['rev'].data
    mnka = tcc['line2'].data
    mnka_err = 1000*tcc['line2err'].data
    cc_diff = 1000*(mnka - mnka_lab)
    ix = np.where(tcc['ccd'] == 4.2)[0]
    ax.errorbar(rtime_cc[ix],cc_diff[ix],yerr=(mnka_err[ix]),fmt="^",\
               markersize=msize,color='black',label=cclab)
    ax2.errorbar(rev_cc[ix],cc_diff[ix],yerr=(mnka_err[ix]),fmt="^",\
               markersize=msize,color='black',label=cclab)
    cclab = ""
    
ax2.axhline(0.0,color='k',ls='dashed')
#ax.set_ylabel(r"Fe Line Energy (keV)")
ax2.set_ylabel(r"$\Delta$ E (obs-lab) (eV)")
ax2.set_xlabel(r"Revolution")
ax2.grid(True)
ax.set_ylim([-40.0,60.0])
ax2.set_ylim([-40.0,60.0])
ax.set_xlim([0.0,20.0])
ax2.set_xlim([0.0,4000.0])
#ax.set_title(f"PN SmallWindow Energy Scale Analysis\n Test #{ver}")
ax.set_title(f"PN LargeWindow Energy Scale Analysis")
plt.legend(numpoints=1)
#plt.savefig('{}/pn_sw_cti49_results_control.png'.format(wdir),dpi=100)
#plt.savefig(f'{wdir}/pn_lw_cti48_results.png',dpi=100)
plt.show()
plt.close()
#