#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 10 10:43:19 2018

Step 03 of the workflow for EPIC processing an XMM observation

Filter the calibrated event lists for high background and generate 
separate Good TimeIntervals for MOS1, MOS2 and PN

@author: ivaltchanov
"""

import os
import subprocess
import sys
import glob

import argparse

import numpy as np 

from astropy.table import Table
from astropy.io import fits
#from astropy.stats import median_absolute_deviation as mad
from astropy.stats import mad_std

import matplotlib as mpl
mpl.use('Agg')

import matplotlib.pylab as plt

import logging
#
# global folders
#
os.environ['SAS_CCFPATH'] = '/xdata/ccf/pub'
#%%
#
# get the arguments
parser = argparse.ArgumentParser(description='XMM SAS workflow step 03 filtering for GTI separately mos1, mos2 and pn')
parser.add_argument('obsid', type=str,
                    help='The OBSID to process')
parser.add_argument('evlist', type=str,
                    help='The input evlist to filter')                    
parser.add_argument('--odfdir', type=str, default=os.getcwd(),
                    help='Where the ODF files are')
#
args = parser.parse_args()
#
#%%
#
odfdir = os.path.abspath(args.odfdir)
obsid = args.obsid
evlist = args.evlist

if (not os.path.isdir(odfdir)):
    print ("The ODF folder {} does not exist! Cannot continue".format(odfdir))
    raise FileNotFoundError
#
#
# crate output folder for PPS and logging
#
# get the SAS version, will be used for the output folder 
try:
    result = subprocess.run(["sasversion -v"], shell=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
    retcode = result.returncode
    if retcode < 0:
        print("sasversion was terminated by signal {} \n {}".format(-retcode,result.stdout.decode()))
    else:
        print("sasversion returned {} \n {}".format(retcode,result.stdout.decode()))
except OSError as e:
    print("Execution of sasversion failed:", e, file=sys.stderr)
#
sasversion = result.stdout.decode().split('[')[1].split(']')[0]
#
ppsDir = os.path.join(odfdir,sasversion)
if (not os.path.isdir(ppsDir)):
    print ("PPS folder {} does not exist. This means step #01 was not run?".format(ppsDir))
    raise FileNotFoundError
#
logging.basicConfig(level=logging.DEBUG,
                    format='%(asctime)s %(levelname)s %(message)s',
                    filename=os.path.join(ppsDir,'proc_step03_sep.log'),
                    filemode='w')
#
os.environ["SAS_ODF"] = odfdir
ccffile = os.path.join(odfdir,'ccf.cif')
if (not os.path.isfile(ccffile)):
    print ("No ccf.cif file found in ODF folder {}. This means step #01 was not run?".format(odfdir))
    raise FileNotFoundError
#
os.environ["SAS_CCF"] = ccffile
logging.info("Set SAS_ODF to {}".format(odfdir))
logging.info("Set SAS_CCF to {}".format(ccffile))
#
os.chdir(ppsDir)
#
task = "evselect"
#
print ("*** RUNNING %s"%task)
if (not os.path.isfile(evlist)):
    logging.error("No ImagingEvts files found in folder {}".format(ppsDir))
    raise FileNotFoundError
#
iev = os.path.basename(evlist)
print ("Generating lightcurve for {}".format(iev))
if ('EMOS1' in iev):
    rate = 'rate_mos1.fits'
    expr = '#XMMEA_EM && (PI>10000) && (PATTERN==0)'
    inst = 'mos1'
elif ('EMOS2' in iev):
    rate = 'rate_mos2.fits'
    expr = '#XMMEA_EM && (PI>10000) && (PATTERN==0)'
    inst = 'mos2'
elif ('EPN' in iev):
    rate = 'rate_pn.fits'
    expr = ' #XMMEA_EP && (PI>10000&&PI<12000) && (PATTERN==0)'        
    inst = 'pn'
   #    
command = 'evselect table={} withrateset=Y rateset={}'.format(iev,rate) + \
    ' maketimecolumn=Y timebinsize=100 makeratecolumn=Y' + \
    ' expression=\'{}\''.format(expr)
try:
    result = subprocess.run(command, shell=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
    retcode=result.returncode
    if retcode < 0:
        logging.warning("evselect was terminated by signal {} \n {}".format(-retcode,result.stdout.decode()))
        print("evselect was terminated by signal {} \n {}".format(-retcode,result.stdout.decode()))
    else:
        print("evselect returned {}, \n {}".format(retcode,result.stdout.decode()))
        logging.info("evselect returned {}".format(retcode,result.stdout.decode()))
except OSError as e:
    print("Execution of evselect failed: {}".format(e))
    logging.error("Execution of evselect failed: {}".format(e))
#
# read the rate tables first
t = {}
time_min = 1.0e10
# get the min time, to be used for relative time
t = Table.read('rate_{}.fits'.format(inst),hdu=1)
time_min = min(np.min(t['TIME']),time_min)
#%%
#
# find the count-rate limit for filtering
# median + 3*MAD
#
medrate = np.median(t['RATE'])
xmad = mad_std(t['RATE'])
ulimit = medrate + 3*xmad
if ('mos' in inst):
    rate_lim = 0.35
else:
    rate_lim = 0.4
#
use_limit = max(ulimit,rate_lim)
logging.info("{}: actual rate limit to use for GTI filtering: {:.3f}".format(inst,use_limit))
#
#%%
#
# now use the derived count-rate limit for the GTI
#
gti_command = 'tabgtigen table=rate_{}.fits'.format(inst) + \
' expression=\'RATE<={}\' gtiset=gti_{}.fits'.format(use_limit,inst)
try:
    result = subprocess.run(gti_command, shell=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
    retcode=result.returncode
    if retcode < 0:
        print("tabgtigen was terminated by signal", -retcode, file=sys.stderr)
        logging.warning("tabgtigen was terminated by signal: {}".format(-retcode))
    else:
        print("tabgtigen returned", retcode, file=sys.stderr)
        logging.info("tabgtigen returned {}, \n {}".format(retcode,result.stdout.decode()))
except OSError as e:
    print("Execution of tabgtigen failed:", e, file=sys.stderr)
    logging.error("Execution of tabgtigen failed: {}".format(e))
#
# check and log the total ONTIME before filtering and the ONTIME with GTI 
hdu = fits.open(evlist)
print ("{} ontime={:.1f}, livetime={:.1f}".format(iev,hdu[1].header['ONTIME'],hdu[1].header['LIVETIME']))
logging.info ("{} ontime={:.1f}, livetime={:.1f}".format(iev,hdu[1].header['ONTIME'],hdu[1].header['LIVETIME']))
#
hdu = fits.open('gti_{}.fits'.format(inst))
mtime = hdu[1].header['ONTIME']
print ("GTI ontime={:.1f}".format(mtime))
logging.info ("GTI ontime={:.1f}".format(hdu[1].header['ONTIME']))
#
#%%
#
# plot the merged GTI on lightcurve

#
fig = plt.figure()
ax = fig.add_subplot(111)
#
# get the good time intervals
#
start = hdu[1].data["START"]
end = hdu[1].data["STOP"]
ngti = len(start)
#
t = Table.read('rate_{}.fits'.format(inst),hdu=1)
reltime = t['TIME'] - time_min
ax.plot(reltime,t['RATE'],label=inst)
ax.set_yscale("log", nonposy='clip')
if ('mos' in inst):
    rate_lim = 0.35
else:
    rate_lim = 0.4
#
ax.axhline(y=rate_lim,color='magenta',ls='solid',label='SAS')
ax.axhline(y=use_limit,color='k',ls='dashed',label='Adopted'.format(use_limit))
ax.set_ylim((0.01,20.0))
ax.grid(True)
ax.legend(loc=2)
ax.set_ylabel("Count rate (cts/s)")
#
# now the GTI bands
#
for jj in np.arange(ngti):
    xx = (start[jj]-time_min,end[jj]-time_min)
    yy1 = (0.01,0.01)
    yy2 = (20.0,20.0)
    ax.fill_between(xx,yy1,yy2,facecolor='yellow',alpha=0.3)
plt.xlabel("Relative time (s)")
ax.set_title("GTI for {}, ONTIME={} s".format(obsid,mtime))
#
png_file = "{}/{}_{}_rates_gti.png".format(ppsDir,obsid,inst)
plt.savefig(png_file,dpi=100)
#plt.show()
logging.info("Saved GTI limits in figure: {}".format(png_file))
plt.close()

#%%
#
# now filter the event lists with the merged GTI
#
if ('_EMOS1' in evlist):
    expr = "#XMMEA_EM && gti(gti_{}.fits,TIME) && (PI>150) && (PATTERN<=12)".format(inst)
    inst = 'mos1'
elif ('_EMOS2' in evlist):
    expr = "#XMMEA_EM && gti(gti_{}.fits,TIME) && (PI>150) && (PATTERN<=12)".format(inst)
    inst = 'mos2'
elif ('_EPN' in evlist):
    expr = "#XMMEA_EP && gti(gti_{}.fits,TIME) && (PI>150) && (PATTERN<=4) && (FLAG==0)".format(inst)
    inst = 'pn'
else:
    raise Exception('Cannot identify the instrument in {}'.format(iev))
#
ev_command = 'evselect table={} withfilteredset=Y filteredset={}_evlist_clean.fits'.format(iev,inst) +  \
' destruct=Y keepfilteroutput=T expression=\'{}\''.format(expr)
#
print ("Filtering {}".format(inst))
try:
    result = subprocess.run(ev_command, shell=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
    retcode=result.returncode
    if retcode < 0:
        print("evselect was terminated by signal", -retcode, file=sys.stderr)
        logging.warning("evselect was terminated by signal {}".format(-retcode))
    else:
        print("evselect returned", retcode, file=sys.stderr)
        logging.info("evselect returned {}, \n {}".format(retcode,result.stdout.decode()))
except OSError as e:
    print("Execution of evselect failed:", e, file=sys.stderr)
    logging.error("Execution of evselect failed: {}".format(e))

#
#%%
print ("All done")        
logging.info ("All done")
