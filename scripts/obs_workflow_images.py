#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 10 10:43:19 2018

Use the GTI filtered event lists to make spectra and images in different 
energy bands

@author: ivaltchanov
"""

import os
import subprocess
import sys

#import argparse


import logging

import argparse
#
# global folders
#
home = os.path.expanduser('~')
sasversion = 'xmmsas_20180504_1731'
os.environ['SAS_CCFPATH'] = '/xdata/ccf/pub'
#
#%%
#
parser = argparse.ArgumentParser(description='XMMSAS workflow step 3 - images')
parser.add_argument('obsid', type=str,
                    help='The OBSID to process')
parser.add_argument('lowe', type=int,
                    help='Low energy of the band in eV (integer)')
parser.add_argument('hie', type=int,
                    help='High energy of the band in eV (integer)')
parser.add_argument('odfdir', type=str, default=os.getcwd(),
                    help='ODF folder')

#parser.add_argument('task', type=str,default='emproc',
#                    help='Which proc to run: emproc or epproc')
#
args = parser.parse_args()
#
obsid = args.obsid
odfdir = args.odfdir
#
# check if ODF exists, if not fail
if (not os.path.isdir(odfdir)):
    raise FileNotFoundError
os.environ['SAS_ODF'] = odfdir
#
wdir = '{}/{}'.format(odfdir,sasversion)
if (not os.path.isdir(wdir)):
    raise FileNotFoundError
#
logging.basicConfig(level=logging.DEBUG,
                    format='%(asctime)s %(levelname)s %(message)s',
                    filename='{}/images.log'.format(wdir),
                    filemode='w')
#
# check if CCF file exists in ODF dir, if not then fail
ccf_file = "{}/ccf.cif".format(odfdir)
if (not os.path.isfile(ccf_file)):
    logging.error("CCF filenot found in SAS_ODF: {}".format(odfdir))
    raise FileNotFoundError
os.environ['SAS_CCF'] = ccf_file
os.chdir(wdir)
#
#%%
# GTI filtered event lists filenames are hardcoded
# from previous step in workflow: obs_workflow_filtering.py
#
pi0 = args.lowe
pi1 = args.hie
print ("Doing band [{},{}]".format(pi0,pi1))
for inst in ['mos1','mos2','pn']:
    evlist = '{}_evlist_clean.fits'.format(inst)
    if (not os.path.isfile(evlist)):
        logging.error("Cannot find cleaned event lists")
    image_name = '{}_image_{}_{}.fits'.format(inst,pi0,pi1)
    #
    if ('mos' in inst):
        expr = 'PI in [{}:{}] &&  (FLAG & 0x766ba000)==0 && PATTERN in [0:12]'.format(pi0,pi1)
    else:
        expr = 'PI in [{}:{}] &&  FLAG==0 && PATTERN in [0:4]'.format(pi0,pi1)
    #    
    command = 'evselect table={} xcolumn=X ycolumn=Y imagebinning=binSize'.format(evlist) +  \
         ' ximagebinsize=80 yimagebinsize=80' + \
         ' expression=\'{}\''.format(expr) +  \
         ' withimageset=true imageset={}'.format(image_name)
    try:
        result = subprocess.run(command, shell=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
        retcode=result.returncode
        if retcode < 0:
            print("evselect was terminated by signal", -retcode, file=sys.stderr)
            logging.warning("evselect was terminated by signal: {}".format(-retcode))
        else:
            print("evselect returned", retcode, file=sys.stderr)
            logging.info("evselect returned {}, \n {}".format(retcode,result.stdout.decode()))
    except OSError as e:
        print("Execution of evselect failed:", e, file=sys.stderr)
        logging.error("Execution of evselect failed: {}".format(e))
#
logging.info ("All done")        