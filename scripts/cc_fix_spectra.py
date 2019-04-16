#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Generate specgroup spectra for calClosed data, prepared for use in XSPEC

@author: ivaltchanov
"""

import os
import subprocess
import sys
from astropy.io import fits

import logging
import argparse
#
# global folders
#
os.environ['SAS_CCFPATH'] = '/xdata/ccf/pub'
#%%
#
parser = argparse.ArgumentParser(description='XMMSAS workflow for CalClosed spectra')
parser.add_argument('obsid', type=str,
                    help='The OBSID to process')
parser.add_argument('subfolder', type=str, default="no_cti",
                    help='The variant \"no_cti\" or \"cti49\"')
parser.add_argument('--odfdir', type=str, default=os.getcwd(),
                    help='Where the ODF files are')
args = parser.parse_args()
#
obsid = args.obsid
inst  = 'pn'
varx = args.subfolder
odfdir = os.path.join(os.path.abspath(args.odfdir),obsid)
skip = False
#
#%%
if (not os.path.isdir(odfdir)):
    print ("The ODF folder {} does not exist! Cannot continue".format(odfdir))
    raise FileNotFoundError
#
ppsDir = os.path.join(odfdir,varx)
#
if (not os.path.isdir(ppsDir)):
    print ("PPS folder {} does not exist! Cannot continue".format(ppsDir))
    raise FileNotFoundError
#
os.environ["SAS_ODF"] = odfdir
ccffile = os.path.join(odfdir,"ccf.cif")
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
#%%
outfile = 'pn_spectrum_grp0.fits'
if (os.path.isfile(outfile) and skip):
    exit
spec = 'spectrum_bin5.fits'
hdu = fits.open(spec)
hdu[1].header['BACKFILE'] = 'NONE'
hdu.writeto(spec,overwrite=True)
hdu.close()
#%%
#
# Wil use already available canned RMF and ARF
#
rmfset = '/home/ivaltchanov/IVAN/PN_SW/epn_e3_sw20_sY9_v17.0.rmf'
arfset = '/home/ivaltchanov/IVAN/PN_SW/cc_ancillary.ds'
#
command = "specgroup spectrumset={} addfilenames=yes rmfset={}".format(spec,rmfset) + \
    " arfset={} groupedset={}".format(arfset,outfile)
print ("*** Link bkg, ARF and RMF filenames to the {} spectrum".format(inst))
try:
    result = subprocess.run(command, shell=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
    retcode=result.returncode
    if retcode < 0:
        print("specgroup was terminated by signal", -retcode, file=sys.stderr)
        logging.warning("specgroup was terminated by signal: {}".format(-retcode))
    else:
        print("specgroup returned", retcode, file=sys.stderr)
        logging.info("specgroup returned {}, \n {}".format(retcode,result.stdout.decode()))
except OSError as e:
    print("Execution of specgroup failed:", e, file=sys.stderr)
    logging.error("Execution of specgroup failed: {}".format(e))

print ("All done")