#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  5 11:50:14 2018

@author: ivaltchanov
"""
import os
import shutil
import subprocess
import sys 

import numpy as np

from astropy.table import Table

wdir = '/home/ivaltchanov/IVAN/PN_SW'

os.chdir(wdir)
#
t = Table.read(f'{wdir}/selected_targets_nocti.csv')
nt = len(t)

imgdir = f'{wdir}/latex/images'
#
for i in np.arange(nt):
    tgt = t['target'][i]
    obsid = '{:010}'.format(t['obsid'][i])
    rev = t['rev'][i]
    ps_file = f'{wdir}/{tgt}/spec_fit/{obsid}_pn_fit_cti49.ps'
    out_file = f'{imgdir}/{tgt}_{obsid}_pn_fit_cti49.ps'
    png_file = f'{imgdir}/{tgt}_{obsid}_pn_fit_cti49.png'
    if (os.path.isfile(ps_file)):
        shutil.copy(ps_file,out_file)
        runcmd = f'convert -density 300 {out_file} {png_file}'
        try:
            result = subprocess.run([runcmd], shell=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
            retcode = result.returncode
            if retcode < 0:
                print("convert was terminated by signal {} \n {}".format(-retcode,result.stdout.decode()))
            else:
                print("convert returned {} \n {}".format(retcode,result.stdout.decode()))
        except OSError as e:
            print("Execution of convert failed:", e, file=sys.stderr)
    else:
        print (f"File {ps_file} not found")
    pass
#