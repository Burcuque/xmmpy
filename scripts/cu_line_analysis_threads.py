
# coding: utf-8

# # Select EPIC PN observations according to some criteria

# In[1]:


import os
import tarfile
import filetype
import subprocess
import sys
import requests
import logging
import threading
from datetime import datetime

import pandas as pd

from astropy.io import fits
from astropy.convolution import convolve, Box1DKernel

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pylab as plt



import seaborn as sns
sns.set(style="white")

plt.rc('text', usetex=False)
plt.rc('font', family='serif')

home = os.path.expanduser('~')
sys.path.append(home + '/IVAN/python')

from xmm_tools import read_pn_obstable, fit_cu_line, run_command

#%%
def proc_ccd(evlist,specfile,selection):
    #
    #
    command = "evselect " + \
        f"table={evlist} energycolumn='PI' withspectrumset=yes" + \
        f" expression='{selection}'" +\
        " withspecranges=yes specchannelmin=0 specchannelmax=20479" +\
        f" spectrumset={specfile} spectralbinsize=5"
    status = run_command(command)
    if (status != 0):
        raise Exception
    return status

# In[2]:
#
# first load the observations from EPICMON account
#

obs_table = read_pn_obstable()
# 
# convert to pandas dataframe
#
obs_tab = obs_table.to_pandas()
# and select all Large Window Mode and duration more thna 100 ks
#
#
texp_limit = 80000.0 # 
xmode = "Large Window"
twork = obs_tab[(obs_tab.perf_dur >= texp_limit) & (obs_tab['mode'].str.contains(xmode))]
nw = len(twork)
print (f"Found {nw} EPIC-PN observations in {xmode} and exposure greater than {texp_limit}")

# In[6]:

# loop over all OBS_IDs
#
outputDir = home + '/IVAN/Cu-line/evlists'
plotDir = home + '/IVAN/Cu-line/plots'
# relative times
time0 = datetime.strptime("2000-01-01T00:00:00","%Y-%m-%dT%H:%M:%S")
#
# file for output
#
fout = open(f"{outputDir}/fit_results_cuka.csv","w")
#
print ("obsid,rev,delta_time,ontime,ccd,lineE,lineE_err,fwhm,fwhm_err,chi2,df",file=fout)
#
doneObs = []
skipDone = False

for ik in range(nw):
    obsid = "{:010}".format(twork.iloc[ik].ac_obs_id)
    print (f"Doing OBS_ID: {obsid} ({ik+1}/{nw})")
    #
    if ((obsid in doneObs)):
        print (f"Skipping as {obsid} is already done")
        continue
    #
    # check if CCD12 is extracted then this OBSID is done
    #
    spec_file = f"{outputDir}/{obsid}_pn_12_spec5.fits"    
    if (skipDone and (os.path.isfile(spec_file))):
        print (f"Skipping as {obsid} is already done")
        continue
    # exposure ID from the file
    expno = twork.iloc[ik].ac_exp_id
    # now tryboth scheduled and unscheduled exposures
    url = "http://nxsa.esac.esa.int/nxsa-sl/servlet/data-action-aio?" + \
        f"obsno={obsid}&expno={expno}&level=PPS&name=PIEVLI&instname=PN"
    #
    ppsFile = f"{outputDir}/{obsid}_evlist.fits.gz"
    #
    if (not os.path.isfile(ppsFile)):
        print (f"No ODF tar file for OBSID {obsid}, downloading from NXSA")
        #with requests.post(url, auth=('ivaltch', 'XXXX')) as r:
        with requests.get(url) as r:
            if (r.status_code == requests.codes.ok):
                with open(ppsFile,"wb") as tmp:
                    tmp.write(r.content)
            else:
                print (f"Cannot download with this URL: {url}. Skipping.")
                continue
    else:
        print (f"PPS event list file {ppsFile} found, will skip downloading it again")
    #
    # now check if it is a tar file, then pick up the first one
    #
    kind = filetype.guess(ppsFile)
    if ('tar' in kind.extension):
        with tarfile.open(ppsFile,'r') as tar:
            for member in tar.getmembers():
                if ('PIEVLI' in member.name):
                    f=tar.extract(member,path=outputDir)
                    ppsFile = f"{outputDir}/{member.name}"
                    break    
    #        
    #
    # now, extract per CCD the spectrum
    #
    threads = []
    for i in range(12):
        ccd = f"{i+1:02}"
        print (f"   Processing CCD # {ccd}")
        mnk_mask =  f"/xdata/xcaldata/XMM/PN/CTI/mask/refmask_{ccd}_mnk_lmap.fits"
        xselect = f"mask({mnk_mask},1,1,RAWX,RAWY) && (CCDNR == {i+1}) && (PATTERN==0) && (PAT_SEQ==0) && #XMMEA_EP"
        spec_file = f"{outputDir}/{obsid}_pn_{ccd}_spec5.fits"
        #
        thread = threading.Thread(target=proc_ccd,args=(ppsFile,spec_file,xselect))
        threads.append(thread)
        #
    for i in range(12):
        threads[i].start()
    for i in range(12):
        threads[i].join()
    #
    print (f"Extraction of per CCD spectra for {obsid} is done")
    # Now plot
    smooth = 7
    fig, axs = plt.subplots(2,6,sharex=True,sharey=True,figsize=(15,10))
    for j in range(12):
        ccd = f"{j+1:02}"
        #slices = glob.glob(f"{specDir}/pn_{ccd:02}_180_spec5.fits")
        slices = f"{outputDir}/{obsid}_pn_{ccd}_spec5.fits"
        if (not os.path.isfile(slices)):
            raise FileNotFoundError
        #
        hdu = fits.open(slices)
        #
        rev = hdu[0].header["REVOLUT"]
        start_time = hdu[0].header['DATE-OBS']
        stime = datetime.strptime(start_time,"%Y-%m-%dT%H:%M:%S")
        delta_time = (stime-time0).total_seconds()/(365.0*24.0*3600.0) # in years
        #
        spec = hdu['SPECTRUM']
        ontime = spec.header["EXPOSURE"]
        channel = spec.data['CHANNEL']*5.0/1000.0
        counts = spec.data['COUNTS']
        y = convolve(counts, Box1DKernel(smooth))
        # now fit a simple linear + Gauss line model
        #
        (fit_out,fit_res) = fit_cu_line(channel,counts,line_c=8.04)
        output = f"{obsid},{rev},{delta_time:.4f},{ontime:.2f},{ccd},{fit_res}"
        print (output,file=fout)
        #
        yfitted = fit_out.eval(x=channel)
        #print (f"Fitting CCD #{ccd:02}",fit_res)
        k = j
        kj = 0
        if (j >= 6):
            k = j - 6
            kj = 1
        if (ccd == "03"):
            axs[kj,k].set_title(f"OBS_ID: {obsid}")
        axs[kj,k].plot(channel,counts,label=f'CCD #{ccd}')
        axs[kj,k].plot(channel,y,label='')
        axs[kj,k].plot(channel,yfitted,color='red',label='')
        axs[kj,k].set_xlim((7,9))
        axs[kj,k].set_ylim((0,100))
        axs[kj,k].grid(True)
        axs[kj,k].legend()
        #axs[kj,k].set_title(f"CCD #{ccd}")
    #
    plt.subplots_adjust(wspace=0, hspace=0)
    #plt.text(-13,-1,'Energy (keV)',ha='center', va='center')
    #plt.text(-36,10,'Counts',rotation='vertical',ha='center', va='center')
    #plt.tight_layout()
    plt.savefig(f"{plotDir}/{obsid}_CuKa_plot.png",dpi=100)
    #plt.show()
    plt.close()
    doneObs.append(obsid)
#
fout.close()
print ("All done")