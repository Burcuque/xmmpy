
# coding: utf-8

# # Select EPIC PN observations according to some criteria

# In[1]:


import os
import tarfile
import subprocess
import sys
import requests
import logging
import threading

import pandas as pd

from astropy.io import fits
from astropy.convolution import convolve, Box1DKernel

import matplotlib.pylab as plt


import seaborn as sns
sns.set(style="white")

plt.rc('text', usetex=False)
plt.rc('font', family='serif')

home = os.path.expanduser('~')
sys.path.append(home + '/IVAN/python')

from xmm_tools import *


# In[33]:


def fit_cu_line(xin,yin, minMax=(7.0,9.0),line_c=8.0):
    # xin is the energy in keV
    m1 = xin >= minMax[0]
    m2 = xin <= minMax[1]
    xw = xin[m1*m2]
    yw = yin[m1*m2]
    i1max = np.argmax(yw)
    y1max = yw[i1max]
    #x1max = xw[i1max]
    #
    poly_mod = PolynomialModel(1,prefix='poly_')
    pars = poly_mod.guess(yw, x=xw)
    #
    gauss1  = GaussianModel(prefix='g1_')
    pars.update( gauss1.make_params())
    pars['g1_center'].set(line_c,min=line_c-0.1,max=line_c+0.1)
    pars['g1_sigma'].set(0.08,min=0.04,max=0.250)
    pars['g1_amplitude'].set(y1max)
    #
    mod = poly_mod + gauss1
    #init = mod.eval(pars, x=x)
    out = mod.fit(yw, pars, x=xw)
    #
    cen = out.params['g1_center'].value
    cen_err = out.params['g1_center'].stderr
    fwhm = out.params['g1_fwhm'].value
    fwhm_err = out.params['g1_fwhm'].stderr
    #peak = out.params['g1_amplitude'].value/expo1
    #peak_err = out.params['g1_amplitude'].stderr/expo1
    chi2 = out.chisqr
    df = len(xw)
    results  = f"{cen:.3f},{cen_err:.3f},{fwhm:.5f},{fwhm_err:.5f},{chi2:.3f},{df}"
    return (out,results)

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
texp_limit = 50000.0 # 
xmode = "Large Window"
twork = obs_tab[(obs_tab.perf_dur >= texp_limit) & (obs_tab['mode'].str.contains(xmode))]
print ("Found {} EPIC-PN observations in {} and exposure greater than {}".format(len(twork),xmode,texp_limit))


# In[5]:


#
# plot the results using altair
import altair as alt

base = alt.Chart(twork,title=r"EPIC-PN {} mode, t > {} ks".format(xmode,texp_limit/1000.0),width=800,height=500)

points = base.mark_point(filled=True, size=50).encode(
    x=alt.X(
        "rev",
        #scale=alt.Scale(domain=(2500,3500)),
        axis=alt.Axis(title='Revolution')
    ),
    y=alt.Y(
        'perf_dur',
        #scale=alt.Scale(type='log',domain=(1,50000)),
        axis=alt.Axis(title="Performed Duration (s)")
    ),
    color='mode',
    tooltip=['rev', 'ac_obs_id','ac_exp_id','filter','target']
)
points


# In[6]:


#
# now download an OBS_ID from NXSA and check for the Cu Ka line (8.04 keV)
#
# select a random OBS_ID from the table
#

obsid = "{:010}".format(twork.sample(n=1).ac_obs_id.values[0])
#print (f"Selected a random OBS_ID: {obsid}")
url = f"http://nxsa.esac.esa.int/nxsa-sl/servlet/data-action-aio?obsno={obsid}&level=PPS&name=PIEVLI&instname=PN"
print (url)


# In[7]:


outputDir = home + '/IVAN/Cu-line'
ppsFile = f"{outputDir}/{obsid}_evlist.fits.gz"
#
if (not os.path.isfile(ppsFile)):
    print (f"No ODF tar file for OBSID {obsid}, downloading from NXSA")
    #with requests.post(url, auth=('ivaltch', 'XXXX')) as r:
    with requests.get(url) as r:
        r.raise_for_status() # ensure we notice bad responses
        with open(ppsFile,"wb") as tmp:
            tmp.write(r.content)
else:
    print (f"PPS event list file {ppsFile} found, will skip downloading it again")
print ("Ready")


# In[8]:


#
# now, extract per CCD the spectrum
#
#
# now loop over all 12 CCDs in EPN, using threads
#
threads = []
for i in range(12):
    ccd = f"{i+1:02}"
    print (f"Processing CCD # {ccd}")
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
print ("Extraction of per CCD spectra is done")


# In[34]:


# Now plot
smooth = 7
fig, axs = plt.subplots(2,6,sharex=True,sharey=True,figsize=(15,10))
for j in range(12):
    ccd = j+1
    #slices = glob.glob(f"{specDir}/pn_{ccd:02}_180_spec5.fits")
    slices = f"{outputDir}/{obsid}_pn_{ccd:02}_spec5.fits"
    if (not os.path.isfile(slices)):
        raise FileNotFoundError
    hdu = fits.open(slices)
    spec = hdu['SPECTRUM']
    channel = spec.data['CHANNEL']*5.0/1000.0
    counts = spec.data['COUNTS']
    y = convolve(counts, Box1DKernel(smooth))
    # now fit a simple linear + Gauss line model
    #
    (fit_out,fit_res) = fit_cu_line(channel,counts,line_c=8.04)
    #
    yfitted = fit_out.eval(x=channel)
    print (f"Fitting CCD #{ccd:02}",fit_res)
    k = j
    kj = 0
    if (j >= 6):
        k = j - 6
        kj = 1
    axs[kj,k].plot(channel,counts,label=f'CCD #{ccd:02}')
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
plt.show();
#plt.savefig(f"{wdir}/{obsid}_{expo}_allccd_plot.png",dpi=100)
#plt.close();

