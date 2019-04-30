#!/usr/bin/env python
# coding: utf-8

# # XCLASS: creating images for targets
# 
# For Matej Kosiba's neural network training. Workflow to download and save images. There will be 3 images per target coordinates: one optical from DSS22 red band, and two XMM-Newton for two bands.
# 
# 1. Set the coordinates of the target.
# 1. Identify the observation id (`OBS_ID`) where the target falls in, using `astroquery.esasky` module.
# 1. Get XMM-Newton pipeline (PPS) images for the PN in all 4 bands that contain the target. Using NXSA URL access.
# 1. Combine bands 2+3 and 4+5 to make two X-ray images in the soft band [0.5-2.0] keV and hard band [2.0-12.0] keV.
# 1. Get the DSS2 Red band image withint the input box size from SkyView using `astroquery.skyview` module.
# 1. Save the optical image.
# 1. Crop the X-ray images in the user selected box (usually the same as the optical).
# 1. Save th ecropped X-ray images
# 1. In parallel, display the three images in the notebook.
# 
# **Note:** the optical and the X-ray cropped images will be saved without any re-normalisation.
# 
import os
import tarfile
import requests

import numpy as np

from astropy.io import fits
from astropy.wcs import WCS
from astropy.visualization import PercentileInterval, ImageNormalize
from astropy.nddata import Cutout2D

from astropy.coordinates import SkyCoord
from astropy import units as u

from astroquery.skyview import SkyView

import pandas as pd
from io import StringIO

import matplotlib.pylab as plt

#%%
def get_xmm_obsid(coord,search_radius=15.0,verbose=False):
    #
    # For an input SkyCoord object coord, find all XMM-Newton observations that 
    # intersect with a cricle with search_radius (default 15 arcmin)
    #
    # Uses TAP access to NXSA
    #
    # returns the OBS_ID
    tap_cat = "v_public_observations"
    rad_deg = search_radius/60.0
    circle = f"circle('ICRS',{coord.ra.value:.4f},{coord.dec.value:.4f},{rad_deg})"
    # 
    if (verbose):
        print (f"Search in {tap_cat} for observations intersecting with {circle}")
    #
    tap_url = "http://nxsa.esac.esa.int/tap-server/tap/sync?REQUEST=doQuery&LANG=ADQL&FORMAT=csv&"
    # now combine two tables to include the observing mode
    # instrument_mode_oid for EPN is from 60 to 69, will only search for those
    #
    # the ADQL query:
    #
    query = f"QUERY=select top 100 t2.observation_id,t2.ra,t2.dec,t2.revolution,t2.duration,t3.mode_friendly_name" + \
    " from v_exposure as t1, v_public_observations as t2, v_instrument_mode as t3" + \
    " WHERE " + \
    "(t1.observation_oid = t2.observation_oid) and " + \
    "(t1.instrument_mode_oid BETWEEN 60 and 69) and " + \
    "(t1.instrument_mode_oid = t3.instrument_mode_oid) and " + \
    f"(1=intersects(observation_fov_scircle,{circle}))"
    xreq = tap_url + query
    with requests.get(xreq) as r:
        jx = r.content.decode()
    #
    if ('ERROR' in jx):
        print ("ADQL query returns error:",jx)
        raise Exception
    buff = StringIO(jx)
    df = pd.read_csv(buff, sep=",")
    if (verbose):
        print (f"Found {len(df)} XMM-Newton EPIC-pn public images in the archive")
    #
    # calculate the separation
    #
    cobs = SkyCoord(df.ra,df.dec,frame='icrs',unit='deg')
    sep = cobs.separation(c)
    # add a new column with the separation in arcmin
    df["offset"] = sep.to(u.arcmin)
    obs_table = df.sort_values("offset")
    obs_table.rename(columns={'observation_id': 'obs_id','mode_friendly_name': 'obs_mode','revolution': 'rev'},inplace=True)
    #
    # now select the closest one, but only if it is with Large Window of Full Fram (or extended Full Frame) mode
    k = 0
    xmode = obs_table.iloc[[k]].obs_mode.values[0]
    while (not ("Large" in xmode or "Full" in xmode)):
        if (k > len(obs_table)):
            print ("Could not find Full Frame or Large Window mode observations.")
            obsid = None
            break
        obsid = f"{obs_table.iloc[[k]].obs_id.values[0]:010}"
        xmode = obs_table.iloc[[k]].obs_mode.values[0]
        k += 1
    #
    if (verbose):
        print (f"*** Will pick up the one where the input target is closer to tha axis:")
        print (obs_table.iloc[[k-1]])    
    return obsid
#%%
def get_xmm_pps_images(obs_id, check_cache=True, tmpdir=".",verbose=False):
    #
    # For a given OBS_ID, download the XMM pipeline images and merge the bands to create
    # a soft-band image in [0.5-2] keV and a hard-band one ni [2-12] keV
    # return both as HDUs
    #
    # Uses URL access to NXSA
    #
    # if check_cache is True then it will check if a tar file with the same OBS_ID already
    # exists in tmpdir and skip downloading it again.
    #
    ftype = "IMAGE"
    extn = "FTZ"
    inst = "PN"
    url = 'http://nxsa.esac.esa.int/nxsa-sl/servlet'
    #req  = f"{url}/data-action-aio?obsno={iobs}&name={ftype}&level=PPS"
    req  = f"{url}/data-action-aio?obsno={obs_id}&extension={extn}&name={ftype}&instname={inst}&level=PPS"
    
    tarFile = f"{tmpdir}/pps_{obs_id}.tar"
    if (check_cache and os.path.isfile(tarFile)):
        print (f"Tar file {tarFile} already exists, will skip downloading it again.")
    else:
        print (f"Tar file {tarFile} does not found: downloading...")        
        with requests.get(req) as r:
            #r.raise_for_status() # ensure we notice bad responses
            if (b'No results' in r.content):
                print ("No PPS products found for {obsid}")
                raise Exception
            else:
                with open(tarFile,"wb") as tmp:
                    tmp.write(r.content)
    if (verbose):
        print (f"{inst} images in all bands saved to {tarFile}")
    # ## Extract the XMM images
    # 
    # Here we extract the images and then co-add the bands to create two images in bands [0.5-2] and [2-12.0] keV. 
    # This means co-adding band 2+3 and 4+5.
    bands = [2,3,4,5]
    maps = {}
    with tarfile.open(tarFile,'r') as tar:
        for xband in bands:
            sband = f'{xband}000'
            for member in tar.getmembers():
                if (sband in member.name):
                    print (f"Extracting {member.name}")
                    f=tar.extract(member,path=tmpdir)
                    hdu = fits.open(f"{tmpdir}/{member.name}")
                    maps[sband] = hdu[0]
    #
    # co-adding 
    #
    sb = maps['2000'].data + maps['3000'].data
    hb = maps['4000'].data + maps['5000'].data
    hdu_soft = fits.PrimaryHDU(sb)
    hdu_soft.header = maps['2000'].header
    hdu_hard = fits.PrimaryHDU(hb)
    hdu_hard.header = maps['4000'].header
    #
    return hdu_soft,hdu_hard

#%%
def get_dss2r_image(coord,out_name="unknown",out_dir="."):
    # Download optical images
    # 
    # Will use DSS2 red images available through the `astroquery.SkyView` service. Can be slow! 
    #
    # returns the filename for the saved DSS2 Red band image
    #
    paths = SkyView.get_images(c,survey=['DSS2 Red'],width=box,height=box, coordinates="ICRS")
    ohdu = paths[0][0]
    if not os.path.isdir(out_dir):
        print ("Please set the out_dir to an aleady existing folder. It is needed.")
        raise FileNotFoundError
    else:
        fout = f"{out_dir}/{out_name}_dss2r.fits"
        ohdu.writeto(fout,overwrite=True)
    return ohdu

#%%
def cutout_image(hdu_in,coord,box,out_name="output.fits",clobber=False):
    #
    # Cut an input image in a smaller one, centred on coord (SkyCoord object)
    # with size box (in arcmin) in both x and y direction.
    #
    # Returns the cutout HDU
    #
    #
    zoomSize = u.Quantity((box,box), u.arcmin)
    wcs = WCS(hdu_in.header)
    cutout = Cutout2D(hdu_in.data, coord, zoomSize, wcs=wcs)
    #
    xdu = fits.PrimaryHDU(cutout.data)
    xdu.header = hdu_in.header
    xdu.header.update(cutout.wcs.to_header())
    xdu.writeto(out_name,overwrite=clobber)
    return xdu

#%%
def plot_xmm_optical(xhdu_soft,xhdu_hard,dss2r):
    #
    # hdu_soft and hdu_hard are the cutout images in both XMM bands
    # optical is the DSS2 red image (already cut out)
    #
    fig = plt.figure(figsize=(10,3),dpi=100)
    #
    pp = 99.5 #
    # labels for the plots
    bands = ["soft [0.5-2.0]","hard [2.0-12.0]"]
    #qbands = ["500_2000","2000_12000"]
    #
    for i,ximage in enumerate([xhdu_soft,xhdu_hard]):
        wcs_cut = WCS(ximage.header)
        ax = fig.add_subplot(1,3,i+1,projection=wcs_cut)
        ax.set_title(f'Band: {bands[i]} keV')
        lon = ax.coords['ra']
        lon.set_axislabel('RA (J2000.0)')
        lon.set_major_formatter('hh:mm:ss.s')
        lat = ax.coords['dec']
        lat.set_major_formatter('dd:mm')
        if (i == 0):
            lat.set_axislabel('Dec (J2000.0)')
        # now normalize the imahe
        norm = ImageNormalize(ximage.data[~np.isnan(ximage.data)], interval=PercentileInterval(pp))
        ax.imshow(ximage.data,norm=norm,cmap=plt.cm.plasma,origin='lower',interpolation='nearest')
        ax.set_autoscale_on(False)
    #
    # now the DSS2 Red image
    pp = 99.9 #
    wcs = WCS(dss2r.header)
    ax = fig.add_subplot(1,3,3,projection=wcs)
    ax.set_title('Band: DSS2 red')
    lon = ax.coords['ra']
    lon.set_axislabel('RA (J2000.0)')
    lon.set_major_formatter('hh:mm:ss.s')
    lat = ax.coords['dec']
    #lat.set_axislabel('Dec (J2000.0)')
    lat.set_major_formatter('dd:mm')
    # now normalize the imahe
    norm = ImageNormalize(dss2r.data[~np.isnan(dss2r.data)], interval=PercentileInterval(pp))
    ax.imshow(dss2r.data,norm=norm,cmap=plt.cm.gray,origin='lower',interpolation='nearest')
    ax.set_autoscale_on(False)
    #
    plt.tight_layout(pad=2)

#%%
home = os.path.expanduser('~')

tmpdir = os.path.join(home,"tmp","XMM_data")
if (not os.path.isdir(tmpdir)):
    print ("Please set the tmpdir to an aleady existing folder. It is needed.")
    raise FileNotFoundError
output_dir = tmpdir = os.path.join(home,"tmp","XMM_data","outputs")
if (not os.path.isdir(output_dir)):
    print ("Please set the output_dir to an aleady existing folder. It is needed.")
    raise FileNotFoundError

#
#GCl in the field of 3c234
#ra = "10h01m17.49s"
#dec= "+28d51m11.4s"
# mrk 1040
#ra = "02h28m14.0s"
#dec= "+31d18m42.0s"
# NGC4151
#ra = "12h10m32.5s"
#dec= "+39d24m20.0s"
ra = 182.6354
dec = 39.4056
sid = "ngc4151" # source ID, arbitrary or from catalogue
#
#c = SkyCoord(ra, dec, frame='icrs')
# if ra,dec are in degrees
c = SkyCoord(ra, dec, frame='icrs',unit='deg')
#
# box size for the region of interest, the images will be cropped to it.
box = 6.0*u.arcmin
#%%

xobsid = get_xmm_obsid(c,search_radius=15.0,verbose=True)

#%%

soft, hard = get_xmm_pps_images(xobsid, check_cache=True, tmpdir=tmpdir,verbose=True)

#%%

xdss2r = get_dss2r_image(c,out_name=f"{sid}",out_dir=output_dir)
#%%
soutfile = f"{output_dir}/{sid}_soft_cutout.fits"
cut_soft = cutout_image(soft,c,box,out_name=soutfile,clobber=True)
houtfile = f"{output_dir}/{sid}_hard_cutout.fits"
cut_hard = cutout_image(hard,c,box,out_name=houtfile,clobber=True)
#
#%%
#
# let me try standartizing the images using sklearn

#from sklearn.preprocessing import scale
#
cut_soft.data = (cut_soft.data - np.mean(cut_soft.data))/np.std(cut_soft.data)
cut_hard.data = (cut_hard.data - np.mean(cut_hard.data))/np.std(cut_hard.data)
xdss2r.data = (xdss2r.data - np.mean(xdss2r.data))/np.std(xdss2r.data)
#%%
#
# save the scaled images 
#
soutfile = f"{output_dir}/{sid}_soft_std.fits"
houtfile = f"{output_dir}/{sid}_hard_std.fits"
ofile = f"{output_dir}/{sid}_dss2r_std.fits"
cut_soft.writeto(soutfile,overwrite=True)
cut_hard.writeto(houtfile,overwrite=True)
xdss2r.writeto(ofile,overwrite=True)

#print (np.min(cut_soft.data),np.max(cut_soft.data))
#print (np.min(cut_hard.data),np.max(cut_hard.data))
#print (np.min(xdss2r.data),np.max(xdss2r.data))

#%%
#
# plot the images just to make sure they are OK.

plot_xmm_optical(cut_soft,cut_hard,xdss2r)
#
print ("All done")






