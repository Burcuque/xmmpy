#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 18 16:38:44 2018

Plot the extracted full spectra for all observations for a target

@author: ivaltchanov
"""

import os
import xspec

import argparse
import glob

import numpy as np 

#from astropy.table import Table
from astropy.io import fits
#from astropy.stats import median_absolute_deviation as mad
#from astropy.stats import mad_std

#import matplotlib as mpl
#mpl.use('Agg')

import matplotlib.pylab as plt

#import logging
#
#%%
def savitzky_golay(y, window_size, order, deriv=0, rate=1):
    r"""Smooth (and optionally differentiate) data with a Savitzky-Golay filter.
    The Savitzky-Golay filter removes high frequency noise from data.
    It has the advantage of preserving the original shape and
    features of the signal better than other types of filtering
    approaches, such as moving averages techniques.
    Parameters
    ----------
    y : array_like, shape (N,)
        the values of the time history of the signal.
    window_size : int
        the length of the window. Must be an odd integer number.
    order : int
        the order of the polynomial used in the filtering.
        Must be less then `window_size` - 1.
    deriv: int
        the order of the derivative to compute (default = 0 means only smoothing)
    Returns
    -------
    ys : ndarray, shape (N)
        the smoothed signal (or it's n-th derivative).
    Notes
    -----
    The Savitzky-Golay is a type of low-pass filter, particularly
    suited for smoothing noisy data. The main idea behind this
    approach is to make for each point a least-square fit with a
    polynomial of high order over a odd-sized window centered at
    the point.
    Examples
    --------
    t = np.linspace(-4, 4, 500)
    y = np.exp( -t**2 ) + np.random.normal(0, 0.05, t.shape)
    ysg = savitzky_golay(y, window_size=31, order=4)
    import matplotlib.pyplot as plt
    plt.plot(t, y, label='Noisy signal')
    plt.plot(t, np.exp(-t**2), 'k', lw=1.5, label='Original signal')
    plt.plot(t, ysg, 'r', label='Filtered signal')
    plt.legend()
    plt.show()
    References
    ----------
    .. [1] A. Savitzky, M. J. E. Golay, Smoothing and Differentiation of
       Data by Simplified Least Squares Procedures. Analytical
       Chemistry, 1964, 36 (8), pp 1627-1639.
    .. [2] Numerical Recipes 3rd Edition: The Art of Scientific Computing
       W.H. Press, S.A. Teukolsky, W.T. Vetterling, B.P. Flannery
       Cambridge University Press ISBN-13: 9780521880688
    """
    import numpy as np
    from math import factorial
    
    try:
        window_size = np.abs(np.int(window_size))
        order = np.abs(np.int(order))
    except ValueError:
        raise ValueError("window_size and order have to be of type int")
    if window_size % 2 != 1 or window_size < 1:
        raise TypeError("window_size size must be a positive odd number")
    if window_size < order + 2:
        raise TypeError("window_size is too small for the polynomials order")
    order_range = range(order+1)
    half_window = (window_size -1) // 2
    # precompute coefficients
    b = np.mat([[k**i for i in order_range] for k in range(-half_window, half_window+1)])
    m = np.linalg.pinv(b).A[deriv] * rate**deriv * factorial(deriv)
    # pad the signal at the extremes with
    # values taken from the signal itself
    firstvals = y[0] - np.abs( y[1:half_window+1][::-1] - y[0] )
    lastvals = y[-1] + np.abs(y[-half_window-1:-1][::-1] - y[-1])
    y = np.concatenate((firstvals, y, lastvals))
    return np.convolve( m[::-1], y, mode='valid')
#%%
os.environ['SAS_CCFPATH'] = '/xdata/ccf/pub'
#%%
#
# get the arguments
#parser = argparse.ArgumentParser(description='Read, fit and plot spectra, saving the results')
#parser.add_argument('obsid', type=str,
#                    help='The OBSID to process')
#parser.add_argument('inst', type=str,
#                    help='Instrument: mos1, mos2 or pn')
#parser.add_argument('--line', type=float, default=6.7,
#                    help='Initial line enery in keV')
#parser.add_argument('--output', type=str, default="output_xspec.csv",
#                    help='Output file, will be appended')
#parser.add_argument('--wdir', type=str, default=os.getcwd(),
#                    help='Where the ODF files are')
#
#args = parser.parse_args()
#args = parser.parse_args("0555471101 pn --line 6.7 --output output_test.csv --wdir /home/ivaltchanov/XMM/CAL/PN_SW/WR140")
#
#%%
#
#
#
wdir = '/home/ivaltchanov/XMM/CAL/PN_SW/ngc5506'
#wdir = '/home/ivaltchanov/XMM/CAL/PN_SW/ngc4151'
#wdir = '/home/ivaltchanov/XMM/CAL/PN_SW/ngc4593'
#wdir = '/home/ivaltchanov/XMM/CAL/PN_SW/ngc3783'
#wdir = '/home/ivaltchanov/XMM/CAL/PN_SW/WR140'
#wdir = '/home/ivaltchanov/XMM/CAL/PN_SW/mrk1048'
#wdir = '/home/ivaltchanov/XMM/CAL/PN_SW/ngc3227'
os.chdir(wdir)

for inst in ["mos1","mos2","pn"]:
    zspec = glob.glob("{}/**/xmmsas_20180620_1732-17.0.0/{}_spectrum_grp0.fits".format(wdir,inst),recursive=True)
    if (len(zspec) > 0):
        fig = plt.figure(figsize=(10,6))
        ax = fig.subplots()
        #ax2 = plt.axes([0.6,0.6,0.25,0.25])
        for j in zspec:
            jpath = os.path.dirname(j)
            os.chdir(jpath)
            hdu = fits.open(j)
            target = hdu[0].header['OBJECT']
            rev = hdu[0].header['REVOLUT']
            submode = hdu[0].header['SUBMODE']
            xfilt = hdu[0].header['FILTER']
            #ontime = hdu[1].header['EXPOSURE']/1000.0 # in ksec
            #
            #title = "{} for {} ({}_{})".format(inst,target,obsid,rev)
            try:
                s = xspec.Spectrum(j)
            except:
                print ("Cannot load {} in XSPEC, will skip".format(j))
                xspec.AllData.clear()
                continue
            #
            s.notice("all")
            xspec.Plot.device = '/null'
            xspec.Plot.setRebin(minSig=40.0,maxBins=10)
            xspec.Plot.xAxis = "keV"
            #xspec.Plot.background = True
            xspec.Plot('data')
            # Get coordinates from plot:
            xx = np.asarray(xspec.Plot.x())
            yy = np.asarray(xspec.Plot.y())
            #ysm = savitzky_golay(yy, 21, 7)
            #bkg = np.asarray(xspec.Plot.backgroundVals())
            yerr = xspec.Plot.yErr()
            #
            # now plot the full spectrum in [0.2,12] keV
            #
            ax.set_xscale('log')
            ax.set_yscale('log')
            ax.errorbar(xx,yy,yerr=yerr,fmt='.',color='b',linestyle='solid')
            # 
            # now an inset
            #
            #ax2.loglog(xx,yy,color='b',linestyle='solid')
            #ax.plot(xx,bkg,'+r')
            xspec.AllData.clear()
            #xspec.Xset.closeLog()
            #
            #break
        ax.set_xlim((0.2,12.0))
        ymax = 0.7
        if (inst == 'pn'):
            ax.set_ylim((0.01,ymax*3))
        else:
            ax.set_ylim((0.01,ymax))
            
        #ax2.set_xlim((5.0,8.0))
        #if (inst == 'pn'):
        #    ax2.set_ylim((0.01,ymax*10))
        #else:
        #    ax2.set_ylim((0.01,0.2))
        ax.set_xlabel("Energy (keV)")
        ax.set_ylabel("Normalized cts/s/keV")
        ax.set_title("{} spectra for {}".format(inst,target))
        #ax.grid(True)
        ax.xaxis.grid(True, which='minor')
        ax.yaxis.grid(True, which='minor')
        ax.fill([5.0,8.0,8.0,5.0], [0.01,0.01,10.0,10.0], 'grey', alpha=0.2, edgecolor='r')
        #ay.errorbar(xx,yy,yerr=yerr,fmt='.',linestyle='solid')
        pngfile = "{}/spec_fit/{}_full_spec.png".format(wdir,inst)
        plt.savefig(pngfile,dpi=100)
        plt.close()
        #plt.show()
    else:
        print ("No {} spectra found".format(inst))
    print ("{} doen".format(inst))
print("all done")