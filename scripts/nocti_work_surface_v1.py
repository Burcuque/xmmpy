#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 20 12:27:08 2018

Plot the nocti results in a plane: time, Energy for AGNs.

The idea, from Norbert, is to incorporate the energy dependence in the 
derivation of the correction curve.

x = E_obs
y = time
z = E_obs/E_lab/(1+z)

2d-linear surface fit to (x,y,z) => then slice (time,z) at y = 6.4 keV and this should 
be the corrected curve 

@author: ivaltchanov
"""
import os

import numpy as np
from numpy.polynomial import polynomial
import scipy

from astropy.table import Table

import matplotlib.pylab as plt
#from matplotlib.patches import Patch
from matplotlib.lines import Line2D
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 unused import
#%matplotlib

plt.rc('text', usetex=False)
plt.rc('font', family='serif')
#%%
def polyfit2d(x, y, f, deg):
    x = np.asarray(x)
    y = np.asarray(y)
    f = np.asarray(f)
    deg = np.asarray(deg)
    vander = polynomial.polyvander2d(x, y, deg)
    vander = vander.reshape((-1,vander.shape[-1]))
    f = f.reshape((vander.shape[0],))
    c = np.linalg.lstsq(vander, f)[0]
    return c.reshape(deg+1)
##
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
targets = ["ngc2992","ngc3227","ngc3783","ngc4151","ngc4593",'ngc5548']
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
wdir = "/xdata/xcaldata/XMM/IVAN/PN_SW"
#wdir = "/home/ivaltchanov/XMM/CAL/PN_SW/{}".format(target)
#wdir = "/home/ivaltchanov/XMM/CAL/PN_SW/WR140"
os.chdir(wdir)

#%%
xx = []
yy = []
zz = []
zz_err = []
#
for target in targets:
    print (f"Adding {target}")
    tab_file = f'{target}/output_xspec_nocti.csv'
    if (not os.path.isfile(tab_file)):%matplotlib qt
        print (f"No nocti results for {target}")
        continue
    tq = Table.read(tab_file)
    iq = np.where(tq['chi2r'].data <= 1.5)[0]
    if (len(iq) < 1):
        print (f"All fit results for {target} have rChi2 > 1.3, discarding")
        continue
    t = tq[iq]
    nt = len(t)
    #
    lineX =  feK/(1.0 + redshift[target]) # redshifted line position
    rev = t['rev'].data
    rtime = t['delta_time'].data
    ff = t['full'].data
    inst = t['inst'].data
    line = t['lineE'].data
    lineErrUp = t['lineE_low'].data
    lineErrLo = t['lineE_up'].data
    lineErr = (lineErrLo + lineErrUp)/2.0
    rchi2 = t["chi2r"].data
    ratio = line/lineX
    ratio_err = lineErr/lineX
    print ("{}: min line energy {}, max line energy {} ".format(target,np.min(line),np.max(line)))
    xx = np.concatenate((xx,line),axis=None)
    yy = np.concatenate((yy,rtime),axis=None)
    zz = np.concatenate((zz,ratio),axis=None)
    zz_err = np.concatenate((zz_err,ratio_err),axis=None)
    #

#%%
#
#
data = np.c_[xx,yy,zz]

# regular grid covering the domain of the data
mn = np.min(data, axis=0)
mx = np.max(data, axis=0)
X,Y = np.meshgrid(np.linspace(mn[0], mx[0], 20), np.linspace(mn[1], mx[1], 20))
XX = X.flatten()
YY = Y.flatten()
    
# best-fit linear plane (1st-order)
A = np.c_[data[:,0], data[:,1], np.ones(data.shape[0])]
C,_,_,_ = scipy.linalg.lstsq(A, data[:,2])    # coefficients
    
# evaluate it on grid
Z = C[0]*X + C[1]*Y + C[2]

# evaluate on 6.4 keV for t in [0,25,1]
t_run = np.linspace(0.0,25.0,num=26)
zeval = C[0]*6.40 + C[1]*t_run + C[2]
#
print (zeval)
msize=10
fig = plt.figure(figsize=(12,8))
#ax = fig.add_subplot(111,projection='3d')
#ax = fig.add_subplot(111)
#plot = ax.scatter(data[:,0],data[:,1],marker='o',c=zz,cmap=plt.cm.jet,vmin=np.min(zz),vmax=1.0)
#cbar = fig.colorbar(plot)
#cbar.set_label('$E_{obs}/E_{lab}$',fontsize=18)
#plt.imshow(Z,origin='lower')
ax = fig.gca(projection='3d')
ax.plot_surface(X, Y, Z, rstride=1, cstride=1, alpha=0.2)
ax.scatter(data[:,0], data[:,1], data[:,2], c='r', s=50)
#plt.imshow(Z, extent=(xx.min(), yy.max(), xx.max(), yy.min()))
#ax = fig.gca(projection='3d')
#ax.plot_surface(X, Y, Z, rstride=1, cstride=1, alpha=0.2)
#ax.scatter(data[:,0], data[:,1], data[:,2], c='r', s=50)
#ax.set_zlabel('Ratio')
#ax.axis('equal')
#ax.axis('tight')

ax.set_xlabel('Energy (keV)')
ax.set_ylabel('Time')

#ax.set_zlabel('E_obs/E_lab')
plt.grid(True)
plt.title('PN SW without long-term CTI, 3-d view')
plt.show()
#
