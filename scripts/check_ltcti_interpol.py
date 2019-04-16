#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  7 11:38:52 2018

Long term CTI check the interpolation at 5.9 and 6.4

@author: ivaltchanov
"""
import numpy as np
from astropy.io import fits

import matplotlib.pylab as plt

#import matplotlib.style
import matplotlib as mpl
mpl.style.use('classic')

#%%
wdir = '/home/ivaltchanov/IVAN/PN_SW'
hdu = fits.open(f"{wdir}/ccfdev/EPN_CTI_0049.CCF")
tt = hdu['LTC_TIMES'].data['TIME'][0]
xx = hdu['LONG_TERM_CTI'].data
sel = np.where(xx['MODE_ID'] == 3)[0]
xsel = xx[sel]
e1 = xsel['ENERGY'][0]
e2 = xsel['ENERGY'][1]
e3 = xsel['ENERGY'][2]
tcoeff1 = xsel['T_COEFF'][0]
tcoeff2 = xsel['T_COEFF'][1]
tcoeff3 = xsel['T_COEFF'][2]
# now calculate the correction
corr1 = np.power((1.0 - tcoeff1)/(1.0 - tcoeff1[0]),190.0)
corr2 = np.power((1.0 - tcoeff2)/(1.0 - tcoeff2[0]),190.0)
corr3 = np.power((1.0 - tcoeff3)/(1.0 - tcoeff3[0]),190.0)
#
# now linearly interpolate to get the expected value at 6.2 keV
#
nt = len(tt)
xp = [e1,e2,e3]
eq = 6.3 # keV
outy = np.zeros(nt)
for i in np.arange(nt):
  yp = [corr1[i],corr2[i],corr3[i]]  
  outy[i] = np.interp(eq,xp,yp)
#
fig = plt.figure(figsize=(12,8))
ax = fig.subplots(nrows=2,ncols=1,sharex=True)
#ax.plot(tt,tcoeff1)
ax[0].plot(tt,corr2,label=r'E$_{lab}$=%.4f keV'%(e2))
ax[0].plot(tt,corr3,label=r'E$_{lab}$=%.4f keV'%(e3))
ax[0].plot(tt,outy,label=r'E$_{lab}$=%.4f keV'%(eq))
#ax[0].xlabel("Time since 2000-01-01:T00:00:00")
ax[0].set_ylabel(r"E$_{obs}$/E$_{lab}$")
ax[0].legend()
ax[0].grid(True)
#
yy3 = e3 - e3/corr3
yy4 = eq - eq/outy
ax[1].plot(tt,1000*(yy3-yy4))
#ax.plot(tt,,label=r'E$_{lab}$=%.4f keV'%(eq))
ax[1].set_xlabel("Time since 2000-01-01:T00:00:00")
ax[1].set_ylabel("E(6.4 keV) - E(6.2 keV) (eV)")
ax[1].grid(True)
plt.tight_layout()
plt.show()
#plt.savefig(f'{wdir}/PN_SW_corrections_work.png',dpi=100)
#plt.close()


