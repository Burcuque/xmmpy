#
#
import matplotlib.pylab as plt
#
from astropy.io import ascii
import numpy as np
import os

home = os.path.expanduser('~')
sim_dir = home + "/XMM/simulations/rev0679"

os.chdir(sim_dir)
#import argparse
#
#parser = argparse.ArgumentParser(description='Plot simulation results')
#parser.add_argument('infile', type=str,
#                    help='The CSV file with simulation results')
#args = parser.parse_args()

#t = ascii.read(args.infile,format='fast_no_header')
tt = {}
tt['g1'] = ascii.read('simulation10_results_1gauss.csv',format='fast_csv',data_start=0,\
  names=["label","cen","err1","err2","chi2","dof"])
tt['g2'] = ascii.read('simulation10_results_2gauss.csv',format='fast_csv',data_start=0,\
  names=["label","cen","err1","err2","chi2","dof"])
tt['g4'] = ascii.read('simulation10_results_4gauss.csv',format='fast_csv',data_start=0,\
  names=["label","cen","err1","err2","chi2","dof"])
#%%
inst = []
texp = []
cen = []
err1 = []
err2 = []
xi2 = []
for j in ["g1","g2","g4"]:
    #
    nt = len(tt[j])
    for i in np.arange(nt):
        k = tt[j]['label'][i]
        lab = k.split(".")[0].split('_')
        inst.append(lab[3])
        texp.append(int(lab[2].replace('ks','')))
        cen.append(tt[j]['cen'][i])
        err1.append(tt[j]['cen'][i] - tt[j]['err1'][i])
        err2.append(tt[j]['err2'][i] - tt[j]['cen'][i])
        xi2.append(tt[j]['chi2'][i])
#
inst = np.asarray(inst)
texp = np.asarray(texp)
cen = np.asarray(cen)
err1 = np.asarray(err1)
err2 = np.asarray(err2)
xi2 = np.asarray(xi2)
#
i1 = np.where(inst == 'mos1')[0]
i2 = np.where(inst == 'mos2')[0]
i3 = np.where(inst == 'pn')[0]
#
fig = plt.figure(figsize=(10,8))
ax = fig.subplots()

ax.errorbar(texp[i1],cen[i1],yerr=(err1[i1],err2[i1]),fmt='o',label='MOS1')
ax.errorbar(texp[i2]+2,cen[i2],yerr=(err1[i2],err2[i2]),fmt='o',label='MOS2')
ax.errorbar(texp[i3]+4,cen[i3],yerr=(err1[i3],err2[i3]),fmt='o',label='PN')
ax.hlines(6.41,10,100,linestyle='dashed',label='Input centre')
ax.set_ylim((6.3,6.5))
ax.set_xlabel("Exposure time (ks)")
ax.set_ylabel("Fe K best fit centre (keV)")
ax.set_title("Fe K simulations with 1, 2 or 4 Gaussian models")
ax.grid()
#
plt.legend(loc=4)
#plt.savefig('figure_sims10.png',dpi=100)
plt.show()
#plt.close()
