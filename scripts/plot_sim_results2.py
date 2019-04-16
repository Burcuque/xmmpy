#
#
import matplotlib.pylab as plt
#
from astropy.io import ascii
import numpy as np
import os

home = os.path.expanduser('~')
#rev = "0679"
rev = "2617"
sim_dir = home + "/XMM/simulations/rev{}".format(rev)

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
#
# now get the min max of spread at each texp for each case
#
xx = np.array([30,50,80,100]) # exposure times used
nx = len(xx)
#
xcen = {'g1_mos1': np.zeros(nx),'g1_mos2': np.zeros(nx),'g1_pn': np.zeros(nx),\
        'g2_mos1': np.zeros(nx),'g2_mos2': np.zeros(nx),'g2_pn': np.zeros(nx),\
        'g4_mos1': np.zeros(nx),'g4_mos2': np.zeros(nx),'g4_pn': np.zeros(nx)}
xcen_up = {'g1_mos1': np.zeros(nx),'g1_mos2': np.zeros(nx),'g1_pn': np.zeros(nx),\
        'g2_mos1': np.zeros(nx),'g2_mos2': np.zeros(nx),'g2_pn': np.zeros(nx),\
        'g4_mos1': np.zeros(nx),'g4_mos2': np.zeros(nx),'g4_pn': np.zeros(nx)}
xcen_down = {'g1_mos1': np.zeros(nx),'g1_mos2': np.zeros(nx),'g1_pn': np.zeros(nx),\
        'g2_mos1': np.zeros(nx),'g2_mos2': np.zeros(nx),'g2_pn': np.zeros(nx),\
        'g4_mos1': np.zeros(nx),'g4_mos2': np.zeros(nx),'g4_pn': np.zeros(nx)}
xi2 = {'g1_mos1': np.zeros(nx),'g1_mos2': np.zeros(nx),'g1_pn': np.zeros(nx),\
        'g2_mos1': np.zeros(nx),'g2_mos2': np.zeros(nx),'g2_pn': np.zeros(nx),\
        'g4_mos1': np.zeros(nx),'g4_mos2': np.zeros(nx),'g4_pn': np.zeros(nx)}
xi2_up = {'g1_mos1': np.zeros(nx),'g1_mos2': np.zeros(nx),'g1_pn': np.zeros(nx),\
        'g2_mos1': np.zeros(nx),'g2_mos2': np.zeros(nx),'g2_pn': np.zeros(nx),\
        'g4_mos1': np.zeros(nx),'g4_mos2': np.zeros(nx),'g4_pn': np.zeros(nx)}
xi2_down = {'g1_mos1': np.zeros(nx),'g1_mos2': np.zeros(nx),'g1_pn': np.zeros(nx),\
        'g2_mos1': np.zeros(nx),'g2_mos2': np.zeros(nx),'g2_pn': np.zeros(nx),\
        'g4_mos1': np.zeros(nx),'g4_mos2': np.zeros(nx),'g4_pn': np.zeros(nx)}
#
for j in ["g1","g2","g4"]:
    labx = tt[j]['label']
    for ik,k in enumerate(xx):
        for l in ["mos1","mos2","pn"]:
            check = 'fit{}_1gauss_{}ks_{}.log'.format(j[::-1],k,l)
            ix = np.where(labx == check)[0]
            if (len(ix) != 10):
                raise Exception
            ttx = tt[j][ix]
            xcen["%s_%s"%(j,l)][ik] = np.mean(ttx['cen'])
            xi2["%s_%s"%(j,l)][ik] = np.mean(ttx['chi2'])
            xcen_down["%s_%s"%(j,l)][ik] = xcen["%s_%s"%(j,l)][ik] - np.min(ttx['cen'])
            xcen_up["%s_%s"%(j,l)][ik] = np.max(ttx['cen']) - xcen["%s_%s"%(j,l)][ik]
            xi2_down["%s_%s"%(j,l)][ik] = xi2["%s_%s"%(j,l)][ik] - np.min(ttx['chi2'])
            xi2_up["%s_%s"%(j,l)][ik] = np.max(ttx['chi2']) - xi2["%s_%s"%(j,l)][ik]

#%%
#
# now plotting
#
fig = plt.figure(figsize=(10,8))
ax = fig.subplots()

for j in ["g1","g2","g4"]:
    ax.errorbar(xx,xcen["%s_mos1"%j],yerr=(xcen_down["%s_mos1"%j],xcen_up["%s_mos1"%j]),fmt='o',color='blue')
    ax.errorbar(xx+2,xcen["%s_mos2"%j],yerr=(xcen_down["%s_mos2"%j],xcen_up["%s_mos2"%j]),fmt='o',color='red')
    ax.errorbar(xx+4,xcen["%s_pn"%j],yerr=(xcen_down["%s_pn"%j],xcen_up["%s_pn"%j]),fmt='o',color='green')
ax.hlines(6.41,10,120,linestyle='dashed',label='Input centre')
ax.set_ylim((6.3,6.5))
ax.set_xlabel("Exposure time (ks)")
ax.set_ylabel("Fe K best fit centre (keV)")
ax.set_title("Fe K simulations with 1, 2 or 4 Gaussian models, Rev: {}".format(rev))
ax.grid()
#
plt.legend(["input","MOS1","MOS2","PN"],loc=4)
#plt.savefig('fig_sims10_rev{}.png'.format(rev),dpi=100)
plt.show()
plt.close()

#%%
# chi-2 plot
fig = plt.figure(figsize=(10,8))
ax = fig.subplots()

fmt = {'g1': 'x', 'g2': 's', 'g4': 'o'}
for j in ["g1","g2","g4"]:
    ax.errorbar(xx,xi2["%s_mos1"%j],yerr=(xi2_down["%s_mos1"%j],xi2_up["%s_mos1"%j]),fmt=fmt[j],color='blue',label='MOS1 %s'%j)
    ax.errorbar(xx+2,xi2["%s_mos2"%j],yerr=(xi2_down["%s_mos2"%j],xi2_up["%s_mos2"%j]),fmt=fmt[j],color='red',label='MOS2 %s'%j)
    ax.errorbar(xx+4,xi2["%s_pn"%j],yerr=(xi2_down["%s_pn"%j],xi2_up["%s_pn"%j]),fmt=fmt[j],color='green',label='PN %s'%j)
ax.set_xlabel("Exposure time (ks)")
ax.set_ylabel("Ch-sqr")
ax.set_title("Fe K simulations with 1, 2 or 4 Gaussian models, Rev: {}".format(rev))
ax.grid()
#
plt.legend(loc=2)
plt.savefig('fig_sims10_rev{}_chi2.png'.format(rev),dpi=100)
#plt.show()
plt.close()



