#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun  1 16:18:01 2018

@author: ivaltchanov
"""

import os
import numpy as np

home = os.path.expanduser('~')
func_file = home + '/Dropbox/Work/XMM/scripts/xspec_sims_func.py'

simdir = home + "/XMM/simulations/rev0679"

exec(open(func_file).read())
#%%
#make_sim_4gauss("mos1",50000,pref="test4g",simdir=simdir,\
#                out_file="simulation_results_4gauss.csv", plot_model=True)
#make_sim_2gauss("mos2",50000,pref="test2g",simdir=simdir,\
#                out_file="simulation_results_2gauss.csv", plot_model=True)
#make_sim_1gauss("pn",50000,pref="test1g",simdir=simdir,\
#                out_file="simulation_results_1gauss.csv", plot_model=True)

#%%
expo_time = [30000,50000,80000,100000] #
nruns = 10

for iexp in expo_time:
    print ("===> Simulation for {} ks".format(int(iexp/1000)))
    for irun in np.arange(nruns):
        print ("    ===> Simulation \#%i"%(irun+1),end="")
        for j in ["mos1","mos2","pn"]:
            print (" {} 4g, ".format(j),end="")
            log_out = make_sim_4gauss(j,iexp,pref="test4g", simdir=simdir)
            parse_logfile(log_out,"simulation10_results_4gauss.csv")
            print (" {} 2g, ".format(j),end="")
            log_out = make_sim_2gauss(j,iexp,pref="test2g", simdir=simdir)
            parse_logfile(log_out,"simulation10_results_2gauss.csv")
            print (" {} 1g, ".format(j),end="")
            log_out = make_sim_1gauss(j,iexp,pref="test1g", simdir=simdir)
            parse_logfile(log_out,"simulation10_results_1gauss.csv")
        print ()
    print ()
#
print ("All done")
