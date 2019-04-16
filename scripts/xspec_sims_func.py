#
# function script to run XPSEC simulations
#
import os
import subprocess
import sys
#
#import numpy as np
#
#%%
def make_sim_4gauss(inst,texp,pref="test", simdir=".", plot_model=False):
    #
    # Simulate 4 gaussians at 6-7 keV, two narrow and two broad lines
    #
    os.chdir(simdir)
    pha_file = "{}_4gauss_{}ks_{}.pha".format(pref,int(texp/1000.0),inst)
    xscript = "{}_4gauss_{}.xcm".format(pref,inst)
    data_file = "{}_4gauss_{}.dat".format(pref,inst)
    # delete the input parameters file if it exists
    if (os.path.isfile(data_file)):
        os.remove(data_file)
    #
    with open(xscript,"w") as fout:
        print ("response 1 %s.rmf"%inst,file=fout)
        print ("arf 1 %s.arf"%inst,file=fout)
        # now the total model
        print ("abund wilm",file=fout)
        print ("model phabs*(powerlaw + gauss + gauss + gauss + gauss) & /*",file=fout)
        # and the parameters
        # phabs 
        print ("newpar 1 9.9e-2",file=fout)
        # powerlaw
        print ("newpar 2 1.7",file=fout)
        print ("newpar 3 1.8e-2",file=fout)
        # line 1, narrow, the strongest
        print ("newpar 4 6.41",file=fout)
        print ("newpar 5 5.0e-2",file=fout)
        print ("newpar 6 3.0e-4",file=fout)
        # line 2, narrow
        print ("newpar 7 6.97",file=fout)
        print ("newpar 8 5.0e-2",file=fout)
        print ("newpar 9 1.0e-4",file=fout)
        # line 3, broad
        print ("newpar 10 6.6",file=fout)
        print ("newpar 11 0.3",file=fout)
        print ("newpar 12 5.0e-4",file=fout)
        # line 4, broad
        print ("newpar 13 6.2",file=fout)
        print ("newpar 14 0.3",file=fout)
        print ("newpar 15 3.0e-4",file=fout)
        #
        print ("fakeit none",file=fout)
        print ("%s.rmf"%inst,file=fout)
        print ("%s.arf"%inst,file=fout)
        print ("y",file=fout)
        print ("\\n",file=fout)
        print (pha_file,file=fout)
        print ("%i, 1.0, 1.0, 0.0"%texp,file=fout)
        print ("save model {}".format(data_file),file=fout)
        if (plot_model):
            mplt_file = "model4g_zoom.ps"
            print ('setplot energy',file=fout)
            print ('setplot comm lab top "4g model"',file=fout)
            print ('setplot comm rescale y 1.0e-4 1.0e-2',file=fout)
            print ('setplot comm rescale x 4.0 10.0',file=fout)
            print ("cpd {}/cps".format(mplt_file),file=fout)
            print ("plot model",file=fout)        
            mplt_file = "data4g_zoom.ps"
            print ('setplot energy',file=fout)
            print ('setplot comm lab top "4g data {}ks, {}"'.format(int(texp/1000),inst),file=fout)
            print ('setplot comm rescale y 0.0 0.5',file=fout)
            print ('setplot comm rescale x 4.0 10.0',file=fout)
            print ("cpd {}/cps".format(mplt_file),file=fout)
            print ("plot data",file=fout)        
        print ("exit",file=fout)
    #
    print ("*** script {} ready".format(xscript))
    #
    # now execute XSPEC
    #
    xspec_comm = "xspec - {}".format(xscript)
    print ("Running XSPEC to simulate a spectrum for {}".format(inst))
    #
    try:
        result = subprocess.run(xspec_comm, shell=True)
        #result = subprocess.run(xspec_comm, shell=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
        retcode = result.returncode
        if retcode != 0:
            print("XSPEC simulation was terminated by signal", -retcode, file=sys.stderr)
            raise Exception
        else:
            print("XSPEC simulation returned", retcode, file=sys.stderr)
    except OSError as e:
        print("Execution of XSPEC simulation failed:", e, file=sys.stderr)
    #
    # NOW fit the simulated PHA file, using a simple model of only phabs*(powerlaw + gauss)
    #
    xscript2 = "fit4g_1gauss_{}.xcm".format(inst)
    log_file = "fit4g_1gauss_{}ks_{}.log".format(int(texp/1000),inst)
    plt_file = "fit4g_1gauss_{}ks_{}.ps".format(int(texp/1000),inst)
    #
    with open(xscript2,"w") as fout:
        print ("data 1:1 {}".format(pha_file),file=fout)
        print ("ignore bad",file=fout)
        # will only fit in energy range [4,10] keV
        print ("ignore 1:**-4.0,10.0-**",file=fout)
        print ("response 1 %s.rmf"%inst,file=fout)
        print ("arf 1 %s.arf"%inst,file=fout)
        print ("abund wilm",file=fout)
        print ('setplot energy',file=fout)
        # now the total model
        print ("model phabs*(powerlaw + gauss) & /*",file=fout)
        # and the parameters
        # phabs 
        print ("newpar 1 9.9e-2",file=fout)
        # powerlaw
        print ("newpar 2 1.7",file=fout)
        print ("newpar 3 1.8e-2",file=fout)
        # line 1, narrow, the strongest
        print ("newpar 4 6.41",file=fout)
        print ("newpar 5 5.0e-2",file=fout)
        print ("newpar 6 3.0e-4",file=fout)
        print ("query yes",file=fout)    
        print ("fit 100".format(data_file),file=fout)    
        print ("query yes".format(data_file),file=fout)
        print ("log {}".format(log_file),file=fout)
        print ("show all",file=fout)
#        print ("error stop 20,,max 3.0, 4-5".format(data_file),file=fout)
        print ('setplot comm lab top "1g fit to 4g model: {} ks, {}"'.format(int(texp/1000),inst),file=fout)
        print ("setplot rebin 50 15 1",file=fout)
        print ("cpd {}/cps".format(plt_file),file=fout)
        print ("plot data ratio",file=fout)
        print ("query yes",file=fout)    
        print ("exit",file=fout)
    #
    print ("*** Fit script {} done".format(xscript2))
    #
    # now execute XSPEC
    #
    xspec_comm = "xspec - {}".format(xscript2)
    print ("Running XSPEC to fit the simulated spectrum {}".format(pha_file))
    try:
        result = subprocess.run(xspec_comm, shell=True)
        #result = subprocess.run(xspec_comm, shell=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
        retcode = result.returncode
        if retcode != 0:
            print("XSPEC fit was terminated by signal", -retcode, file=sys.stderr)
            raise Exception
        else:
            print("XSPEC fit returned", retcode, file=sys.stderr)
    except OSError as e:
        print("Execution of XSPEC fit failed:", e, file=sys.stderr)
    return log_file
#
#%%
def make_sim_2gauss(inst,texp,pref="test", simdir=".", plot_model=False):
    #
    # Simulate 2 gaussians at 6.41 narrow and 6.2 broad
    #
    os.chdir(simdir)
    pha_file = "{}_2gauss_{}ks_{}.pha".format(pref,int(texp/1000.0),inst)
    xscript = "{}_2gauss_{}.xcm".format(pref,inst)
    data_file = "{}_2gauss_{}.dat".format(pref,inst)
    # delete the input parameters file if it exists
    if (os.path.isfile(data_file)):
        os.remove(data_file)
    #
    with open(xscript,"w") as fout:
        print ("response 1 %s.rmf"%inst,file=fout)
        print ("arf 1 %s.arf"%inst,file=fout)
        # now the total model
        print ("abund wilm",file=fout)
        print ("model phabs*(powerlaw + gauss + gauss) & /*",file=fout)
        # and the parameters
        # phabs 
        print ("newpar 1 9.9e-2",file=fout)
        # powerlaw
        print ("newpar 2 1.7",file=fout)
        print ("newpar 3 1.8e-2",file=fout)
        # line 1, narrow, the strongest
        print ("newpar 4 6.41",file=fout)
        print ("newpar 5 5.0e-2",file=fout)
        print ("newpar 6 3.0e-4",file=fout)
        # line 2, broad
        print ("newpar 7 6.2",file=fout)
        print ("newpar 8 0.3",file=fout)
        print ("newpar 9 3.0e-4",file=fout)
        #
        print ("fakeit none",file=fout)
        print ("%s.rmf"%inst,file=fout)
        print ("%s.arf"%inst,file=fout)
        print ("y",file=fout)
        print ("\\n",file=fout)
        print (os.path.basename(pha_file),file=fout)
        print ("%i, 1.0, 1.0, 0.0"%texp,file=fout)
        print ("save model {}".format(data_file),file=fout)
        if (plot_model):
            mplt_file = "model2g_zoom.ps"
            print ('setplot energy',file=fout)
            print ('setplot comm lab top "2g model"',file=fout)
            print ('setplot comm rescale y 1.0e-4 1.0e-1',file=fout)
            print ('setplot comm rescale x 4.0 10.0',file=fout)
            print ("cpd {}/cps".format(mplt_file),file=fout)
            print ("plot model",file=fout)        
            mplt_file = "data2g_zoom.ps"
            print ('setplot energy',file=fout)
            print ('setplot comm lab top "2g data {}ks, {}"'.format(int(texp/1000),inst),file=fout)
            print ('setplot comm rescale y 0.0 1.0',file=fout)
            print ('setplot comm rescale x 4.0 10.0',file=fout)
            print ("cpd {}/cps".format(mplt_file),file=fout)
            print ("plot data",file=fout)        
        print ("exit",file=fout)
    #
    print ("*** script {} ready".format(xscript))
    #
    # now execute XSPEC
    #
    xspec_comm = "xspec - {}".format(xscript)
    print ("Running XSPEC to simulate a spectrum for {}".format(inst))
    #
    try:
        result = subprocess.run(xspec_comm, shell=True)
        #result = subprocess.run(xspec_comm, shell=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
        retcode = result.returncode
        if retcode != 0:
            print("XSPEC simulation was terminated by signal", -retcode, file=sys.stderr)
            raise Exception
        else:
            print("XSPEC simulation returned", retcode, file=sys.stderr)
    except OSError as e:
        print("Execution of XSPEC simulation failed:", e, file=sys.stderr)
    #
    # NOW fit the simulated PHA file, using a simple model of only phabs*(powerlaw + gauss)
    #
    xscript2 = "fit2g_1gauss_{}.xcm".format(inst)
    log_file = "fit2g_1gauss_{}ks_{}.log".format(int(texp/1000),inst)
    plt_file = "fit2g_1gauss_{}ks_{}.ps".format(int(texp/1000),inst)
    #
    with open(xscript2,"w") as fout:
        print ("data 1:1 {}".format(pha_file),file=fout)
        print ("ignore bad",file=fout)
        # will only fit in energy range [4,10] keV
        print ("ignore 1:**-4.0,10.0-**",file=fout)
        print ("response 1 %s.rmf"%inst,file=fout)
        print ("arf 1 %s.arf"%inst,file=fout)
        print ("abund wilm",file=fout)
        # now the total model
        print ("model phabs*(powerlaw + gauss) & /*",file=fout)
        # and the parameters
        # phabs 
        print ("newpar 1 9.9e-2",file=fout)
        # powerlaw
        print ("newpar 2 1.7",file=fout)
        print ("newpar 3 1.8e-2",file=fout)
        # line 1, narrow, the strongest
        print ("newpar 4 6.41",file=fout)
        print ("newpar 5 5.0e-2",file=fout)
        print ("newpar 6 3.0e-4",file=fout)
        print ("query yes",file=fout)    
        print ("fit 100".format(data_file),file=fout)    
        print ("query yes".format(data_file),file=fout)
        print ("log {}".format(log_file),file=fout)
        print ("show all",file=fout)
#        print ("error stop 20,,max 3.0, 4-5".format(data_file),file=fout)
        print ("setplot energy",file=fout)
        print ('setplot comm lab top "1g fit to 2g model: {} ks, {}"'.format(int(texp/1000),inst),file=fout)
        print ("setplot rebin 50 15 1",file=fout)
        print ("cpd {}/cps".format(plt_file),file=fout)
        print ("plot data ratio",file=fout)
        print ("query yes",file=fout)    
        print ("exit",file=fout)
    #
    print ("*** Fit script {} done".format(xscript2))
    #
    # now execute XSPEC
    #
    xspec_comm = "xspec - {}".format(xscript2)
    print ("Running XSPEC to fit the simulated spectrum {}".format(pha_file))
    try:
        result = subprocess.run(xspec_comm, shell=True)
        #result = subprocess.run(xspec_comm, shell=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
        retcode = result.returncode
        if retcode != 0:
            print("XSPEC fit was terminated by signal", -retcode, file=sys.stderr)
            raise Exception
        else:
            print("XSPEC fit returned", retcode, file=sys.stderr)
    except OSError as e:
        print("Execution of XSPEC fit failed:", e, file=sys.stderr)
    return log_file
#
#%%
def make_sim_1gauss(inst,texp,pref="test", simdir=".", plot_model=False):
    #
    # Simulate 1 gaussian at 6.41 keV
    #
    os.chdir(simdir)
    pha_file = "{}_1gauss_{}ks_{}.pha".format(pref,int(texp/1000.0),inst)
    xscript = "{}_1gauss_{}.xcm".format(pref,inst)
    data_file = "{}_1gauss_{}.dat".format(pref,inst)
    # delete the input parameters file if it exists
    if (os.path.isfile(data_file)):
        os.remove(data_file)
    #
    with open(xscript,"w") as fout:
        print ("response 1 %s.rmf"%inst,file=fout)
        print ("arf 1 %s.arf"%inst,file=fout)
        # now the total model
        print ("abund wilm",file=fout)
        print ("model phabs*(powerlaw + gauss) & /*",file=fout)
        # and the parameters
        # phabs 
        print ("newpar 1 9.9e-2",file=fout)
        # powerlaw
        print ("newpar 2 1.7",file=fout)
        print ("newpar 3 1.8e-2",file=fout)
        # line 1, narrow, the strongest
        print ("newpar 4 6.41",file=fout)
        print ("newpar 5 5.0e-2",file=fout)
        print ("newpar 6 3.0e-4",file=fout)
        #
        print ("fakeit none",file=fout)
        print ("%s.rmf"%inst,file=fout)
        print ("%s.arf"%inst,file=fout)
        print ("y",file=fout)
        print ("\\n",file=fout)
        print (os.path.basename(pha_file),file=fout)
        print ("%i, 1.0, 1.0, 0.0"%texp,file=fout)
        print ("save model {}".format(data_file),file=fout)
        if (plot_model):
            mplt_file = "model1g_zoom.ps"
            print ('setplot energy',file=fout)
            print ('setplot comm lab top "1g model"',file=fout)
            print ('setplot comm rescale y 1.0e-4 1.0e-2',file=fout)
            print ('setplot comm rescale x 4.0 10.0',file=fout)
            print ("cpd {}/cps".format(mplt_file),file=fout)
            print ("plot model",file=fout)        
            mplt_file = "data1g_zoom.ps"
            print ('setplot energy',file=fout)
            print ('setplot comm lab top "1g data {}ks, {}"'.format(int(texp/1000),inst),file=fout)
            print ('setplot comm rescale y 0.0 1.0',file=fout)
            print ('setplot comm rescale x 4.0 10.0',file=fout)
            print ("cpd {}/cps".format(mplt_file),file=fout)
            print ("plot data",file=fout)        
        print ("exit",file=fout)
    #
    print ("*** script {} ready".format(xscript))
    #
    # now execute XSPEC
    #
    xspec_comm = "xspec - {}".format(xscript)
    print ("Running XSPEC to simulate a spectrum for {}".format(inst))
    #
    try:
        result = subprocess.run(xspec_comm, shell=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
        retcode = result.returncode
        if retcode != 0:
            print("XSPEC simulation was terminated by signal", -retcode, file=sys.stderr)
            raise Exception
        else:
            print("XSPEC simulation returned", retcode, file=sys.stderr)
    except OSError as e:
        print("Execution of XSPEC simulation failed:", e, file=sys.stderr)
    #
    # NOW fit the simulated PHA file, using a simple model of only phabs*(powerlaw + gauss)
    #
    xscript2 = "fit1g_1gauss_{}.xcm".format(inst)
    log_file = "fit1g_1gauss_{}ks_{}.log".format(int(texp/1000),inst)
    plt_file = "fit1g_1gauss_{}ks_{}.ps".format(int(texp/1000),inst)
    #
    with open(xscript2,"w") as fout:
        print ("data 1:1 {}".format(pha_file),file=fout)
        print ("ignore bad",file=fout)
        # will only fit in energy range [4,10] keV
        print ("ignore 1:**-4.0,10.0-**",file=fout)
        print ("response 1 %s.rmf"%inst,file=fout)
        print ("arf 1 %s.arf"%inst,file=fout)
        print ("abund wilm",file=fout)
        # now the total model
        print ("model phabs*(powerlaw + gauss) & /*",file=fout)
        # and the parameters
        # phabs 
        print ("newpar 1 9.9e-2",file=fout)
        # powerlaw
        print ("newpar 2 1.7",file=fout)
        print ("newpar 3 1.8e-2",file=fout)
        # line 1, narrow, the strongest
        print ("newpar 4 6.41",file=fout)
        print ("newpar 5 5.0e-2",file=fout)
        print ("newpar 6 3.0e-4",file=fout)
        print ("query yes",file=fout)    
        print ("fit 100".format(data_file),file=fout)    
        print ("query yes".format(data_file),file=fout)
        print ("log {}".format(log_file),file=fout)
        print ("show all",file=fout)
#        print ("error stop 20,,max 3.0,4-5".format(data_file),file=fout)
        print ("setplot energy",file=fout)
        print ('setplot comm lab top "1g fit to 1g model: {} ks, {}"'.format(int(texp/1000),inst),file=fout)
        print ("setplot rebin 50 15 1",file=fout)
        print ("cpd {}/cps".format(plt_file),file=fout)
        print ("plot data ratio",file=fout)
        print ("query yes",file=fout)    
        print ("exit",file=fout)
    #
    print ("*** Fit script {} done".format(xscript2))
    #
    # now execute XSPEC
    #
    xspec_comm = "xspec - {}".format(xscript2)
    print ("Running XSPEC to fit the simulated spectrum {}".format(pha_file))
    try:
        result = subprocess.run(xspec_comm, shell=True)
        #result = subprocess.run(xspec_comm, shell=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
        retcode = result.returncode
        if retcode != 0:
            print("XSPEC fit was terminated by signal", -retcode, file=sys.stderr)
            raise Exception
        else:
            print("XSPEC fit returned", retcode, file=sys.stderr)
    except OSError as e:
        print("Execution of XSPEC fit failed:", e, file=sys.stderr)
    return log_file
#
#%%
def parse_logfile(logfile,outfile):
    #
    # Parse the output and save the results in a csv file
    #
    with open(logfile,'r') as fin:
        xline = fin.readline()
        while xline:
            split_line = xline.split()
            if ('LineE' in xline):
                cen = split_line[6] # the line best fit centroid
                cen_err = split_line[8] # the line best fit centroid
                print (cen)
            if ('Reduced' in xline):
                xi2r = split_line[4] # the reduced chi-square
                dof = split_line[6] # the reduced chi-square
#            if ('Confidence Range' in xline):
#                xline = fin.readline()
#                found = False
#                while (xline):
#                    if ('#     4' in xline):
#                        split_line = xline.split()
#                        cen_lo = split_line[2] # lower 1-sigma
#                        cen_hi = split_line[3] # high 1-sigma
#                        found = True
#                        break
#                    else:
#                        xline = fin.readline()
#                if (not found): 
#                    cen_lo = cen - cen_err # in case of problem
#                    cen_hi = cen + cen_err # 
#            #
#            else:
#                cen_lo = cen - cen_err # in case of problem
#                cen_hi = cen + cen_err #
#            #                
            xline = fin.readline()
    #
    with open(outfile,'a') as fout:
       print ("{},{},{},{},{},{}".format(logfile,cen,cen_err,cen_err,xi2r,dof),file=fout)
       #print ("{},{},{},{},{},{}".format(logfile,cen,cen_lo,cen_hi,xi2r,dof),file=fout)
    #
    return

