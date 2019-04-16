#
# master script to run XPSEC simulations
#
import os
import subprocess
import sys
import argparse
#
parser = argparse.ArgumentParser(description='XSPEC simulations for Fe K')
parser.add_argument('instrument', type=str,
                    help='The instrument to simulate, can be "mos1", "mos2" or "pn"')
parser.add_argument('expo', type=int,
                    help='The exposure time in seconds')
parser.add_argument('--prefix', type=str, default="test",
                    help='Prefix of the output spectra')
args = parser.parse_args()

home = os.path.expanduser('~')
simdir = "{}/Dropbox/Work/XMM/simulations".format(home)

inst = args.instrument
texp = args.expo
pref = args.prefix
#
pha_file = "{}/{}_4gauss_{}ks_{}.pha".format(simdir,pref,int(texp/1000.0),inst)
xscript = "{}/{}_4gauss_{}.xcm".format(simdir,pref,inst)
data_file = "{}/{}_4gauss_{}.dat".format(simdir,pref,inst)

# delete the input parameters file if it exists
if (os.path.isfile(data_file)):
    os.remove(data_file)
#
with open(xscript,"w") as fout:
    print ("response 1 %s.rmf"%inst,file=fout)
    print ("arf 1 %s.arf"%inst,file=fout)
    # now the total model
    print ("abund wilm",file=fout)
    print ("model phabs*(powerlaw + bbody + gauss + gauss + gauss + gauss) & /*",file=fout)
    # and the parameters
    # phabs 
    print ("newpar 1 9.9e-2",file=fout)
    # powerlaw
    print ("newpar 2 2.4",file=fout)
    print ("newpar 3 1.8e-2",file=fout)
    # the blackbody
    print ("newpar 4 2.0",file=fout)
    print ("newpar 5 2.4e-4",file=fout)
    # line 1, narrow, the strongest
    print ("newpar 6 6.41",file=fout)
    print ("newpar 7 5.0e-2",file=fout)
    print ("newpar 8 3.0e-4",file=fout)
    # line 2, narrow
    print ("newpar 9 6.97",file=fout)
    print ("newpar 10 5.0e-2",file=fout)
    print ("newpar 11 1.0e-4",file=fout)
    # line 3, broad
    print ("newpar 12 6.6",file=fout)
    print ("newpar 13 0.3",file=fout)
    print ("newpar 14 5.0e-4",file=fout)
    # line 4, broad
    print ("newpar 15 6.2",file=fout)
    print ("newpar 16 0.3",file=fout)
    print ("newpar 17 3.0e-4",file=fout)
    #
    print ("fakeit none",file=fout)
    print ("%s.rmf"%inst,file=fout)
    print ("%s.arf"%inst,file=fout)
    print ("y",file=fout)
    print ("\\n",file=fout)
    print (os.path.basename(pha_file),file=fout)
    print ("%i, 1.0, 1.0, 0.0"%texp,file=fout)
    print ("save model {}".format(os.path.basename(data_file)),file=fout)
    print ("quit",file=fout)
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
xscript2 = "fit_1gauss_{}.xcm".format(inst)
log_file = "fit_1gauss_{}ks_{}.log".format(int(texp/1000),inst)
plt_file = "fit_1gauss_{}ks_{}.ps".format(int(texp/1000),inst)
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
    print ("newpar 2 2.4",file=fout)
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
    print ("error stop 20,,4-5".format(data_file),file=fout)
    print ("cpd {}/cps".format(plt_file),file=fout)
    print ("plot data ratio",file=fout)
    print ("query yes",file=fout)    
    print ("exit",file=fout)
#
print ("*** Fit script {} done".format(os.path.basename(xscript2)))
#
# now execute XSPEC
#
xspec_comm = "xspec - {}".format(xscript2)
print ("Running XSPEC to fit the simulated spectrum {}".format(pha_file))
try:
    result = subprocess.run(xspec_comm, shell=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
    retcode = result.returncode
    if retcode != 0:
        print("XSPEC fit was terminated by signal", -retcode, file=sys.stderr)
        raise Exception
    else:
        print("XSPEC fit returned", retcode, file=sys.stderr)
except OSError as e:
    print("Execution of XSPEC fit failed:", e, file=sys.stderr)
#
# Now parse the output and save the results in a csv file
#
with open(log_file,'r') as fin:
    xline = fin.readline()
    while xline:
        if ('LineE' in xline):
            cen = xline.split()[6] # the line best fit centroid
        if ('Reduced' in xline):
            xi2r = xline.split()[4] # the reduced chi-square
            dof = xline.split()[6] # the reduced chi-square
        if ('Confidence' in xline):
            xline = fin.readline()
            cen_lo = xline.split()[2] # the reduced chi-square
            cen_hi = xline.split()[3] # the reduced chi-square
        xline = fin.readline()
#
out_file = "simulation_results.csv"
with open(out_file,'a') as fout:
    print ("{},{},{},{},{},{},{},{}".format(log_file,inst,int(texp/1000.0),cen,cen_lo,cen_hi,xi2r,dof),file=fout)
#

