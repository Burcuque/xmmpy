#
# parse the logfile
#

log_file = "fit_1gauss_50ks_pn.log"

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
print ("{}, {}, {}, {}, {}, {}".format(log_file,cen,cen_lo,cen_hi,xi2r,dof))
#
