#
# Make a Latex table for the AGNs
#
import os
import numpy as np
from astropy.table import Table
#

redshift = {'ngc4151': 0.003262, 'ngc3227': 0.00386, 'mrk1048': 0.0427, 'ngc3783': 0.009755,\
            'ngc4593': 0.008344, 'ngc5506': 0.00589, 'mcg-5-23-16': 0.008226, 'ngc3516': 0.008816,\
            'ngc5548': 0.01627, 'ngc2992':  0.007296, 'ngc1566': 0.005036}

wdir = "/xdata/xcaldata/XMM/IVAN/PN_SW"
os.chdir(wdir)

fout = open("agn_table.tex",'w')
for x in redshift.keys():
    print (f"Doing {x}")
    tx = Table.read(f'{x}/output_xspec.csv')
    nx = len(tx)
    obs_list = tx['obsid']
    revs = tx['revol']
    print ("{} & {} & ".format(x,redshift[x]),end="",file=fout)
    for j in np.arange(nx):
        if (('pn' in tx['instrument'][j]) and ('Small' in tx['submode'][j])):
            print("{}\\_{}, ".format(revs[j],obs_list[j]),end="",file=fout)
    print("\\\\",file=fout)
#
fout.close()

    

