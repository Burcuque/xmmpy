#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  1 12:11:13 2019

Check the $TMPDIR when processing on the grid

@author: ivaltchanov
"""
import os
import shutil

from fnmatch import fnmatch, filter
#
def include_patterns(*patterns):
    """Factory function that can be used with copytree() ignore parameter.

    Arguments define a sequence of glob-style patterns
    that are used to specify what files to NOT ignore.
    Creates and returns a function that determines this for each directory
    in the file hierarchy rooted at the source directory when used with
    shutil.copytree().
    """
    def _ignore_patterns(path, names):
        keep = set(name for pattern in patterns
                            for name in filter(names, pattern))
        ignore = set(name for name in names
                        if name not in keep and not os.path.isdir(os.path.join(path, name)))
        return ignore
    return _ignore_patterns
#
home = os.path.expanduser('~')
#
wdir = home + '/IVAN/Cu-line'
os.chdir(wdir)

obsid = '0745110101'
srcdir = os.path.join(wdir,obsid)

#
# will copy this folder to the grid node, only the ODF files
#
fout = open(f'{wdir}/check_content_node.txt','w')

if ('TMPDIR' in os.environ.keys()):
    tmpdir = os.environ["TMPDIR"]
    destdir = os.path.join(tmpdir,obsid)
    print (f"Destination is {destdir}")
    #
    os.chdir(wdir)
    shutil.copytree(obsid,destdir,ignore=shutil.ignore_patterns('pn*','nocti*', 'ccf*','cti*', 'PN_LW*'))
    #shutil.copytree(obsid,destdir,ignore=include_patterns('*.FIT', '*.ASC','*.SUM'))
    c = os.listdir(path = tmpdir)
    print (c,file=fout)
    c = os.listdir(path = destdir)
    print (c,file=fout)
    check = os.path.join(destdir,"nocti")
    if (os.path.isdir(check)):
        c = os.listdir(path = check)
        print (c,file=fout)
    fout.close()
    # creating a folder with a file
    #
#    ppsdir = f"{tmpdir}/pps"
#    nodename = os.path.basename(tmpdir)
#    if (not os.path.isdir(ppsdir)):
#        print (f"Creating a folder {ppsdir}")
#        os.mkdir(ppsdir)
#    fname = f"{ppsdir}/{nodename}_test.dat"
#    fout = open(fname,"w")
#    print ("Saving to a file")
#    print (os.environ,file=fout)
#    fout.close()
#    #
#    # now copy back the file to home
#    #
#    copy_folder = f"{home}/tmp/XMM_data"
#    shutil.copy(fname,copy_folder)
#    print (f"Copy of {fname} to {copy_folder} finished")    
else:
    print ("$TMPDIR not in environment")
print ("End")
#
