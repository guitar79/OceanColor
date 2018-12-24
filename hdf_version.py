#!/usr/bin/env python

import os
import subprocess

#ifile = '/rsclass/data/netCDF4_examples/A2003060134000.L2_LAC_AT108.nc' 
ifile = '/rsclass/data/netCDF4_examples/A2014240170500.L2_LAC_OC.x.hdf' 
 

p= subprocess.Popen(["file", ifile], stdout=subprocess.PIPE)
output, err = p.communicate()
print (output)


version_4_check= "version 4" in output
version_5_check= "version 5" in output


if version_4_check == 1:
    ftype= 'hdf4'
if version_5_check == 1:
    ftype= 'hdf5'

print (ftype)

#return, ftype
    

