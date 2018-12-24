#!/usr/bin/env python
from my_general_read_utilities import *
import matplotlib.pyplot as plt
import numpy as np


sstFile = '/Users/md875/data/sst/ghrsst/20100609-JPL_OUROCEAN-L4UHfnd-GLOB-v01-fv01_0-G1SST.nc'
sst = read_ghrsst_pf_sst(sstFile)
print '\nSST File shape: ', sst.shape

newshape = ( sst.shape[0] / 100 , sst.shape[1] / 100 )

slices = [ slice(None, None, old/new) for old, new in zip(sst.shape, newshape) ]
new_sst = sst[slices]

print '\nnew SST shape', new_sst.shape

#image output
mycmap = plt.get_cmap('spectral')
mycmap.set_bad('k')
plt.imshow( np.flipud(new_sst), cmap = mycmap, vmin = -2, vmax=35 )
plt.show()
