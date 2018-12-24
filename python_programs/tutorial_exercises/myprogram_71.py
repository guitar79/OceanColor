#!/usr/bin/loc python
from my_hdf_cdf_utilities import *
import numpy as np
import matplotlib.pyplot as plt
sst_file = '/Users/md875/data/dataD/AQUA_SMI/A2002185.L3m_DAY_NSST_sst_4km.nc'
#read in unscaled SST data
dn_sst_array = read_hdf_prod( sst_file, 'sst')
print dn_sst_array

bad_loc = np.where(dn_sst_array == 0)

# use scale info from hdf_prod_scale to convert to SST
# converting DN value to geophysical vals (oC)
scale = 0.000717185; intercept = -2.0; 
sst = dn_sst_array * scale + intercept

# setting bad to nan
sst[ bad_loc[0], bad_loc[1] ] = np.nan

# read in quality file associated with SST data
qual = read_hdf_prod( sst_file, 'qual_sst' )

# NaN out the bad quality data
# missing data set = -1 , poor quality > 2
bad_qual = np.where( (qual == -1 ) | ( qual >2 ) )
sst[ bad_qual[0] , bad_qual[1] ] = np.nan

# Image the SST data
mycmap = plt.get_cmap('spectral')
mycmap.set_bad('k')
plt.imshow(sst, cmap = mycmap, vmin = -2, vmax = 35)
plt.show()