#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors

#read the SST file 
fname1 = '/Users/md875/data/tutorial_data/sst/2004.mon06.f774x843'
f1 = open(fname1)
data1 = np.fromfile( fname1, dtype = np.float32 )
f1.close()

# row & col no (change this)
rowNo = 843.0 ; colNo = 774.0;

# reshp to 2D matrix
data1Array = data1.reshape([rowNo,colNo])

##subsetting
#bookkeeping & converting vars

# limit/extent  (draw grid)   (change this)
northLimit = 62.0;				 westLimit = -82.0;  eastLimit = -48.0; 
southLimit = 25.0; 

# automatically calculated
del_lat = abs(northLimit-southLimit);    	del_lon = abs(westLimit-eastLimit);

# lat/lon for subset img (change this)
northLimit_subset = 46.0;				 westLimit_subset = -72.0;  eastLimit_subset = -63.0; 
southLimit_subset = 37.0; 

# conversion automatically calculated
row_subset_North = (northLimit - northLimit_subset) * (rowNo-1) / del_lat ;   	col_subset_West = (westLimit_subset - westLimit) * (colNo-1) / del_lon ;  
row_subset_South = (northLimit - southLimit_subset) * (rowNo-1) / del_lat ;   	col_subset_East = (eastLimit_subset - westLimit) * (colNo-1) / del_lon ;

# subsetting operation
subset_img = data1Array[ row_subset_North:row_subset_South , col_subset_West:col_subset_East ]

##display imgs

#color scale
colorScale = plt.get_cmap('spectral')
colorScale.set_bad('k')

#create & display figure
plt.figure(1)
plt.subplot(221)
plt.imshow(data1Array, cmap = colorScale, vmin = 0.0, vmax = 30.0)
plt.colorbar()

plt.subplot(222)
plt.imshow(subset_img, cmap = colorScale, vmin = 0.0, vmax = 30.0)
plt.show()

