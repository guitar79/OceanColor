#!usr/bin/env python
from my_hdf_cdf_utilities import *
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors
from mpl_toolkits.basemap import Basemap

fname = '/Users/md875/data/tutorial_data/S1998148172338_chlor_a.f999x999'
f = open(fname)
mapped_data = np.fromfile(f, dtype = np.float32)
f.close()
mapped_data = mapped_data.reshape([999,999])

#set lan lon ( a priori)

north = 46.0; west = -72.0 ;
south = 37.0; east = -63.0 ;

mapPara =  Basemap(projection = 'cyl',\
llcrnrlon = west, llcrnrlat = south, \
urcrnrlon = east, urcrnrlat = north, \
resolution = 'h')
					
# llcrnr = lower left corner , resolution = 'h' --> high

#set color palette
colorScale = plt.get_cmap('spectral')
colorScale.set_bad('k')

# flip array upside down - for mapping only
mapped_data = np.flipud(mapped_data)

# display the mapped img in map projection window
mapPara.imshow(mapped_data, cmap = colorScale, \
				vmin = 0.01, vmax = 30.0, \
				norm =  matplotlib.colors.LogNorm())

# draw coastline, lat lon grid & colorbar
mapPara.drawcoastlines()
mapPara.fillcontinents(color = 'grey', lake_color = 'white')

#labels (left, right, top,bottom)
parallel = np.arange(south, north, 2.)
mapPara.drawparallels(parallel, labels = [True, False, False, False])

meridian = np.arange( west, east, 2.)
mapPara.drawmeridians(meridian, labels = [False, False, True, False])

#save the map imgs a png
plt.savefig('/Users/md875/Desktop/test.png', bbox_inches = 'tight')

# show the plot to monitor
plt.show()

plt.close()
