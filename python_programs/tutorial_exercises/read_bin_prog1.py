#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors


#open files & save data into arrays
dataFile = open('/Users/md875/data/tutorial_data/S1998148172338_chlor_a.f999x999')
data1 = np.fromfile(dataFile, dtype = np.float32)
dataFile.close()

#reshape 2D data to the 2D matrixs
data1 = data1.reshape([999,999])

print 'shape of the data1 array is \n\n', data1.shape

# define color scale
myChlMap = plt.get_cmap('spectral')  # loads color rainbow palette
myChlMap.set_bad('k')  #set NaN values to display as black

#creates & display the figure to your monitor with data in logscale
plt.figure(1)
plt.imshow(data1, myChlMap, vmin = 0.01, vmax = 20.0, norm = matplotlib.colors.LogNorm())
plt.colorbar()
plt.show()