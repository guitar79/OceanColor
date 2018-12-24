#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors
import glob

# file & data
fileNameList = glob.glob('/Users/md875/data/tutorial_data/8day_files/*S*')

rowDim = 999; colDim = 999;

#initializing
cumulativeDataArray = np.zeros( [1, rowDim * colDim] )
cumulativeCount = np.zeros( [1,rowDim * colDim] )

# looping for cumulating/summation for each file
for eachFile in fileNameList:
	f1 = open(eachFile)
	data1 = np.fromfile( f1, dtype = np.float32 )
	cumulativeDataArray += np.nan_to_num(data1)     #returns 0 for NAN and high numerical value for infinite and ultimately returns the original no if neither
	cumulativeCount += (~(np.isnan(data1))).astype(int)

avgDataArray = cumulativeDataArray/cumulativeCount
avgDataMatrix = avgDataArray.reshape([ rowDim, colDim ])

## plotting color map
#color scale
colorScale = plt.get_cmap('spectral')
colorScale.set_bad('k')

#create & display figure
plt.figure(1)
plt.imshow(avgDataMatrix, cmap = colorScale, vmin = 0.01, vmax = 20.0, norm = matplotlib.colors.LogNorm() )
plt.colorbar()
plt.show()

