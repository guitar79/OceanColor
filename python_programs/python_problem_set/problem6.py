#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors
import os
import glob

fileList = glob.glob('/Users/md875/data/tutorial_data/dly_swf1998/*1998*')

#dataSize
rowDim = 999; colDim = 999;

#initializing
cumulativeDataArray = np.zeros( [1, rowDim * colDim] )
cumulativeCount = np.zeros( [1,rowDim * colDim] )

# cumulating based on condition
for eachFile in fileList:
	fileBaseName = os.path.basename(eachFile)
	julian_date =  int(fileBaseName[5:8])
	if (julian_date >= 130 & julian_date <= 178):
		f1 = open(eachFile)
		data1 = np.fromfile( f1, dtype = np.float32 ) ; 	f1.close() ;
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

