#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors

file1 = '/Users/md875/data/tutorial_data/S1998148172338_chlor_a.f999x999'
file2 = '/Users/md875/data/tutorial_data/S1998148172338_chlor_a.f999x999'

# open binary files & save data into arrays
f1 = open(file1)
data1 = np.fromfile(f1, dtype = np.float32)  #assumes data were previously written out to hdd as 32bit floating pt no
f1.close()

f2 = open(file2)
data2 = np.fromfile(f2 , dtype = np.float32)
f2.close()

#reshaping data into 2D matrixs
data1 = data1.reshape([999,999])   #a priori tht the data was converted from 999x999 to an array

#reshaping data into 2D matrixs
data2 = data2.reshape([999,999])   #a priori tht the data was converted from 999x999 to an array

#printing the array size
print 'array size: ', data1.shape

#averaging arrays
total_array = np.nan_to_num(data1) + np.nan_to_num(data2)  # <- nan to num replaces NANs with 0 and infinity with large nos
total_counts = (~np.isnan(data1)).astype(int) + (~np.isnan(data2)).astype(int)
avg_array = total_array/total_counts  # avg betw 2 data1 & data2

#imaginig 2D data to the monit

#define color scale
myColorScale = plt.get_cmap('spectral')
myColorScale.set_bad('k')  #set nan to black

#create & display fig to monitor
plt.figure(1)
plt.imshow(data1, cmap = myColorScale, vmin = 0.01, vmax = 20.0, norm = matplotlib.colors.LogNorm())  # <-changing the color scale to log instead of rectifying the variable
plt.colorbar()


plt.figure(2)
plt.imshow(data2, cmap = myColorScale, vmin = 0.01, vmax = 20.0, norm = matplotlib.colors.LogNorm())


plt.figure(3)
plt.imshow(avg_array, cmap = myColorScale, vmin = 0.01, vmax = 20.0, norm = matplotlib.colors.LogNorm())
plt.show()