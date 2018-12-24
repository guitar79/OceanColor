#!/usr/bin/env python
from my_hdf_cdf_utilities import *
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors

#data variable
fName = '/Users/md875/data/tutorial_data/A2010241173000.L2_LAC_OC.x.hdf'
prodName = 'chlor_a'
slope  =1.0
intercept = 0.0

mydata = read_hdf_prod(fName, prodName)

mydata = slope * mydata + intercept

print mydata.shape

mycmap = plt.get_cmap('spectral')
mycmap.set_bad('k')

plt.figure(1)

plt.imshow(mydata, cmap = mycmap, vmin = 0, vmax = 20)

plt.colorbar()
plt.show()
