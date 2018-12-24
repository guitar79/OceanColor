#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors
import os
from my_hdf_cdf_utilities import *

#read the chl file 
fname1 = '/Users/md875/data/tutorial_data/sst/2004.mon06.f774x843'
f1 = open(fname1)
data1 = np.fromfile( fname1, dtype = np.float32 )
f1.close()
