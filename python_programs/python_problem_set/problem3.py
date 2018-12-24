#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors

# file & data
bathFileName = '/Users/md875/data/tutorial_data/bathymetry.i999x999'
fBath= open(bathFileName)
bathData = np.fromfile( fBath, dtype = np.int16 )
fBath.close()
bathArray = bathData.reshape([999,999])

chlFileName = '/Users/md875/data/tutorial_data/S1998148172338_chlor_a.f999x999'
fChl= open(chlFileName)
chlData = np.fromfile( fChl, dtype = np.float32 )
fChl.close()
chlArray = chlData.reshape([999,999])

# input min & max depth
minDep = input('\nProvide min depth: ')
maxDep = input('\nProvide max depth: ')

# selecting our satisfied regions
satisfied_regions_indices = np.where( (bathArray>=minDep) & (bathArray<=maxDep) )  #returns a matrix/array with indices that satisfies our condition

# mean and SD for selected regions
meanChl = np.nanmean(chlArray[satisfied_regions_indices])
SDChl = np.nanstd(chlArray[satisfied_regions_indices])

print meanChl
print SDChl









