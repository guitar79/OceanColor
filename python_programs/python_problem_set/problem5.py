#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors
import glob

# file & data
fileNameList = glob.glob('/Users/madhur/data/tutorial_data/wkly_swf1998/*')

rowNo = 999; colNo = 999;

#initializing
cumulativeDataArray = np.zeros( [1, rowNo * colNo] )
meanOfWeeks = []; weekNo = [];  index = 0;

##subsetting
# limit/extent  (draw grid)   (change this)
northLimit = 46.0;				 westLimit = -72.0;  eastLimit = -63.0; 
southLimit = 37.0; 

# automatically calculated
del_lat = abs(northLimit-southLimit);    	del_lon = abs(westLimit-eastLimit);

# lat/lon for subset img (change this)
northLimit_subset = 40.0;				 westLimit_subset = -69.0;  eastLimit_subset = -66.0; 
southLimit_subset = 38.0; 

# conversion automatically calculated
row_subset_North = (northLimit - northLimit_subset) * (rowNo-1) / del_lat ;   	col_subset_West = (westLimit_subset - westLimit) * (colNo-1) / del_lon ;  
row_subset_South = (northLimit - southLimit_subset) * (rowNo-1) / del_lat ;   	col_subset_East = (eastLimit_subset - westLimit) * (colNo-1) / del_lon ;

# looping for cumulating/summation for each file
for eachFile in fileNameList:
	index += 1; weekNo.append(index);
	f1 = open(eachFile)
	data1 = np.fromfile( f1, dtype = np.float32 )
	data1Array = data1.reshape([rowNo, colNo])
		
	# subsetting operation
	subset_img = data1Array[ row_subset_North:row_subset_South , col_subset_West:col_subset_East ]

	#averaging
	
	meanOfWeeks.append(np.nanmean(subset_img))
	
	#debugging
	#print 'subset img: ', subset_img
	#print 'meanOfWeeks', meanOfWeeks
	
#plotting
plt.plot(weekNo, meanOfWeeks)
plt.show()
plt.close()	



