#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors
import os
import glob

# reading a tab delimited text files of stn info
fish_fname = '/Users/md875/data/fish_track_data.txt'
f1 = open( fish_fname, 'r' )
all_lines = f1.readlines(); 		# reads in the entire ascii file as one very long stream of continuous ascii characters as a list
f1.close()

# managing header and other lines
header_line = all_lines[0].strip()   #1st line is the header and also removes (strips out) \n
all_lines = all_lines[1:]  #remote header line from list of lines

# limit/extent of original img (draw grid)   (change this)
northLimit = 46.0;				 westLimit = -72.0;  eastLimit = -63.0; 
southLimit = 37.0; 
# automatically calculated
del_lat = abs(northLimit-southLimit);    	del_lon = abs(westLimit-eastLimit);

#file mgmt of chl
fname_chl =   '/Users/md875/data/tutorial_data/mon_swf2006/A20060602006090_mar_chlor_a_f999x999.avg'
f1_chl = open(fname_chl)
data1_chl = np.fromfile( f1_chl, dtype = np.float32 )
f1_chl.close()

# row & col no (change this)
rowNo = 999.0 ; colNo = 999.0;
#convert data in 2D matrix
data1Array_chl = data1_chl.reshape( [ rowNo, colNo ] )

# opening a new file for text output
new_fish_fname = '/Users/md875/data/fish_track_with_chl_info.txt'
f2 = open(new_fish_fname, 'w')
f2.write( header_line + '\tchl\n' )

#for loop that copy & paste from fish file to new file along with chl value
for i in range( len( all_lines ) ):
	
	#copying from fish file
	one_line = all_lines[i].strip()
	one_split_line = one_line.split('\t')	
	year = one_split_line[0]; 	jday = one_split_line[1] ; 	hr = one_split_line[2]
	lat = float(one_split_line[3]) ; lon = float(one_split_line[4]) ;
	fish_id = one_split_line[5]
	
	## retrieving corresponding values of chl	
	irow = ( northLimit - lat ) * (rowNo-1)/del_lat					#convert input lat into row
	icol = ( lon - westLimit ) * (colNo-1) / del_lon				#convert input lon into col 

	station_chl = data1Array_chl[irow, icol]
	
	#writing into the file
	f2.write( one_line + '\t' + str(station_chl) + '\n' )

f2.close()



	
	