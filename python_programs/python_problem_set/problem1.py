#!/usr/bin/env python
import numpy as np

latUserInput = float( input('\n Enter your lat between (25, 62): ') )
lonUserInput = float ( input('\n Enter your lon between (-48, -82): ') )

#read the SST file 
fname1 = '/Users/md875/data/tutorial_data/sst/2004.mon06.f774x843'
f1 = open(fname1)
data1 = np.fromfile( f1, dtype = np.float32 )
f1.close()

# row & col no (change this)
rowNo = 843.0 ; colNo = 774.0;

#convert data in 2D matrix
data1Array = data1.reshape( [rowNo,colNo] )

##subsetting
#bookkeeping & converting vars

# limit/extent of original img (draw grid)   (change this)
northLimit = 62.0;				 westLimit = -82.0;  eastLimit = -48.0; 
southLimit = 25.0; 
# automatically calculated
del_lat = abs(northLimit-southLimit);    	del_lon = abs(westLimit-eastLimit);

irow = ( northLimit - latUserInput ) * (rowNo-1)/del_lat			#convert input lat into row
icol = ( lonUserInput - westLimit ) * (colNo-1) / del_lon				#convert input lon into col 

SST_value = data1Array[irow, icol]

print 'The SST value for lon = %f'%lonUserInput, ' & lat = %f'%latUserInput, \
' which is equivalent x = %f'%icol, ' & y = %f'%irow, \
'is %f'%SST_value 