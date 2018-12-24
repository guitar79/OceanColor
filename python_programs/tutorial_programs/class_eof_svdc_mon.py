#!/usr/bin/env python

import glob
import numpy as np
from subprocess import call
import sys, os
import math

from my_hdf_cdf_utilities import *
from my_general_utilities import *
from my_general_read_utilities import *

import scipy.misc
import scipy.stats
from scipy.linalg import svd      

import matplotlib.pyplot as plt

from pylab import *  
        
#This is an example program that is used to show how eof analysis is performed
#on a time series of satellite images using the singular value decomposition (svd) method.
#This program can also serve as a template that can be modified for future needs of each 
#user. Note that the write statements near the bottom of this program are commented out.  
#If you want to use this program for your own work, you will want to uncomment the write 
#statements.
	

#Note: it is VERT IMPORTANT to set any missing data equal to NAN for this program to work 
 
#------------------------------------------------------------------------------------- 
 
data_dir=  '/rsclass/data'   #<------------<<<


# CHOOSE PRODUCT TO ANALYZE
#----------------------------------------------------
prod_to_analyze= 'chl'    # choices are 'sst' or 'chl'  #<------------<<<
 

  

# 1 SET YOUR OWN OUTPUT DIRECTORY PATH
# ---------------------------------------------------------------------------	
#1 -- When running this program for your own research, the following outdir
#     should be set to YOUR OWN directory path on YOUR computer. 
#---------------------------------------------------------------------------- 	 
if prod_to_analyze == 'sst':
    outdir= data_dir + '/' + 'sst_eof_output'		
    if not os.path.exists(outdir):
        os.makedirs(outdir) 			

if prod_to_analyze == 'chl':
    outdir= data_dir + '/' + 'chl_eof_output_mon'		
    if not os.path.exists(outdir):
        os.makedirs(outdir) 




# 2 READ IN DATA FILES FROM INPUT DATA DIRECTORY
# ---------------------------------------------------------------------------
if prod_to_analyze == 'sst': fname= glob.glob(data_dir + '/' + 'global_reynolds_oi/oiv2*')                		
if prod_to_analyze == 'chl': fname= glob.glob(data_dir + '/' + 'monthly_aqua_global_9km/*')  

sfname= glob.glob(data_dir + '/' + 'monthly_aqua_global_9km/S*9km.nc')  
afname= glob.glob(data_dir + '/' + 'monthly_aqua_global_9km/A*9km.nc') 
fname=  sfname+afname
nfiles= len(fname) 



# 3. SET THE 2-D DIMENSIONS OF YOUR INPUT DATA 
#----------------------------------------------------
# dimension of SST data to be read in for eof analysis...
if prod_to_analyze == 'sst':
    data_xdim=360
    data_ydim=180 

# dimension of CHL data to be read in for eof analysis...
if prod_to_analyze == 'chl':
    data_xdim=4320
    data_ydim=2160 


			
# 4. SET THE MAGNITUDE OF THE SPATIAL RESOLUTION  REDUCTION YOU 
#    WANT e.g, A FACTOR OF 10 REDUCED THE RESOLUTION BY 10-FOLD
#    this is used to reduced high spatial frequency noise and
#    speed up eof computation.
#---------------------------------------------------- 
if prod_to_analyze == 'sst': reduce_resolution_factor= 1 #for 1-deg resolution sst, no reduction is needed...
if prod_to_analyze == 'chl': reduce_resolution_factor= 10




# 5. The monthly time series of SST files goes from 1995 to 2004
#    This next statement creates a vector of numbers to use to label the
#    x-axis (time) for plotting the PC timeseries (done nearthe end of this program).  
#
#    NOTE: this vector must be of the same length as nfiles found in file_search above
#    I wrote the following formula very specifically for monthly data that began
#    in 1995. Your research case may be very different (i.e., differnent years and
#    different time intervals - like weekly or yearly instead of monthly) so you 
#    will adjust this next ine accordingly for your own purposes with other data sets.  
#	
# NOTE ---> YOU WILL NEED TO CHANGE THIS VECTOR FOR YOU OWN SPECIFIC
# DATA TIME SERIES!!!!!
#---------------------------------------------------------
if prod_to_analyze == 'sst': time_points= 1995 + np.arange(nfiles)/12.0	 # use for monthly global sst starting in 1998	
#if prod_to_analyze == 'chl': time_points= 1998 + np.arange(nfiles)       # use for annual global chl starting in 1998
if prod_to_analyze == 'chl': time_points= 1 + np.arange(nfiles)  		

				
	 
# 6. informaton about lat/lon bounds of the images and pixel resolution
#    this is used only for mapping continents etc..
#------------------------------------------------------ 	 
south=  -90.
north=   90.
west=  -180. 
east=   180.


 
	

# 7.  Read in the orginal sst files and reduce the resolution if needed to reduce  
#     fhigh requency noise AND to allow faster eof/svd computation, before storing in 
#     the data cube (time series of iamges) for subsequent eof/svd analysis.
#
# NOTE:--->  I have written a function called rebin_down_nan that acts 
#            similar to rebin expect that it ignores NANs.  In both cases (rebin and 
#            rebin_down_nan) you MUST rebin by whole factors of 2.  In the cases where
#            your own data have x-dimensions and y-dimensions that are not divisible by
#            factors of 2, you need to use CONGRID to temporarly increse/decrease the dimensions
#            to the closest whole value in size that IS divisible by factors of 2. For
#            example:  if the orginal dimensions of an image were 347 by 495 you need to
#            to use CONGIRD to increse/decrease the image size to 346 by 496 before 
#            you use rebin.  It the current SST case this was not necseeary, but it will
#            be the needed for many cases you encounter....
#  				
#  --------------------------------------------------------------------------------------------------------	
xdim=       data_xdim/reduce_resolution_factor
ydim=       data_ydim/reduce_resolution_factor
data_cube=  np.zeros((nfiles, ydim, xdim), dtype=float)	



# FOR 9KM GLOBAL MONTHLY SST USE THIS (AND COMMENT OUT GLOABL CHL BELOW)

if prod_to_analyze == 'sst':
    for i in range(len(fname)):
        geophys=read_reynolds_oi(fname[i])
        geophys=scipy.misc.imresize(geophys, (ydim,xdim), interp='bilinear',mode='F')       
        data_cube[i,:,:]= geophys



#  FOR 9KM GLOBAL ANNUAL CHL USE THIS (AND COMMENT OUT GLOABL SST ABOVE)
#----------------------------------------------------------------------------------------			
if prod_to_analyze == 'chl':
    for i in range(len(fname)):
        geophys= read_hdf_prod(fname[i],'chl_ocx')
	geophys=flipud(geophys)                      #flips image to right side up
	geophys = np.roll(geophys,data_xdim/2,1)     #shifs lon from  -180 to 180 ---> 0-360
	bad_locations1= np.where(geophys < 0.01) 
	bad_locations2= np.where(geophys > 64.0) 	
	geophys[bad_locations1[0],bad_locations1[1]]=np.nan  
	geophys[bad_locations2[0],bad_locations2[1]]=np.nan 
        #geophys=scipy.misc.imresize(geophys, (ydim,xdim), interp='bilinear',mode='F')  
	geophys= rebin_down_nan(geophys,ydim,xdim)
	data_cube[i,:,:]= geophys
	
	
    # temporally average  together avg_period number of months.
    avg_period= 3
    
    #first trim the file list to accomodate whole number of t_interval months working back from the most recent
    new_time_ponints= nfiles/avg_period
    files_to_cut= np.mod(nfiles,  avg_period)
    
    new_data_cube= np.zeros((new_time_ponints, ydim, xdim), dtype=float)	   
    for t in range(new_time_ponints):
        istart_time= files_to_cut + t*avg_period
        iend_time=   files_to_cut + t*avg_period + avg_period
        new_data_cube[t,:,:]=  np.nansum(data_cube[istart_time:iend_time,:,:],axis=0)/np.nansum(np.isfinite(data_cube[istart_time:iend_time,:,:]),axis=0)
        
    data_cube= new_data_cube
    nfiles= new_time_ponints
    time_points= 1 + np.arange(nfiles)  
#--------------------------------------------------------------------------------------------------------	
 	
 


# NOTE:
#---------------------------------------------------------------------- 
# The rest of the program does the EOF analysis on the data cube and  
# for the most part does not need to be modified 
#
# >>>---> BUT SEE SECTION BELOW ABOUT<----<<< "using the approach described in 
# the book by w. j. emery and r. e. thomson data analysis methods in physical 
# ocenography (pp328-334)..."  
#---------------------------------------------------------------------- 


# NOTE
#----------------------------------------------------------------------
# The follow code unwraps the xdim and ydim of the 3-D data cube into 
# a "very wide" 2-D array dimension. zdim (time) reamins the same and 
# you end up with a 2-D matrix (called here the eof matrix) where the first 
# dimentions is time and the second dimension is very wide and composed of 
# all the spatial locations in the orginal 2-D image (pixels now laid 
# side by side in second dimension).This 2-D matrix (eof_matrix) is what 
# the python eof analysis function takes as input.
#---------------------------------------------------------------------- 



# find pixesl locatons where all the data is present for each time point.
# if even one time point at a given location is missing, that locaton will
# be excluded from the eof analysis...
#---
flattend_array= np.sum(np.isfinite(data_cube),axis=0)  #2D Array formed after sum
good_locations= np.where(flattend_array  == nfiles)    # locations in a 2D that have all good data for all nfiles...

 
eof_matrix= np.zeros((nfiles,len(good_locations[0])), dtype=float)
for i in range(nfiles): 
    result= data_cube[i,good_locations[0],good_locations[1]]       
    eof_matrix[i,:]=result 
 

# zero out the original data cube since all the information is now
# contained in eof_matrix
# ---
data_cube=  0
geophys=    0

 



# The following removes the temporal mean for a given pixel locaton from the 
# time series for that same location to produce a time series of anomalies 
# -----------------------------------------------------  
avg_geophys= np.sum(eof_matrix,axis=0)/nfiles                       
for i in range(nfiles): eof_matrix[i,:]= eof_matrix[i,:]-avg_geophys

 


#################################################################
# OPTIONAL MODIFICATION TO THE INPUT DATA PRIOR TO EOF ANALYSIS
# using the approach described by w. j. emery and r. e. thomson
# data analysis methods in physical ocenography (pp328-334)...
#################################################################

# UNCOMMENT THE LINES BELOW TO APPLY THE DATA MODIFICATION...   
##################################################################

# 1) remove the linear trend from the time series at each location:
#  ----------------------------------------------------    
for i in range(len(good_locations)):    
    slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(np.arange(nfiles), eof_matrix[:,i])   
    eof_matrix[:,i]= eof_matrix[:,i] - (slope*np.arange(nfiles) + intercept)
 
 
# 2) NORMALIZE EACH TIME SERIES BY DIVIDING BY THE RESPECTIVE TEMPORAL STANDARD DEVIATION 
#    normalize the time series at each location by dividing by the
#    standard deviation of the data at the repective location: OPTIONAL
#    But probaby need for CHL since it it LOG Distributed...  THINK
#    of dividing by stdev and looking for veraibaility in terms of a 
#    percent change rather than an abolute change  (for CHL this means
#    doublings in the subtropics (typically lower chl) have the same 
#    weight as doubling at hight lat (typicaly higher chl)
#  ----------------------------------------------------      
if prod_to_analyze == 'chl':
    stdv_geophys=  np.asarray(map(math.sqrt, np.sum(eof_matrix**2.0,axis=0)))/(nfiles-1) 
    for i in range(nfiles): eof_matrix[i,:]=  eof_matrix[i,:]/stdv_geophys


##################################################################



#----------------------------------------------------------------------
#----------------------------------------------------------------------

print '\nBegin SVDC Function...\n'  

# Transpose is taken becaue svd function expects a matrix that has vectors 
# that are trasnposed relatove to hwo the data were orginally read into
# the eof_martix.
# ---
eof_matrix=np.transpose(eof_matrix)
U, S, V = svd(eof_matrix,full_matrices=False)	
				
	
variances= 100.0*S**2/np.sum(S**2)

print 'nfile ---> ', nfiles
print ' '
print range(nfiles)
print 'V shape ---> ', V.shape
	 
s_matrix= np.zeros((nfiles, nfiles))

print 's_matric shape ---> ', s_matrix.shape
print 'S shape ---> ', S.shape

for i in range(nfiles): s_matrix[i,i]=S[i] 
svdc_time_series= np.dot(s_matrix,V)        
 
print V.shape	
	
# Sort eof/svd index order so that each eof mode variance explained goes from
# from mode with the highes amount of orginal variance explained to mode with 
# lease amount of variance explained so later in the program the modes can be 
# out put with mode 1 having the highest explained variance...
# ---
sorted_var_index= np.argsort(-variances)  #added - sign to sort will sort larget to small
	

print '\npercent variance explained by each eof mode'
print '--------------------------------\n'
for m in range(nfiles):  print m+1, ':  ', variances[sorted_var_index[m]]
 
print '--------------------------------\n'
 
 
# write out the variance explained for each mode to a text file
# --- 
var_fname= outdir + '/' + 'variance_explained.txt' 
var_info = open(var_fname, 'w')
var_info.write('eof_mode' + '\t' + 'variance' + '\n')
for m in range(nfiles): var_info.write(str(m+1) + '\t' + str(variances[sorted_var_index[m]]) + '\n')
var_info.close()

 

	 
#----------------------------------------------------------------------	 
#NOTE---> Typically you would not try to inerpret more the first 3 eof modes, 
#         BUT if you want more modes displayed and saved to the hard drive,  
#         then increase the upper of the following for loop above 3.
#--------------------------------------------------------
for z in range(3):
   
    
    working_index= sorted_var_index[z]				 					 			
    pc_number= str(z+1)				
    temp_img=np.zeros((ydim,xdim)) +  np.nan    
    temp_img[good_locations[0],good_locations[1]]= U[:,working_index]
				
    
    pc_normalize_factor= np.nanmax(abs(temp_img))*(-1.0)     #Use +/- (1.0) to change the color coding
				                             #of weights 			 
			
    temp_img=temp_img/pc_normalize_factor   
    svdc_time_series[working_index,:]= svdc_time_series[working_index,:]*pc_normalize_factor 
	 		 
    
  # Write out the engenfunction maps and the PC timeseries plot to the outdir  
  # ---  
    mycmap = get_cmap('jet')
    mycmap.set_bad('k')	
    
    png_fname= outdir + '/' + 'spatial_engenfunction_mode_' + str(z+1) + '.png'    
    eign= plt.imshow(flipud(temp_img), cmap=mycmap, vmin=-1, vmax=1)
    cb= plt.colorbar(orientation='horizontal', shrink=.90, aspect=40, pad=.1)  
    plt.savefig(png_fname)
    plt.close()   
    
    png_fname= outdir + '/' + 'pc_timeseries_mode_' + str(z+1) + '.png'
    data=svdc_time_series[working_index,:]    
    plt.plot(time_points,data)
    plt.ylabel('Principal Component Value')
    plt.xlabel('Year')
    plt.savefig(png_fname)
    plt.close()
    
     