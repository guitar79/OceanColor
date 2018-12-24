#! /usr/bin/env python

import sys

sys.path.insert(0, 'utilities')
sys.dont_write_bytecode = True

import batch_L12
import batch_L23



 
#           To execute, edit the directories and variables below,
#         and call this python script from the command line by typing:
#
#                         python Batch_Proc.py
  




###############################################################################
# ---> Define the Beginning and Desired Ending Satellite Data
# ---> Processing Levels                 
###############################################################################



# Beginning Processing Level: '1' for Level-1 files, '2' for Level-2 files
#-----------------------------------------------------------------------------
Input_Level ='2'


# Final Processing Level: '2' for Level-2 files, '3' for Level-3 mapped files
#------------------------------------------------------------------------------
Final_level = '3'

 
 
 

#############################################################################
# ---> Setup the input and output directories and lat/lon bounds 
# ---> of mapped output
#############################################################################


#  Location of the Level-1 Files:
#  -------------------------------------------------------------------------
#
#l1a_dir = '/Users/byunghyun/Downloads/2003.4.week'
l1a_dir = '/home/guitar79/OC/SeaWiFS_2003-2006_l2/l1'




#  Location where Level-2 Files will be written when processed from 
#  Level-1 to Level-2.  Or... Location where Level-2 Files were downloaded 
#  from Ocean Color Web
#  -------------------------------------------------------------------------
#
#l2_dir = '/rsclass/netcdf_data_2016/seawifs_12_may27'
#l2_dir = '/Users/byunghyun/Downloads/2003.4.week/l2'
l2_dir = '/home/guitar79/OC/SeaWiFS_2003-2006_l2/2003/l2'
	


#  Location where L3map Files and PNG Files will be Written:
#  -------------------------------------------------------------------------
#
#binmap_dir = '/Users/byunghyun/Downloads/2003.4.week/l3'
binmap_dir = '/home/guitar79/OC/SeaWiFS_2003-2006_l2/2003/l3'




#  Setup Latitude and Longitude for Level-3 file generation
#  Note: Latitude-Longitude Order MUST BE: 'south, west, north, east'
#  -------------------------------------------------------------------------
#
#latlon = '32.31,126.74,45.00,135.00'           # Lat/Lon Order is S,W,N,E
latlon = '32.00,126.00,45.00,135.00'           # Lat/Lon Order is S,W,N,E



# -------------------------------------------------------------------------  
# ------------------------------------------------------------------------- 
 
 
 
 
 
 
############################################################################# 
#  ------------>  Optional L1 -> L2 Processing Variables   <---------------                              
#############################################################################
 

# LEVEL-2 OC PRODUCTS TO BE GENERATED FROM LEVEL-1 RAW DATA VIA L2GEN 
# ------------------------------------------------------------------
# The ocean color (and/or sst) list of products to be derived in L1 to L2 via
# l2gen function in seadas is set with the varable: 'prod_list_L12' 
#
# Set: prod_list_L12 = 'OC_suite' to produce a standard suite of OC products 
# (see batch_l2.py) for current suite. Otherwise set prod_list_L12 to a narrow
# list of specific OC products for example: prod_list_L12 ='chlor_a,pic,poc'
#
# For the case of modis and virrs you also have the option to set sst output
# Options Are: prod_list_L12 = 'sst' or 'none' (future plans will allow sst4)
# -----------------------------    
prod_list_L12 =    'sst'  # Option is 'OC_suite' or comma separated list of products
#prod_list_L12 = 'chlor_a,pic'
prod_list_L12_sst= 'sst'      # Option is 'sst' or 'none'



# High Resolution Product processing for "MODIS" only" --- OPTIONS: 'on' or 'off'
# ---------------------------
hires = 'off'           



# NOTE: --->  Processing hires bands is a 'slow' process if the lat/lon limits 
#exceed 3 or 4-deg lat/lon. If you are running modis be sure to extract the
#the L1A files with tight bounds at the time of ordering the data.
#
#For FRS Meris data, extraction of L1B files is not possible. In this case
#use the function:bulk_extract_meris(l1a_dir, extract_dir, extract_latlon) 
#that is contained in my_general_utilities.py to first extract the L1B. This 
#will make new extracted files (without removing originals). 
#l1a_dir= dir of unxtracted files: extract_dir= director to place new extraced
#files: extract_latlon= comma delimed string of tight lat/lon bounds swne 
# for example '36.0,-72.0, 38.0,-70.0'   Then use extracted files any normal
#Level-1 file for processing... 



# Use of Short Wave Infrared --- Options: 'on' or 'off'
# ------------------------------------------------------
swir_onoff = 'off' 
 

# There are about 1 billion different options you can choose to
# alter how l2gen processes your data and all and 99.999% are outside the
# scope of this class, but open a terminal window and type l2gen to see 
# all 1 billion options.  If some options interest you, then go to 
# batch_L12.py and manually add them to: l2gen CALL([...,...,...,]).
#-------------------------------------------------------------------------- 



###############################################################################
#   ------------>  Optional L2 -> L3 Processing  Variables  <---------------                            
###############################################################################

# l2bin -Spatial Binning Resolution (units are km) 
# Options Are: 1,2,4,9,36 
#-------------------------------
space_res = '1'       


# l3bin -Temporal Binning (Averaging Period) daily, weekly, or monthly
# Options Are: DLY, WKY, or MON  
# ------------------------------
time_period = 'MON'

# Binning or Straight Mapping Statistics Output, on or off. 
# This tunes on file output of variance and numer of pixel in binning process
# Default is no, options are 'yes' or  'no'
# -----------------------------
stats_yesno = 'no'




# Force Staraight Map...
# ----------------------------------
straight_map= 'yes'  # options are yes or no...



# Quality Flags to be used for color products. #Default is 'standard' as set 
# prescibed in the SeaDAS Installation Directory. Otherwise provide your own
# list of named flags to check (a single string of comma separated) 
# Same thing goes sst products...
#------------------------------
color_flags = 'standard'  
sst_flags = 'standard'    
                           



# If you have l2 files with more products than you want to map right now, then
# you can choose to limit which products in the l2 are mapped to Level 3. 
# Options are:   prod_list_L23 ='all'  for all products in the L2 file to be 
# mapped or else a specific list of products, for example..  
# prod_list_L23 ='all''chlor_a,pic,poc,cdom_index'
# ----------------------------
#prod_list_L23 = 'chlor_a,pic,sst'
prod_list_L23 = 'chlor_a'



# Projection for standard mapped image. 'RECT' only option for now  (= CYCL) 
# ---------------------------- 
smi_proj = 'platecarree'  #(cylindrical proj. using center lon of region.

 
  
    
###############################################################################




# ---> Don't touch the stuff below...
################################################################################
################################################################################
################################################################################

if Input_Level == '1' and Final_level == '3':
    batch_L12.batch_proc_L12(l1a_dir, l2_dir, prod_list_L12, prod_list_L12_sst, swir_onoff, hires)
    batch_L23.batch_proc_L23(l2_dir, binmap_dir, prod_list_L23, space_res,time_period, color_flags, sst_flags, latlon, smi_proj, stats_yesno, straight_map)
elif Input_Level =='1' and Final_level == '2':
   # batch_L12.batch_proc_L12(l1a_dir, l2_dir, prod_list_L12, prod_list_L12_sst, NO2_onoff, swir_onoff, hires)
    batch_L12.batch_proc_L12(l1a_dir, l2_dir, prod_list_L12, prod_list_L12_sst, swir_onoff, hires)
elif Input_Level == '2' and Final_level == '3':
    batch_L23.batch_proc_L23(l2_dir, binmap_dir, prod_list_L23, space_res,time_period, color_flags, sst_flags, latlon, smi_proj, stats_yesno, straight_map)
else:
    print '#####  Please specify different input and output levels  #####'
    sys.exit()
