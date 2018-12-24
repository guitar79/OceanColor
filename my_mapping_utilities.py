#!/usr/bin/env python

from numpy import *
from pylab import *
 

import pyproj
import pyresample as pr


from my_hdf_cdf_utilities import *
from my_general_utilities import *


import sys, os 
import subprocess  
import math 
import shutil

from matplotlib import *
import matplotlib.pyplot as plt
import matplotlib.colors
from mpl_toolkits.basemap import *

import scipy.misc


   
def compute_cntr_latlon(south,west,north,east):
#-------------------------------------------------------------------------------        

    center_lat = 0.0
     
    if west > 0.0 and east > 180.0:
        center_lon = west + abs(east - west)/2.0
    elif west > 0.0 and east < 0.0:
        delta_w = 180.0 - abs(west)
        delta_e = 180.0 - abs(east)
        
        if delta_w > delta_e:
            center_lon = west + (delta_w + delta_e)/2.0
        else: center_lon = east - (delta_w + delta_e)/2.0
    else:
        center_lon = west + abs(west - east)
   
    return center_lat, center_lon
    


def map_l2_to_cyl(l2_lon, l2_lat, l2_data, l2_resolution, map_coords):
#------------------------------------------------------------------------------- 
    
    # if lat and lon 2D arrays are not the same dimesions as the geophy data dimensions
    # then lineratly intereloplte lat lon values to the same 2d dimensions as the data...
    # ---
    ydim_data,  xdim_data    = l2_data.shape    
    ydim_latlon, xdim_latlon = l2_lon.shape   
       
    if ydim_latlon != ydim_data or xdim_latlon != xdim_data:
   
        l2_lon=scipy.misc.imresize(l2_lon, l2_data.shape, interp='bilinear',mode='F')
        l2_lat=scipy.misc.imresize(l2_lat, l2_data.shape, interp='bilinear',mode='F')
         
     
    # NOTE: l2_lon, l2_latm l2_data are 2D arrays read in from an l2 files read 
    # in the main program
    # l2_resolution is the spatial resolution of the input l2 data (in meters)
    # map_swne and map_resolution are the output map bounds (in degree)
    
    #NOTE.....May have to mask lat on where data bad... ALSO Check Fill_Values...
        
    #---------------------------------------------------------------------------
    # convert delta lat and delta lon of desired output map into the number of 
    # out map pixels (xdim, ydim) needed so accomodate the spatial resolution 
    # of the l2_data.  The determination of xdim, ydim is is based on the fact
    # that there are 111km per degree of lat (and lon at the equator).
    # The pix_per_deg for lon well vary from 111.0 at the equator to
    # less than this at higher latitudes. The reduction goes as the cos(lat)..
            
    #note: for 10000 meter resolution pix_per_deg = 111.0
    #      for 500 meter resolution pix_per_deg =   222.0
    #      for 250 meter resolution pix_per_deg =   444.0
        
    
    north= float(map_coords.north)
    west=  float(map_coords.west) 
    south= float(map_coords.south) 
    east=  float(map_coords.east) 
    
  
    lat_0, lon_0 =  compute_cntr_latlon(south,west,north,east)  
   
    pix_per_deg =   111.0*(1000.0/l2_resolution)        
    lon_scale_fac=  math.cos((math.pi/360.0)*abs(lat_0))        
    
    ydim= math.ceil(abs(north-south)*pix_per_deg)
                          
    if west > 0.0 and east < 0.0:
        
        delta_w= 180.0-abs(west)
        delta_e= 180.0-abs(east) 
        xdim= math.ceil((delta_w + delta_e)*pix_per_deg*lon_scale_fac)

    else:
        xdim= math.ceil(abs(west-east)*pix_per_deg*lon_scale_fac)    
    
            
    #convert expected mapped lat/lon bounds from degrees to meters for an
    #cylindrical (cyl) coordinate system also known as an eqirectangular 
    #coordinate system(eqc)
    p = pyproj.Proj('+proj=eqc +lat_ts=0 +lat_0=0 +lon_0=0 +x_0=0 +y_0=0 +a=6378137 +b=6378137 +units=m')   
    west_m, south_m, = p(west, south)   
    east_m, north_m  = p(east, north)   

     
    #set up variabl names for subsequent use in the mapping call...
    area_id =   'global'  #can be any name
    area_name = 'Global'  #cn be any name
    proj_id =   'cyl'     #can be any name
    proj4_args = '+proj=eqc +lat_ts=0 +lat_0=0 +lon_0=0 +x_0=0 +y_0=0 +a=6378137 +b=6378137 +units=m'
    area_extent= (west_m, south_m, east_m, north_m)   #obtaind abve from call to pyproj.Proj
    
  
    area_def = pr.utils.get_area_def(area_id, area_name, proj_id, proj4_args, xdim, ydim, area_extent)
    swath_def = pr.geometry.SwathDefinition(l2_lon, l2_lat)
    result = pr.kd_tree.resample_nearest(swath_def, l2_data, area_def, \
                        radius_of_influence=5000, fill_value=-32767.0)      
     
     #original radius of influence 5000 works for 1km resoluton swath
     #so for 90 meter maybe try 10x high or lower 
     
    fill_locations= where(result == -32767.0)
    if len(fill_locations[0]) != 0 : result[fill_locations]= nan
     
     
    return result


 


def write_png_with_basemap(png_fname, geophys_img, product, latlon, proj_name):
#---------------------------------------------------------------------------- 
  
    # CREATES A PNG OF AN EXISTING BINARY MAPPED FILE 
    # inputs:
    # png_fname ==> output png file
    # geophys_img ==> input 2D array of values
    # product ==> products (such as chlor_a)
    # latlon ==> object of type Coords()
    # proj_name ==> projection type (such as 'RECT' or 'CYL')
   
    rootname = os.path.basename(png_fname)    
    #center_lat, center_lon = get_cntr_latlon(latlon)#not used
    rotation_deg = 0.0
    size_info = geophys_img.shape # returns (# rows, # cols)
    product = product.strip() #trim leading or trailing space
    
    global_smi_chk = '.main' in png_fname
    if global_smi_chk is True: center_lon = 0.0    
   
    if proj_name is 'Equidistant Cylindrical' or proj_name is 'RECT': proj_name = 'Cylindrical'
    
    #-------------------------------------------------------------------#
    
    # reset scale_type between successive calls to write_png_with_basemap...
    scale_type= ''
   
    ocsswroot_path= os.getenv('OCSSWROOT')
    seadas_lut_resources = ocsswroot_path + '/run/data/common/palette'   
    print ( '\nseadas_lut_resources $OCSSWROOT... see my_map_utilities (def write_my_png) for more detailes  ---> ',  seadas_lut_resources )
     
    local_resources= os.getenv('LOCAL_RESOURCES')
    print ( 'local_resources $LOCAL_RESOURCES... see my_map_utilities  (def write_my_png) ---> ',  local_resources )
                     
   
    prod_min_max_table = local_resources + '/png_min_max_settings/prod_min_max_tab_delimted_txt'
    prod_min_max_info =  get_prod_min_max(prod_min_max_table, product)
    
    
    print ( '\nprod_min_max_info -----------> ', prod_min_max_info )
    print ( 'png fname to be created -----> ', png_fname )
    

    low_limit = prod_min_max_info[1]
    upper_limit = prod_min_max_info[2]
    scale_type = prod_min_max_info[3]

  
                      
    if low_limit == '': low_limit = -999.
    else: low_limit = float(low_limit)  
    
    if upper_limit == '': upper_limit = 999.
    else: upper_limit = float(upper_limit)

                      
    if scale_type == '': scale_type = 'LIN'       
 

    

      
    # specify the color palette to use...
    if product[:3] != 'sst':
        mycmap = custom_cmap(seadas_lut_resources + '/brs_chl.pal')          
    if product[:3] == 'sst':
        mycmap = custom_cmap(seadas_lut_resources + '/brs_sst.pal')     
          
    # specify color map with NaN's as black...    
    mycmap.set_bad('k')
 
              
    
    # Set Basemap Resolution to Low (l) for global maps, 
    # high (h) for most nomral regional maps and 
    # full (f) for landsat images       
    # ---
    png_root = os.path.basename(png_fname)
    
    if global_smi_chk is True:    
        m = Basemap(projection='cyl',llcrnrlon=latlon.west, \
            llcrnrlat=latlon.south, urcrnrlon=latlon.east, \
            urcrnrlat=latlon.north, resolution = 'l')                                     
                                           
    if global_smi_chk is False and png_root[0] != 'L':    
        m = Basemap(projection='cyl',llcrnrlon=latlon.west, \
            llcrnrlat=latlon.south, urcrnrlon=latlon.east, \
            urcrnrlat=latlon.north, resolution = 'h')    
    
    if global_smi_chk is False and png_root[0] == 'L':    
        m = Basemap(projection='cyl',llcrnrlon=latlon.west, \
            llcrnrlat=latlon.south, urcrnrlon=latlon.east, \
            urcrnrlat=latlon.north, resolution = 'f')      
    
   
       
    geophys_img= flipud(geophys_img) # <----<<< flip array upside down only for mapping...         
    
    if low_limit == 0 and scale_type == 'LOG':
        print ( '\n[WARNING]!!!   =======> You are using LOG scale AND vmin=ZERO <=====' )
        print ( '\nThis will break your m.show function in the write_png_with_basemap function (with my_mapping_utilities)' )
        print ( '\nGo the python_programs/local_procesing_resources/png_min_max_settings/prod_min_max_tab_delimted_txt' )
        print ( '\nand change either the ZERO minmium or change LOG scale to LIN....' )
    
   
    # Draw Colorbar..
    if scale_type == 'LIN': m.imshow(geophys_img, cmap=mycmap, vmin=low_limit, vmax=upper_limit)                                                                     
    if scale_type == 'LOG': m.imshow(geophys_img, cmap=mycmap, vmin=low_limit, vmax=upper_limit,norm=matplotlib.colors.LogNorm()) 
    m.colorbar() 
   
   
    
    
    # Draw and fill continents..
    m.drawcoastlines()
    m.fillcontinents(color='grey',lake_color='white')
    
 
       
    # Setup and draw Lat/Lon Labels..
    delta_lon=abs(latlon.west - latlon.east)
    delta_lat=abs(latlon.north-latlon.south)  
    
    if delta_lon < 4:                      lon_space=.25
    if delta_lon >=  4 and delta_lon < 8:  lon_space=1
    if delta_lon >=  8 and delta_lon < 16: lon_space=2
    if delta_lon >= 16:                    lon_space=5  
    
    if delta_lat < 4:                      lat_space=.25
    if delta_lat >=  4 and delta_lat < 8:  lat_space=1
    if delta_lat >=  8 and delta_lat < 16: lat_space=2
    if delta_lat >= 16:                    lat_space=5  
 
    #labels = [left,right,top,bottom]
    parallels = arange(latlon.south,latlon.north,lat_space)   
    m.drawparallels(parallels,linewidth=0, labels=[True,False,False,False])

    meridians = arange(latlon.west,latlon.east,lon_space)
    m.drawmeridians(meridians,linewidth=0,labels=[False,False,True,False])

    
     
    plt.savefig(png_fname)
    
    
    #cb.remove()
    m.drawparallels(parallels,linewidth=0, labels=[False,False,False,False])
    m.drawmeridians(meridians,linewidth=0,labels=[False,False,False,False])
    plt.close()
   
    
    
    ##-------------- write png info file
    png_dir = os.path.dirname(png_fname)
    png_info = open(png_dir + '/README_SCALE', 'w')
    
    png_info.write('prod: ' + product)
    png_info.write('\nMIN: ' + str(low_limit))
    png_info.write('\nMAX: ' + str(upper_limit))       
    png_info.write('\nSCALE: ' + scale_type)
    png_info.write('\nPROJECTION: ' + proj_name + '\n')
    png_info.write('\nLatitude - Longitude Limits...')
    png_info.write('\nnorth: ' + str(latlon.north))
    png_info.write('\nsouth: ' + str(latlon.south))
    png_info.write('\nwest: ' + str(latlon.west))
    png_info.write('\neast: ' + str(latlon.east) + '\n')
    png_info.write('\nNOTE ----> If you do not like the min/max colorbar scaling of the png images, they can be changed...\n' +
                    '---------------------------------------------------------------------------------------------------\n' +
                    '\nFIRST, You can use SeaDAS Interactive (GUI) mode on one of the netcdf4 (mapped or smi) files\n' +
                    'and carefully apply different scales min/max and log/lin settings to find the best settings for your batch products.\n\n' +
                    'SECOND, you can go to the following text file to set up your new min/max scale settings that are specific\n' +
                    'for your product OR for your specific geographic location or season:\n\n===> ' + prod_min_max_table + 
                    '\n\n\nNote: To add a new product min/max to the current list of products to scale, append a line to the the prod_min_max_scales_tbl.txt file\n' +
                    'that looks identical to the line above, then then change the product name to the scale info a appropriate for your product.\n' +
                    '\n\nTHIRD, Once a new scaling has been set, open the ~/python_programs/utilities/Batch_PNG_Rescale.py script and follow the instructions\n'
                    'for generating new png images with the new min/max scaling.  The new png images will replace the old png images with the old scaling.\n\n' +
                    
                    'By The Way.. if you do not like the spacing the lat/lon labels or the colorbar style see: my_mapping_utilities.py\n' +
                    'under def write_png_with_basemap: scroll down to #-Setup and draw Lat/Lon Labels.. and #-add colorbar.  You will need to\n' +
                    'have some working knowledge of mapplotbib basemap functions - google things to find more info\n\n')
    
    png_info.close()
    
 
      
                
        
def mask_from_l2flags(l2_file, l2flag_names_to_check):
#----------------------------------------------------------------------------------------------------
 
     
    sst_mask_chk= 'SST' in l2_file
    
    
    if sst_mask_chk == True:
       
       print ( '\nusing -- sst_qual -- to check for bad or missing sst data...' )
       print ( '\nNOTE:  this is standard obpg qualitely level from l2 to binned.' )
       print ( '       To increase or decrease the qual level, go to mask_from_l2flags' )
       print ( '       found within my_mapping_utilities...\n' )
       
       qual_img= read_hdf_prod(l2_file, 'qual_sst') 
       mask_img= ones(qual_img.shape)
       
       bad_data_locations= where((qual_img == -1) | (qual_img > 2))  # See: http://oceancolor.gsfc.nasa.gov/forum/oceancolor/topic_show.pl?dln=13190;pid=26263
       mask_img[bad_data_locations[0],bad_data_locations[1]]=nan
       
    
    
    if sst_mask_chk == False:
        
        # Begin with an array mask_img of all ones and then whenever there is an 
        # occurance of a tripped flag in the l2flag_img, the mask_img at that location
        # is set to zero.  In the end the mask has ones where no flag tripped and
        # zeros where is tripped....
       
        print ( 'l2flag_names_to_check ----->  ', l2flag_names_to_check )
        print ( '\n' )
    
        l2flag_img= read_hdf_prod(l2_file, 'l2_flags')
    
        flags_to_check= l2flag_names_to_2lflag_number(l2_file, l2flag_names_to_check)    
        x,y =       l2flag_img.shape
        mask_img =   ones((x,y))
        ### 2to3
        #full_flags = arange(32,dtype=int32) + 1L
        #deciml_flag =  [2L]**(full_flags - 1L)
        full_flags = arange(32,dtype=int32) + 1
        deciml_flag =  [2]**(full_flags - 1)

    
   
        for i in range(len(flags_to_check)):   
         
            flags_on=    where(full_flags == flags_to_check[i])
            dmflag=      deciml_flag[flags_on]                  
            match_flag=  where((l2flag_img & dmflag[0]) == dmflag[0])    
      
            if len(l2flag_img[match_flag]) != 0:
                mask_img[match_flag]= nan         
  
    return mask_img
      
 
 
 


def l2flag_names_to_2lflag_number(l2_file, l2flag_names_used):    
#--------------------------------------------------------------------------------------------------
 
    
   # flag_names= asarray(['ATMFAIL', 'LAND', 'PRODWARN', 'HIGLINT', 'HILT', 'HISATZEN', 'COASTZ', 'SPARE', 'STRAYLIGHT', \
   #             'CLDICE', 'COCCOLITH', 'TURBIDW', 'HISOLZEN', 'SPARE', 'LOWLW', 'CHLFAIL', 'NAVWARN', 'ABSAER', 'SPARE', 'MAXAERITER', \
   #             'MODGLINT', 'CHLWARN', 'ATMWARN', 'SPARE', 'SEAICE', 'NAVFAIL', 'FILTER', 'SPARE', 'BOWTIEDEL', 'HIPOL', 'PRODFAIL', 'SPARE'],dtype='|S8') 
      
    
    flag_names=  get_l2hdf_full_l2flags_names(l2_file)              
    flag_number= arange(1,len(flag_names)+1)  

    
    #Convert one long strong of flag names into a python "List" and then vector of names
    l2flag_names_used_list= l2flag_names_used.split(',')
    l2flag_names_used_vec=  asarray(l2flag_names_used_list, dtype='|S8') #vector form
    


    # OBPG made a name change from 'HISOLZ' to 'HISOLZEN' ---  I think this fix is no longer needed... bcm 4/19/2016
    # bad_solar_zenith_name= where(l2flag_names_vec == 'HISOLZ')
    # if l2flag_names_vec[bad_solar_zenith_name] != 0: l2flag_names_vec[bad_solar_zenith_name]= 'HISOLZEN'
         

    l2flag_number_used= arange(len(l2flag_names_used_vec))     
  
     
    
    for i in range(len(l2flag_names_used_list)): 
        
       good_locations= where(flag_names == l2flag_names_used_vec[i])                
       l2flag_number_used[i]= flag_number[good_locations]

  
    return l2flag_number_used
  

  
  
def resolution_from_sat_fname(fname):
#-------------------------------------------------------------------------     
  
# Expect input string of the form: fname= 'A20004hhmmss...'
     
    rootname = os.path.basename(fname)       
    sat_id= rootname[0]
  
    
    mer_frs_chk= 'FRS' in fname    #returns 1 if true or 0 if false   
    hkm_chk=     'HKM' in fname   
    qkm_chk=     'QKM' in fname
     
    
    if sat_id == 'S': 
         sat_resoulution = 1000.0         

                  
    if ((sat_id == 'A') or (sat_id == 'T')):
        if hkm_chk == 1:
            sat_resoulution = 500.0
        elif qkm_chk == 1:
            sat_resoulution = 250.0
        else:
           sat_resoulution = 1000.0
         
         
    if sat_id == 'M':
        if mer_frs_chk == 1:
            sat_resoulution = 300.0
        else:
            sat_resoulution = 1000.0
            
            
    if sat_id == 'V':
        sat_resoulution = 750.0
        
        
    if sat_id == 'H':        
        sat_resoulution = 90.0
           
    
    if sat_id == 'L':        
        sat_resoulution = 30.0
        
            
                    
    return sat_resoulution




def get_product_list(ifile, requested_prod):   
#-------------------------------------------------------------------------        
    # Get color products from HDF and match them to input products
    color_prod = []
            
    #get all products available to HDF
    all_prod = asarray(get_l2hdf_prod(ifile))
    
    print ( 'requested product ---> ', requested_prod )
    print ( 'all products --------> ', all_prod )
      
    if requested_prod[0] == 'all':
        color_prod = all_prod   
    
    else:
        #for each user-specified product, put in color_prod if also availible to HDF file
        for pr in requested_prod:
            good_prod_indx = where( all_prod == pr )
            if len(good_prod_indx) != 0: 
                color_prod = concatenate((color_prod, all_prod[good_prod_indx]), axis=0)
            
    # get rid of any empty products
    print ( color_prod )
    #color_prod = color_prod[where(color_prod != '')]
    color_prod = color_prod[where(color_prod != '')]
    return color_prod
    
    
     


def get_sds7_default_l2flags(sensor_id, prod_type):
#-------------------------------------------------------------------------          

    #  input:
    #  sensor_id ==> single character, beginning of satellite name, e.g. 'A', 'M', ...
    #  prod_type ==> products ('color' or 'sst')
    #
    # output: 
    # comma-separated string of flags ('ATMFAIL,LAND,HILT,HISATZEN,STRAYLIGHT,CLDICE,...')
     
    sensor_subdir= ['czcs', 'hmodisa', 'hmodist', 'meris', 'modis', 'modisa', 'modist', 'ocrvc', 'seawifs', 'viirsn', 'hico']
    
    if sensor_id == 'S' and prod_type == 'color':
        #2to3
        #fname = '$OCSSWROOT/run/data/seawifs/l2bin_defaults.par'
        fname = '$OCSSWROOT/share/seawifs/l2bin_defaults.par'
    if sensor_id == 'A' and prod_type == 'color':
         fname = '$OCSSWROOT/run/data/hmodisa/l2bin_defaults.par'
    if sensor_id == 'T' and prod_type == 'color':
         fname = '$OCSSWROOT/run/data/hmodist/l2bin_defaults.par'
    if sensor_id == 'M' and prod_type == 'color':
         fname = '$OCSSWROOT/run/data/meris/l2bin_defaults.par'
    if sensor_id == 'V' and prod_type == 'color':
         fname = '$OCSSWROOT/run/data/viirsn/l2bin_defaults.par'
    if sensor_id == 'C' and prod_type == 'color':
         fname = '$OCSSWROOT/run/data/czcs/l2bin_defaults.par'
    if sensor_id == 'H' and prod_type == 'color':
         fname = '$OCSSWROOT/run/data/hico/l2bin_defaults.par'
    
    
    if sensor_id == 'L' and prod_type == 'color':
         #fname = '$OCSSWROOT/run/data/oli/l2bin_defaults.par'
         #fname = '$OCSSWROOT/run/data/hmodisa/l2bin_defaults.par'
         fname = '$OCSSWROOT/share/modis/l2bin_defaults.par'
         
        

    if sensor_id == 'A' and prod_type == 'sst':
        ##2to3 
        #fname = '$OCSSWROOT/run/data/modisa/l2bin_defaults_SST.par'
         fname = '$OCSSWROOT/share/modis/l2bin_defaults_SST.par'
         
    if sensor_id == 'T' and prod_type == 'sst':
        #fname = '$OCSSWROOT/run/data/modist/l2bin_defaults_SST.par'
        fname = '$OCSSWROOT/share/modis/l2bin_defaults_SST.par'

    if sensor_id == 'V' and prod_type == 'sst':
         fname = 'cc'  

    result = os.popen('more ' + fname + ' | grep  flaguse').read()

    result = result.replace("\n", "")

    flags = result.split('=')

    if prod_type == 'sst': 
        flags[1] = flags[1] + ',CLDICE,NAVWARN,SEAICE,SSTWARN,SSTFAIL' 
     
    #return flags[1]
    return flags[0]
    
    

    
#
# input: 
#   ct ==> IDL color parameter file
# ouput: python color map
#
def custom_cmap(ct):
    f = open(ct)
    color_array = array([])
    for line in f.readlines():
    	color_array = append(color_array, line)
    f.close()
    
    #2to3
    #color_array = map( lambda i: map( lambda j: round(float(j)/255.0, 2), i.split() ), color_array[0:] )
    color_array = [[round(float(j)/255.0, 2) for j in i.split()] for i in color_array[0:]]
    color_int_array = arange(0.0,1.0,round(1.0/(len(color_array)-1), 3))
    color_int_array[-2] = 1.0    
    
    #2to3
    #red_array = map( lambda i: i[0], color_array)
    #green_array = map( lambda i: i[1], color_array)
    #blue_array = map( lambda i: i[2], color_array)
    red_array = [i[0] for i in color_array]
    green_array = [i[1] for i in color_array]
    blue_array = [i[2] for i in color_array]
   

    red_tup = [[color_int_array[i], red_array[i+1], red_array[i]] for i in range(0,len(color_int_array)-1)]
    green_tup = [[color_int_array[i], green_array[i+1], green_array[i]] for i in range(0,len(color_int_array)-1)]
    blue_tup = [[color_int_array[i], blue_array[i+1], blue_array[i]] for i in range(0,len(color_int_array)-1)]

   

    cdict = {'red': red_tup,		
    		 'green': green_tup,	
    		 'blue': blue_tup}
    		 

    return matplotlib.colors.LinearSegmentedColormap('colormap',cdict,256)


#
# input: 
#   fname ==> text file with table of values (~/idl_pros/dly_wkl_mon_qcmasked_pros/png_min_max_settings/prod_min_max_tab_delimted_txt)
#   prod ==> product (such as chlor_a)
# 
# output: vector of values [product, minimum, maximum, scale type]
#
def get_prod_min_max(fname, prod):
    print ( '\ngetting product min/max for color scaling png output using:' )
    print ( 'min/max table: ', fname )
    print ( 'searching for product: ', prod )
    print ( '\n' )

    f = open(fname)
    l = f.readlines()
    m = l[0].split('\r')
    try:
        #ind = next(idx for idx, string in enumerate(m) if prod in string)
        ind = next(idx for idx, string in enumerate(m) if prod in string and len(string.split()[0]) == len(prod))
        a = m[ind].split()      
    except:
        a = []
    if len(a) == 0:
        print ( '#'*50 )
        print ( ' requested geophysical product not found in prod_min_max_table: ' )
        print ( ' table file ===> ', fname )
        print ( ' to add a new product open max table and insert a new product line ===USING TABS TO SEPARATE COLUMN VALUES===' )
        print ( ' new product to add ===> ', prod )
        print ( '#'*50, '\n' )
        return_vec = ['','','','']
    else:
        
        prod_vec =  a[0]
        minval_vec = a[1]
        maxval_vec = a[2]
        scaletype_vec = a[3]
        print ( 'getting min/max for =====> ', prod_vec, '\n\n' )
        return_vec= [prod_vec, minval_vec, maxval_vec, scaletype_vec]

        
    return return_vec

                            

def str_map_gen(l2_file_list, ofname, prod, proj_type, input_coords, space_res, named_flags_2check, stats_yesno):
    
      
      north= float(input_coords.north) 
      west=  float(input_coords.west) 
      south= float(input_coords.south) 
      east=  float(input_coords.east) 
      
      
      # get spatial resoluton of l2_file (in meters) from file name 
      # (all files are assumed to be the same type of file/resoltion)   
      # ---          
      swath_resolution= resolution_from_sat_fname(l2_file_list[0])
   
      lat_0, lon_0 =  compute_cntr_latlon(south,west,north,east)     
      pix_per_deg =   111.0*(1000.0/swath_resolution)        
      lon_scale_fac=  math.cos((math.pi/360.0)*abs(lat_0))        
    
      ydim= int(math.ceil(abs(north-south)*pix_per_deg))
                          
      if west > 0.0 and east < 0.0:
        
          delta_w= 180.0-abs(west)
          delta_e= 180.0-abs(east) 
          xdim= int(math.ceil((delta_w + delta_e)*pix_per_deg*lon_scale_fac))

      else: xdim= int(math.ceil(abs(west-east)*pix_per_deg*lon_scale_fac))   
    
     
      sumx=  zeros((ydim,xdim), dtype=float) 
      sumxx= zeros((ydim,xdim), dtype=float) 
      nobs=  zeros((ydim,xdim), dtype=float)
      

      for ifile in l2_file_list:
          
          swath_data=    read_hdf_prod(ifile, prod)
          slope_intercept= get_l2hdf_slope_intercept(ifile, prod)
          
          
          print ( '\n>>>>---  str map slope_intercept -----> ', slope_intercept )
          
          swath_data= swath_data*slope_intercept[0] + slope_intercept[1]  # if no scaling found, assumed: slope=1, interecept=0.
          
          swath_qcmask=  mask_from_l2flags(ifile, named_flags_2check)
          swath_data=    swath_qcmask*swath_data                          # mask has 1's for valid and NaNs where not valid...
          
          swath_lon=  read_hdf_prod(ifile,"longitude")
          swath_lat=  read_hdf_prod(ifile,"latitude")
                    
          mapped_data= map_l2_to_cyl(swath_lon, swath_lat, swath_data, swath_resolution, input_coords)    
               
          sumx  += nan_to_num(mapped_data)
          sumxx += nan_to_num(mapped_data)**2.0 
          nobs  += (~isnan(mapped_data)).astype(int)
          
      data_avg= sumx/nobs
      data_var= sumxx/nobs - data_avg**2.0
               
      if len(l2_file_list) == 1: stats_yesno = 'no'  #if only a single file (i.e., DLY), then force stats to "no"
      
                        
      write_netcdf4_map(ofname, prod, proj_type, input_coords, space_res, named_flags_2check, data_avg, data_var, nobs, stats_yesno)
          
      

     

# create png file with coastlines from an mapped (.nc) file using mapplotlib's basemap 
# ---
# create png file with coastlines from an mapped (.nc) file using mapplotlib's basemap 
# ---
def png_rescale(indir, product):
   
    
    png_dir= indir + '/png'    
    os.system('mkdir -p ' + 'png_dir')  
    fname_list= glob.glob(indir + '/*.nc')
    
    
    for ifile in fname_list:
        
        if ifile.find('.smi.') != -1:
            mappping_approach= 'binmap' #string .smi. found (it did not retrun a -1 when is made the string search)
        
        if ifile.find('.map.') != -1:
            mappping_approach= 'str_map' #string .map. found (it did not retrun a -1 when is made the string search)
    
    
        if mappping_approach == 'binmap' : prod_img = asarray(read_hdf_prod(ifile, product))            #read -smi.nc file
        if mappping_approach == 'str_map': prod_img = asarray(read_hdf_prod(ifile, product + '-mean'))  #read -map.nc file
         
         
        if mappping_approach == 'binmap': 
            slope_intercept= get_l3mapgen_slope_intercept(ifile, product)
            prod_img =   prod_img*slope_intercept[0] + slope_intercept[1]
    
    
        bad_locations = where(asarray(prod_img) == nan)   
        if bad_locations[0] is not -1:
            prod_img[bad_locations] = -32767.0
         
    
        if mappping_approach == 'binmap'  : proj_name = get_smi_projection(ifile)
        if mappping_approach == 'str_map' : proj_name = read_hdf_prod(ifile, 'map_projection')
    
        if mappping_approach == 'binmap'  : extracted_coords = get_hdf_latlon(ifile)
        
        if mappping_approach == 'str_map' : 
            map_bounds =  read_hdf_prod(ifile, 'map_bounds_swne')
            extracted_coords = map_coords.map_coords()
            extracted_coords.south= map_bounds[0]
            extracted_coords.west=  map_bounds[1] 
            extracted_coords.north= map_bounds[2] 
            extracted_coords.east=  map_bounds[3]

    
        if not os.path.exists(png_dir):
            os.makedirs(png_dir)
        png_ofile = png_dir + '/' + os.path.basename(ifile)[:-7] + '.png'
        
   
        # call to png generating function
        write_png_with_basemap(png_ofile, prod_img, product, extracted_coords, proj_name)
    
        if os.path.exists(png_ofile):
            print ( '\nwrote file ', png_ofile, '\n\n' )
        else: print ( 'could not generate png!!!' )
    
             
      
    
             
      
  
       
    

    
        
            
                
                    
                        
                            
                                
                                    
                                        
                                            
                                                
                                                    
                                                        
                                                            
                                                                
                                                                    
                                                                        
                                                                            
                                                                                
                                                                                    
                                                                                        
                                                                                                
