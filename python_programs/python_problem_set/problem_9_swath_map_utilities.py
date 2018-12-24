#!/usr/bin/env 

import numpy  as np
from my_hdf_cdf_utilities import *
from my_general_utilities import *

import pyproj
import pyresample as pr



#-------------------------------------------------------------------------------  
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
       
   
   
   
#------------------------------------------------------------------------------- 
def map_l2_swath(l2_lon, l2_lat, l2_data, l2_resolution, map_coords):
#------------------------------------------------------------------------------- 
    
    # if lat and lon 2D arrays are not the same dimesions as the geophy data dimensions
    # then linearly intereloplte lat lon values to the same 2d dimensions as the data...
    # ---
    ydim_data,  xdim_data    = l2_data.shape    
    ydim_latlon, xdim_latlon = l2_lon.shape   
       
    if ydim_latlon != ydim_data or xdim_latlon != xdim_data:
   
        l2_lon=scipy.misc.imresize(l2_lon, l2_data.shape, interp='bilinear',mode='F')
        l2_lat=scipy.misc.imresize(l2_lat, l2_data.shape, interp='bilinear',mode='F')
             
    
    south= map_coords[0]
    west=  map_coords[1] 
    north= map_coords[2] 
    east=  map_coords[3] 
      
    #---------------------------------------------------------------------------
    # convert delta lat and delta lon of desired output map into the number of 
    # out map pixels (xdim, ydim) needed so accomodate the spatial resolution 
    # of the l2_data.  The determination of xdim, ydim is is based on the fact
    # that there are 111km per degree of lat (and lon at the equator).
    # The pix_per_deg for lon well vary from 111.0 at the equator to
    # less than this at higher latitudes. The reduction goes as the cos(lat)..
            
    # note: for 10000 meter resolution pix_per_deg = 111.0
    #       for 500 meter resolution pix_per_deg =   222.0
    #       for 250 meter resolution pix_per_deg =   444.0
    #---------------------------------------------------------------------------     
     
 
  
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

     
    #set up variable names for subsequent use in the mapping call...
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
    #so for 90 meter maybe try 10x high or lower ???
     
     
    return result


