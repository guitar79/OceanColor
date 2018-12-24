#!/usr/bin/env python

import numpy as np
from my_hdf_cdf_utilities import *


from matplotlib import *
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import *
from my_hdf_cdf_utilities import *
from my_mapping_utilities import *
from my_general_utilities import *
import map_coords
from pylab import *
import math
import os




def read_ghrsst_pf_sst(fname):
#---------------------------------------------------------------
    
    root_name= os.path.basename(fname)
    pf_check= "AVHRR_Pathfinder" in root_name
    
    if pf_check == 1: 
        sst= read_hdf_prod(fname, 'sea_surface_temperature').squeeze()         
    else: sst= read_hdf_prod(fname, 'analysed_sst').squeeze() 
    
      
    sst_data=  sst.data     
    sst_mask=  sst.mask 
        
    masked_locations= np.where(sst_mask)    
    sst_data[masked_locations]=np.nan   
    
    sst_data= sst_data*.01   
    sst_data= sst_data - 273.15
             
    if pf_check == 1: sst_data= flipud(sst_data)  

    return sst_data

 
 
def read_quick_scat(fname):
#---------------------------------------------------------------    
   
    rainqc=0
   
    asc_wind_U= read_hdf_prod(fname,'asc_avg_wind_vel_u')    
    asc_wind_V= read_hdf_prod(fname,'asc_avg_wind_vel_v')
    
    des_wind_U= read_hdf_prod(fname,'des_avg_wind_vel_u')  
    des_wind_V= read_hdf_prod(fname,'des_avg_wind_vel_v') 
   
    asc_rain_flag= read_hdf_prod(fname,'asc_rain_flag')
    des_rain_flag= read_hdf_prod(fname,'des_rain_flag')
      
    
    # flag bad data due to rain at the prescribed qc level maximum
    # flage level is 7.  Best data is 7 worst data is zero ... 
    # --- set bad data to zero ---
    asc_wind_U[np.where(asc_rain_flag != rainqc)]=0
    asc_wind_V[np.where(asc_rain_flag != rainqc)]=0
        
    des_wind_U[np.where(des_rain_flag != rainqc)]=0
    des_wind_V[np.where(des_rain_flag != rainqc)]=0
            
   
    # average assending and decending wind for u and v respectively
    # slope = 0.01
    data_cnt_asc= asc_wind_U != 0
    data_cnt_des= des_wind_U != 0
    total_cnt = 1*data_cnt_asc + 1*data_cnt_des       #1* converts True,False to 1s, 0s
    wind_u= (asc_wind_U + des_wind_U)*0.01/total_cnt 
    wind_u[where(total_cnt == 0)]= np.nan   
     
    data_cnt_asc= asc_wind_V != 0
    data_cnt_des= des_wind_V != 0   
    total_cnt = 1*data_cnt_asc + 1*data_cnt_des       #1* converts True,False to 1s, 0s
    wind_v= (asc_wind_V + des_wind_V)*0.01/total_cnt 
    wind_v[where(total_cnt == 0)]= np.nan
    
   
    wind_uv = np.zeros((2,720,1440))
    wind_uv[0,:,:] =  np.roll(wind_u, 720,1) 
    wind_uv[1,:,:] =  np.roll(wind_v, 720,1)  
 
                                            #np.roll --->  shifts lon from 0-360 --->  -180 to 180
 
    return wind_uv              
 
 
 
 
def read_ssmi_speed(filename): 
#---------------------------------------------------------------   
 
    #Note: binarydata= (2,5,720,1440, dtype=int8)     
    #Full Data Array that will be read in from the hard drive
    #1440 and 720 are the xpixel and ypixel number
    #5 is the number of geophyical variable with wind = var 1
    #2 is the day/evening pass for each of the five variables
    
    #Note:  multipliers to change binary DN data to real 
    #geophysical values  xscale=np.asarray([6,.2,.3,.01,.1])                        


    f1 = open(filename, 'rb')				   	
    data = np.fromfile(f1, dtype=np.uint8) 
    f1.close()	
         
    #binarydata= data.reshape(1440,720,5,2)  
    binarydata= data.reshape(2,5,720,1440)                                  
    
    speed_am= binarydata[0,1,:,:]      #1=speed variable, 0=morning
    speed_pm= binarydata[1,1,:,:]      #1=speed variable, 1=evening
  
  
    bad_am_locations= np.where(speed_am > 250)
    bad_pm_locations= np.where(speed_pm > 250)
       
    speed_am=  speed_am*0.2   #scale factor for wind speed
    speed_pm=  speed_pm*0.2   #scale factor for wind speed
  
      
    speed_am[bad_am_locations]=  np.nan
    speed_pm[bad_pm_locations]=  np.nan
    
    total_array  = nan_to_num(speed_am) + nan_to_num(speed_pm)
    total_counts = (~isnan(speed_am)).astype(int) + (~isnan(speed_pm)).astype(int)
    speed_am_pm  = total_array/total_counts
    
    speed_am_pm = np.roll(speed_am_pm,720,1)
   
    return speed_am_pm 

 
 
def read_aviso_duacs_2014_madt_h(ifile):
#---------------------------------------------------------------    
    ssa = read_hdf_prod(ifile,'adt').squeeze()
    ssa = np.where(ssa == -2147483647, np.nan, ssa)    #set missing data to NAN
    ssa = ssa*0.00010000                               #muliply  by scale factor to get geophysical values (meters)
    ssa = np.roll(ssa,720,1)                           #Shift to -180 to 180 lon     
    print ( '\naviso adt data have units of meters...\n' )
    return ssa


def read_aviso_duacs_2014_madt_uv(ifile):    
#---------------------------------------------------------------       
    uv = np.zeros((2,720,1440))
    uv[0,:,:] = read_hdf_prod(ifile,'u')
    uv[1,:,:] = read_hdf_prod(ifile,'v')
    uv = np.where(uv == -2147483647, np.nan, uv)      #set missing data to NAN
    uv = uv*0.00010000                                #muliply  by scale factor to get geophysical values (m/s)
    uv = np.roll(uv,720,2)                            #Shift to -180 to 180 lon     
    print ( '\naviso uv data have units of meters/sec...\n' )
    return uv
   
    
def read_reynolds_oi(ifile):
#---------------------------------------------------------------     
    f1 = open(ifile, 'rb')				   	
    data = fromfile(f1, dtype=np.float32)   
    f1.close()					   	

    data=data.newbyteorder()   #swap the byte order
    sst=  data[8:64800+8]      #skip 8 hearder pieces of information
    
    input_sst_dir = os.path.dirname(ifile)
        
    land_file= input_sst_dir + '/' + 'lstags.onedeg.dat'
    f1 = open(land_file, 'rb')				   	
    land = fromfile(f1, dtype=np.int32)   
    f1.close()				 
    
    land=land.newbyteorder() 
    land_locations= np.where(land == 0)
    bad_locations1= np.where(sst < 1.8)
    bad_locations2= np.where(sst > 35.0)
    
    sst[land_locations]=np.nan
    sst[bad_locations1]=np.nan
    sst[bad_locations2]=np.nan
   
    sst=sst.reshape(180,360)
    
    return sst
     


