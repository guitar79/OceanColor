#!/usr/bin/env python

 
import subprocess
from netCDF4 import Dataset   
from pyhdf.SD import *

import numpy as np 
import numpy.ma as ma

import map_coords


def hdf_cdf_version(ifile):
#---------------------------------------------------------------    
    p= subprocess.Popen(["file", ifile], stdout=subprocess.PIPE)
    output, err = p.communicate()
  
    version_4_check= "version 4"   in output
    version_5_check= "version 5"   in output
    version_netcdf_check= "NetCDF" in output
    
      
    if version_4_check == 1:
       ftype= 'hdf4'       
    if version_5_check == 1:
       ftype= 'hdf5'      
    if version_netcdf_check == 1:
       ftype= 'hdf5'
   ftype= 'hdf5'
          
    return ftype
    


        
def read_hdf_prod(ifile,prod):
#---------------------------------------------------------------    
    ftype= hdf_cdf_version(ifile)
        
    if ftype == 'hdf4':
       
       DATAFIELD_NAME=prod
        
       f= SD(ifile,SDC.READ)     
       d1 = f.select(DATAFIELD_NAME)
        
       data= d1[:,:]    
       
       d1.endaccess()
       f.end()
        
       return data
       
           
    if ftype == 'hdf5':
       
             
        f = Dataset(ifile, 'r')        
        f.set_auto_maskandscale(False)  # NOTE: setting this to False turnes OFF the automatic (on the fly) application of scale_factor and off_set when netdcf4 data are read
                                        # scale_factor and off_set when netdcf4 data are read in.  This was done to ensure consistent application of slope intercerpt 
                                        # which now will have to be done manually after reading things into the main program not matter how the netcdf data were written out.
        
        print ( '\n----------------------------------------------------------------------------------------------------' )
        print ( 'Reading netCDF4 data (using -- read_hdf_prod -- fuction) with automatic mask and scale turned OFF!!'  )
        print ( 'This means you are REQUIRED to manually apply any scale_factor and offset to the NetCDF data after' )
        print ( 'reading it in with this function...' )
        print ( '----------------------------------------------------------------------------------------------------\n' )
        
        
        group_names= f.groups.keys()  #get group names (e.g. geophtsical_data_sets or navigation)
              
    
        if len(group_names) != 0:
                
              #if smi format has one group called 'processing_control', but the variables at top level (under no group) so read as if no groups        
              if group_names[0] == 'processing_control':
                                           
                  p =  f.variables[prod]                                                                                
                  if len(p.shape) == 1:data=  p[:] 
                  if len(p.shape) == 2:data=  p[:,:] 
                  if len(p.shape) == 3:data=  p[:,:,:]   
                  if len(p.shape) == 4:data=  p[:,:,:,:]                    
            
                  f.close()             
                  return data 
         
              
              #for all other regular l2 netcdf4 files with groups and varibles under groups
              else:          
                 
                  for grp_name in group_names:                               
                      var_name= f.groups[grp_name].variables.keys()  #get names of objects within each group (e,g. chlor_a or longitude)     
                      for vn in var_name:                           
                          if vn == prod:                             
                              
                              p =  f.groups[grp_name].variables[vn]                       
                              if len(p.shape) == 1:data=  p[:] 
                              if len(p.shape) == 2:data=  p[:,:] 
                              if len(p.shape) == 3:data=  p[:,:,:]    
                              if len(p.shape) == 4:data=  p[:,:,:,:]             
                                         
                              f.close()                
                              return data 

        
        # netcdf4 file has NO groups so variable are at top level and not under a group...  
        else:
            var_name= f.variables.keys()           
            
            p =  f.variables[prod]                                          
            if len(p.shape) == 1:data=  p[:] 
            if len(p.shape) == 2:data=  p[:,:] 
            if len(p.shape) == 3:data=  p[:,:,:]    
            if len(p.shape) == 4:data=  p[:,:,:,:]    
            
            f.close()                            
            return data
            
    
        
    
def hdf_prod_info(ifile):
#---------------------------------------------------------------    
    ftype= hdf_cdf_version(ifile)
    
    if ftype == 'hdf4':
    
        f= SD(ifile,SDC.READ)
        dsets= f.datasets()
        dsNames = dsets.keys()
        dsNames= sorted(dsNames)
        f.end()
        print ( '\nVariable Name:' )
        print ( '----------------------------------' )
        for vn in dsNames:
            print ( vn )
         
            
    if ftype == 'hdf5':
        
        f = Dataset(ifile, 'r') 
        group_names= f.groups.keys()
                
        
        # Check to see if netcdf4 file has groups (e.g. geophtsical_data_sets or navigation)              
        if len(group_names) != 0:
                     
            
            # If the only group (the first group) is 'processing_control' then it is an smi formated file.
            # this check on smi format is because the variables are NOT under a group and instead are at top level so skip this... 
            if group_names[0] != 'processing_control':
                         
                print ( '\nGroups and Variables within Groups' )
                print ( '----------------------------------' )
                for grp_name in group_names:        
                    print ( '\nGroup Name= ',  grp_name )
                    var_name= f.groups[grp_name].variables.keys()
                    var_name= sorted(var_name)
                    for vn in var_name:
                        print ( '   ', vn )

                f.close()
                
    
            else:
                # necdf4 file has groups, but the only group is 'processing_control' and that is a smi file and the variable are at the top level
                var_name= f.variables.keys()
                var_name= sorted(var_name)
                for vn in var_name:
                    print ( '   ', vn )
        
                f.close() 
       
        
        # netcdf4 file has NO groups then variable are at top level and not under a group...                
        else:
            var_name= f.variables.keys()
            var_name= sorted(var_name)
            for vn in var_name:
                print ( '   ', vn )
        
            f.close() 




def hdf_prod_scale(ifile, prod):
#---------------------------------------------------------------
    

    ftype= hdf_cdf_version(ifile)
    
        
    if ftype == 'hdf4':
           
        f= SD(ifile)     
        d1 = f.select(prod)    
    
        d1Attr= d1.attributes()   
        attNames= d1Attr.keys() 
        attNames.sort()
    
        print ( '\n' )
        print ( '-'*50 )
                           
        for nm in attNames:
             t= d1Attr[nm]
             print ( nm, t )
         
        d1.endaccess()
        f.end()    
    
        print ( '-'*50 )
        print ( '\n' )
           


    if ftype == 'hdf5':
        
        f = Dataset(ifile, 'r')
        group_names= f.groups.keys()    #get group names (e.g. geophtsical_data_sets or navigation)
 
        print ( '\n' )
        print ( '-'*60 )
 
        
        if len(group_names) != 0:
 
             
            # If the only group (the first group) is 'processing_control' then it is an smi formated file.
            # this check on smi format is because the variables are NOT under a group and instead are at top level so skip this... 
            if group_names[0] != 'processing_control':                
                for grp_name in group_names:        
                    var_name= f.groups[grp_name].variables.keys()  #get names of objects within each group (e,g. chlor_a or longitude)    
                    for vn in var_name:
                        if vn == prod:
                             p =  f.groups[grp_name].variables[vn]    
                             print (  p  )
                         
            
            # necdf4 file has groups, but the only group is 'processing_control' and that is a smi file and the variable are at the top level             
            else:
               var_name= f.variables.keys()  #get names of objects within each group (e,g. chlor_a or longitude)                    
               for vn in var_name:
                   if vn == prod:
                       p =  f.variables[vn]    
                       print (  p  )
           
               f.close()
     
        
        
        # netcdf4 file has NO groups then variable are at top level and not under a group...  
        else:
        
            var_name= f.variables.keys()  #get names of objects within each group (e,g. chlor_a or longitude)                    
            for vn in var_name:
                if vn == prod:
                    p =  f.variables[vn]    
                    print (  p  )   
           
                f.close()
        
        
        
        print ( '-'*60 )
        print ( '\n' )




def write_netcdf4_map(ofile, prod, proj_type, map_coords, space_res, named_flags_2check, data_avg, data_var, nobs, stats_yesno):

    ydim, xdim = data_avg.shape   #note that data_var, nobs have the same shape...
     
   

    # groups
    root_grp = Dataset(ofile, 'w', format='NETCDF4')
    fcstgrp = root_grp.createGroup('Mapped_Data_and_Params')


    # dimensions
    fcstgrp.createDimension('map_coord_dim',  4)
    fcstgrp.createDimension('resolution_dim', 1)
    fcstgrp.createDimension('projection_dim', 1)
    fcstgrp.createDimension('l2flags_dim', 1)
    
    fcstgrp.createDimension('lon_dim', xdim)
    fcstgrp.createDimension('lat_dim', ydim)


    # variables
    projections_type =  fcstgrp.createVariable('map_projection', 'S4',('projection_dim',))
    map_bounds_swne =   fcstgrp.createVariable('map_bounds_swne',  'f8', ('map_coord_dim',))
    space_resolution =  fcstgrp.createVariable('map_resolution', 'f8', ('resolution_dim',))
    l2_flags =          fcstgrp.createVariable('l2_flags_applied',   'S4', ('l2flags_dim',))    
    
    geophys_mean =      fcstgrp.createVariable(prod + '-mean', 'f8', ('lat_dim', 'lon_dim',))
    if stats_yesno == 'yes': 
        geophys_var =   fcstgrp.createVariable(prod + '-var',  'f8', ('lat_dim', 'lon_dim',))
        geophys_nobs =  fcstgrp.createVariable(prod + '-nobs', 'f8', ('lat_dim', 'lon_dim',))

    
    # data
    projections_type[:] = np.asarray([proj_type])
    map_bounds_swne[:] =  [float(map_coords.south), float(map_coords.west), float(map_coords.north), float(map_coords.east)] 
    space_resolution[:] = [float(space_res)]
    l2_flags[:] =         np.asarray([named_flags_2check])
    
    geophys_mean[:,:] =   data_avg
    if stats_yesno == 'yes': 
        geophys_var[:,:] =   data_var
        geophys_nobs[:,:] =   nobs

 
    root_grp.close()

    print ( '\n\nwrote mapped output file: ', ofile )
    print ( '\n\n' )





def write_generic_2D_netcdf4(ofile, data_2d):

    ydim, xdim = data_2d.shape    
   
    root_grp = Dataset(ofile, 'w', format='NETCDF4')
    fcstgrp = root_grp.createGroup('Data')

    fcstgrp.createDimension('x_dim', xdim)
    fcstgrp.createDimension('y_dim', ydim)
  
    geophys_data =   fcstgrp.createVariable('data', 'f8', ('y_dim', 'x_dim',))
    geophys_data[:,:] =  data_2d
   
    root_grp.close()

    print ( '\n\nwrote generic netcdf output file: ', ofile )
    print ( '\n\n' )



def get_l2hdf_prod(ifile):
#---------------------------------------------------------------    
    master_prod_list = ['angstrom','aot_862','aot_865','aot_869','cdom_index','chlor_a','ipar','Kd_490','nflh','par','pic','poc',
                        'Rrs_410','Rrs_412','Rrs_413','Rrs_443','Rrs_486','Rrs_488','Rrs_490','Rrs_510','Rrs_531','Rrs_547','Rrs_551',
                        'Rrs_555','Rrs_560','Rrs_620','Rrs_665','Rrs_667','Rrs_670','Rrs_671','Rrs_681','Rrs_645','Rrs_859','Rrs_482','Rrs_561','Rrs_655','adg_giop',
                        'adg_gsm','adg_qaa','aph_giop','aph_gsm','aph_qaa','arp','a_giop','a_gsm','a_qaa','bbp_giop','bbp_gsm','bbp_qaa',
                        'bb_giop','bb_gsm','bb_qaa','BT','calcite_2b','calcite_3b','cfe','chlor_oc2','chlor_oc3','chlor_oc4','chl_clark','chl_ocx',
                        'chl_gsm','chl_octsc','evi','flh','ipar','Kd_lee','Kd_morel','Kd_mueller','Kd_obpg','KPAR_lee','KPAR_morel','ndvi',
                        'poc_clark','poc_stramski_490','tsm_clark','Zeu_morel','Zhl_morel','Zphotic_lee','Zsd_morel', 'chl_oc2', 'sst','sst4']
                        
    prod_list = []
    
      
    ftype= hdf_cdf_version(ifile)
        
    
    if ftype == 'hdf4':
   
        f= SD(ifile,SDC.READ)
        dsets= f.datasets()
        dsNames = dsets.keys()
        dsNames= sorted(dsNames)
        f.end()
        
        full_var_name= np.asarray(dsNames)         
        
        bad_names= np.asarray(['elat','slat','clat','elon','slon','clon','k_no2','cntl_pt_cols', \
                   'k_oz','tilt','cntl_pt_rows','latitude','vcal_gain','csol_z','longitude', \
                   'vcal_offset','day','msec','wavelength','detnum','mside','year','l2_flags','F0', 'Tau_r', 'aw', 'bbw', \
                   'scan_ell','sen_mat', 'sun_ref', 'tilt_flags', 'tilt_ranges','nflag','ntilts','orb_vec','alt_ang','att_ang'])
        
        for vn in full_var_name:
            test_index= np.where(bad_names == vn)
            if len(bad_names[test_index]) == 0:
               prod_list.append(vn)                                             
         
        print ( '\nfull prod list inside of hdf4 get_l2hdf_prod...  ' )
        print ( prod_list )
        
        return prod_list  
            
              
    if ftype == 'hdf5':               
    
        f = Dataset(ifile, 'r') 
        group_names= f.groups.keys()
  
        var_name= f.groups['geophysical_data'].variables.keys()
        var_name= sorted(var_name)
        var_name= np.asarray(var_name)
         
        full_list_indx= np.where(var_name != 'l2_flags')
        prod_list= var_name[full_list_indx]     
                                                                         
        f.close()
        
        print ( '\nfull prod list inside of hdf5 get_l2hdf_prod...  ' )
        print ( prod_list )
       
        return prod_list




def get_smi_projection(file):
#---------------------------------------------------------------        
    f = Dataset(file, 'r') 
    subgroup= f.groups["processing_control"].groups 
    input_params= subgroup['input_parameters']
    
    projection= input_params.projection
    
    f.close()
    
    return projection
    



def get_l3mapgen_prod_list(ifile):
#---------------------------------------------------------------        
    
    f = Dataset(ifile) 
    group_names= f.groups.keys()
                    
    var_names= f.variables.keys()
    
    f.close() 
    
    return var_names
    
                
def get_hdf_latlon(file):                        
#---------------------------------------------------------------        
    
    extracted_coords = map_coords.map_coords()
    
    f = Dataset(file, 'r') 
    subgroup= f.groups["processing_control"].groups 
    input_params= subgroup['input_parameters']
    
    #used for smigen files...
    #extracted_coords.south= float(input_params.latsouth)
    #extracted_coords.west=  float(input_params.lonwest)
    #extracted_coords.north= float(input_params.latnorth)
    #extracted_coords.east=  float(input_params.loneast)
    
   #used for l3mapgen files...
    extracted_coords.south= float(input_params.south)
    extracted_coords.west=  float(input_params.west)
    extracted_coords.north= float(input_params.north)
    extracted_coords.east=  float(input_params.east)
    
    f.close()
    
    return extracted_coords
   



def get_l2hdf_slope_intercept(ifile, prod):
#---------------------------------------------------------------        

    slope_inter= np.asarray([1.0, 0.0])        
    
    
    ftype= hdf_cdf_version(ifile)       
    
    if ftype == 'hdf4':
       
        f= SD(ifile)     
        d1 = f.select(prod)    
    
        d1Attr= d1.attributes()   
        attNames= d1Attr.keys() 
        attNames.sort()
                            
        for nm in attNames:
            if nm == 'slope': slope_inter[0]= float(d1Attr[nm])
            if nm == 'intercept': slope_inter[1]= float(d1Attr[nm])
            
            if nm == 'scale_factor': slope_inter[0]= float(d1Attr[nm])
            if nm == 'add_offset': slope_inter[1]= float(d1Attr[nm])
                    
        d1.endaccess()
        f.end()                      
        return slope_inter       
              
                  
    if ftype == 'hdf5':  

        f = Dataset(ifile, 'r')
        group_names= f.groups.keys()    #get group names (e.g. geophtsical_data_sets or navigation)
     
 
        for grp_name in group_names:        
            var_name= f.groups[grp_name].variables.keys()  #get names of objects within each group (e,g. chlor_a or longitude)    
            for vn in var_name:
                if vn == prod:
                    p =  f.groups[grp_name].variables[vn]                    
                    try: slope_inter= np.asarray([float(p.scale_factor), float(p.add_offset)])               
                    except: print ( '\nDid not find slope intercept valules in l2 file. Using as default: slope = 1.0 and interecept = 0.0\n' )
          
        f.close()         
        return slope_inter     

 
    
def get_l3mapgen_slope_intercept(ifile, prod):
#---------------------------------------------------------------        

    slope_inter= np.asarray([1.0, 0.0])        
    
    f = Dataset(ifile, 'r')    
    
    # netcdf4 file from l3mapget has NO groups then variable are at top level and not under a group...  
    var_name= f.variables.keys()  #get names of objects within each group (e,g. chlor_a or longitude) 
    for vn in var_name:
        if vn == prod:
            p =  f.variables[vn]  
            try: slope_inter= np.asarray([float(p.scale_factor), float(p.add_offset)])               
            except: print ( '\nDid not find slope intercept valules in l2 file. Using as default: slope = 1.0 and interecept = 0.0\n' )
    
    f.close()         
    return slope_inter 
    
              
                      
            
def get_l2hdf_full_l2flags_names(ifile): 
       
    f = Dataset(ifile, 'r')
    group_names= f.groups.keys()    #get group names (e.g. geophtsical_data_sets or navigation)
     
 
    for grp_name in group_names:        
        var_name= f.groups[grp_name].variables.keys()  #get names of objects within each group (e,g. chlor_a or longitude)    
        for vn in var_name:
            if vn == 'l2_flags':   
                p =  f.groups[grp_name].variables[vn]                                
                flag_names= p.flag_meanings               
           
    f.close()
    flag_names_list= flag_names.split(' ') ##list form 
    flag_names_vec=  np.asarray(flag_names_list, dtype='|S8') #vector form  
            
    return flag_names_vec 
             
             
  
