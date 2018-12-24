#! /usr/bin/env python

import glob
import numpy as np
from subprocess import call
import sys, os

from my_mapping_utilities import *
from my_hdf_cdf_utilities import *
from my_general_utilities import *
import map_coords


sys.path.insert(0, 'utilities')




indir= '/rsclass/data/aqua_l1a_small_batch'
l1a_fnames= glob.glob(indir + '/*L1A*.hdf')             #get list of files contained in indir...
root_name = [os.path.basename(i) for i in l1a_fnames]   #list of file names split from path...

map_coords.north= 46.
map_coords.south= 37.
map_coords.east= -63.
map_coords.west= -72.

odir=  '/rsclass/data/aqua_l1a_small_batch_output_feb24.v2'  #make output directory...
if not os.path.exists(odir):
        os.makedirs(odir)

pngdir= '/rsclass/data/aqua_l1a_small_batch_output_feb24.v2/png' #make output directory...
if not os.path.exists(pngdir):
        os.makedirs(pngdir)

prod_list= 'chlor_a'
 
          
            
              
                  
for i in range(len(l1a_fnames)):
         
   
    root_name = os.path.basename(l1a_fnames[i])   
    root_name_trim= root_name[0:14]
   
    l2_fname=  odir   + '/' + root_name_trim + '.L2_OC'   
    l2b_fname= odir   + '/' + root_name_trim + '.l2b'
    l3b_fname= odir   + '/' + root_name_trim + '.l3b'
    smi_fname= odir   + '/' + root_name_trim + '.smi'  
    png_fname= pngdir + '/' + root_name_trim + '.png'
     
    
    #print '\ngenerating geolocation file from modis L1A standard resolution bands...'    
    #fname_geo = odir +'/' + root_name_trim + '.GEO'             
    #call('modis_GEO.py -v -o ' + fname_geo + ' ' + l1a_fnames[i], shell=True)
    
     
    fname_geo = odir +'/' + root_name[0:14] + '.GEO'
    call('modis_GEO.py -v -o ' + fname_geo + ' ' + l1a_fnames[i], shell=True)
    if (not os.path.exists(fname_geo)):
       call('modis_GEO.py -v --refreshDB -o ' + fname_geo + ' ' + l1a_fnames[i], shell=True)                  
    if (not os.path.exists(fname_geo)):  
        seadas_home= os.environ['OCSSWROOT']
        call('rm %s/run/var/ancillary_data.db' % seadas_home , shell=True)       
        call('modis_GEO.py -v -o ' + fname_geo + ' ' + l1a_fnames[i], shell=True)  
        
 
    print '\ngenerating Level-1B file from Level-1A file and Geolocation File...'   
    fname_l1b = odir +'/' + root_name_trim  + '.L1B_LAC'
    call('modis_L1B.py -v -o ' + fname_l1b + ' ' + l1a_fnames[i]  + ' ' + fname_geo, shell=True)
    
    
   
    
    print '\nchecking for best ancillary data (Met and Ozone) ' + \
            'on this computer or retrieving them from web if needed...'    
    call('getanc.py -v ' + l1a_fnames[i] , shell=True)
    fname_ancil_list = root_name + '.anc'    
    
    # generate level-2 file from level-1b file and geo file and ancillary files...
    # -----------------------------------------------------------------------------         
    print '\n >=====> generating level-2 OC data using l2gen...'
    print 'name of ancil file ---> ', fname_ancil_list
    print '\n'
                     
    call(['l2gen',
        'ifile='      + fname_l1b,
        'ofile1='     + l2_fname,
        'l2prod1='    + prod_list,             
        'geofile='    + fname_geo,
        'par=' + fname_ancil_list,
        'resolution=' + '1000'])                            
                                      
    os.remove(fname_ancil_list)
    os.remove(fname_geo)
    os.remove(fname_l1b)
    os.remove(indir +'/' + root_name_trim  + '.L1B_HKM.x.hdf')
    os.remove(indir +'/' + root_name_trim  + '.L1B_QKM.x.hdf')
    
              
   
    named_flags_2check= 'LAND'     
    space_res= '1'        
    product= 'chlor_a'
    
    # make bl2 file
    #----------------------------             
    call(['l2bin', 
        'infile='  + l2_fname, 
        'ofile='   + l2b_fname, 
        'l3bprod=' + product, 
        'resolve=' + str(space_res).strip(),
        'flaguse=' + named_flags_2check])
    

    # make ascii file containing bl2 file name(s) for use by l3b_gen  
    #------------------------------  
    
    
    print '\nl2b_filename in the list -----> ', l2b_fname
    print '\n'
    
    ascii_l2bin_file_list = 'l2bin_file_list.txt'
    f = open(ascii_l2bin_file_list, 'w')
    f.write(l2b_fname + '\n')
    f.close()
  
 
    #make bl3 files (l3bin)  
    #----------------------------            
    call(['l3bin', 
        'in=' + ascii_l2bin_file_list,
        'out=' + l3b_fname,
        'out_parm=' + product,                          
        'latsouth=' + str(map_coords.south),
        'latnorth=' + str(map_coords.north),
        'lonwest='  + str(map_coords.west), 
        'loneast='  + str(map_coords.east), 
        'noext='    + '1'])
    


       
    #make standard mapped (l3mapgen) ...
    #-------------------------------           
    call(['l3mapgen',
          'ifile=' + l3b_fname,
          'ofile=' + smi_fname,
          'prod=' + product,
          'deflate=' + '4',          
          'scale_type=' + 'linear',
          'projection=' + 'platecarree',
          'resolution=' + space_res + 'km',
          'interp=' + 'nearest',
          'west=' + str(map_coords.west).strip(),
          'east=' + str(map_coords.east).strip(),
          'north=' + str(map_coords.north).strip(),
          'south=' + str(map_coords.south).strip()])
    
            
                 
    
    # create png file with coastlines using mapplotlib
    #-------------------------------------------------   
    prod_img= read_hdf_prod(smi_fname, product)               
    proj_name= get_smi_projection(smi_fname)  
    extracted_coords= get_hdf_latlon(smi_fname)     
      
    
    
    
    
    
    write_png_with_basemap(png_fname, prod_img, product, map_coords, proj_name)
    
    if os.path.exists(png_fname):
        print '\nwrote file ', png_fname, '\n\n'
    else: print '\ncould not generate png!!!'
    
    os.remove(ascii_l2bin_file_list)
    os.remove(l2_fname)
    os.remove(l2b_fname)
    os.remove(l3b_fname)   
 

 
 

 