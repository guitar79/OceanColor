#! /usr/bin/env python

import glob
import numpy as np
from subprocess import *
import sys, os
import datetime

sys.dont_write_bytecode = True
sys.path.insert(0, 'utilities')

#import general_utilities

from my_hdf_cdf_utilities import *
import my_general_utilities as general_utilities

#------------------------------------------------------------------------------   
# lists for standard products produced in the L1 to L2 processing step...
# note that is a specific prod_list is prescribed in the Batch_Proce.py then 
# the OC_suite below is ignored and processing will produce the products
# precribed in Batch_Proce.py
#------------------------------------------------------------------------------


def get_oc_suite(satellite_name):
  
    if satellite_name == 'modis':
        prod_list = 'aot_869,angstrom,Rrs_412,Rrs_443,Rrs_488,Rrs_531,Rrs_547,Rrs_667,chlor_a,Kd_490,pic,poc,cdom_index,ipar,nflh'           
    
    elif satellite_name == 'seawifs':
        prod_list = 'aot_865,angstrom,Rrs_412,Rrs_443,Rrs_490,Rrs_510,Rrs_555,Rrs_670,chlor_a,Kd_490,pic,poc,cdom_index,par'
               
    elif satellite_name == 'meris':
        prod_list = 'aot_865,angstrom,Rrs_413,Rrs_443,Rrs_490,Rrs_510,Rrs_560,Rrs_620,Rrs_665,Rrs_681,chlor_a,Kd_490'            
      
    elif satellite_name == 'viirs':
        prod_list = 'chlor_a,Kd_490,aot_862,angstrom,Rrs_410,Rrs_443,Rrs_486,Rrs_551,Rrs_671,pic,poc,par'            
    
    elif satellite_name == 'oli':
        prod_list = 'chlor_a,Kd_490,Rrs_443,Rrs_482,Rrs_561,Rrs_655'      

    return prod_list




#-------------------------------
#   MODIS PROCESSING STEPS
#-------------------------------

def modis_level12(file_name, root_name, prod_list, prod_list_sst, color_l2_file_fname, sst_l2_file_fname, aerosol_corr_type, hires):
    
    
    # -----------------------------------------------------------------------------------------       
    # generate geolocation file (modis_GEO.py)
    # note: cycle up to 5 times because of bug in modis_GEO.py for Mac Mavericks OS
    # -----------------------------------------------------------------------------------------
   
    print ('\n >=====>  generating geolocation file from modis L1A standard resolution bands...')
    
    fname_geo = root_name[0:14] + '.GEO'
    call('modis_GEO.py -v -o ' + fname_geo + ' ' + file_name, shell=True)
    if (not os.path.exists(fname_geo)):
       call('modis_GEO.py -v --refreshDB -o ' + fname_geo + ' ' + file_name, shell=True)                  
    if (not os.path.exists(fname_geo)):  
        seadas_home= os.environ['OCSSWROOT']
        call('rm %s/run/var/ancillary_data.db' % seadas_home , shell=True)       
        call('modis_GEO.py -v -o ' + fname_geo + ' ' + file_name, shell=True)

    # -----------------------------------------------------------------------------------------    
    # generate L1B from L1A (modis_L1B.py)
    # note: cycle up to 5 times because of bug in modis_L1B.py for Mac Mavericks OS
    # -----------------------------------------------------------------------------------------
    
    print ('\n >=====>  generating Level-1B file from Level-1A file and Geolocation File...')
    
    if hires == 'off':
    
        fname_l1b = root_name[0:14] + '.L1B_LAC'
        cnt = 0
        while (not os.path.exists(fname_l1b)) and cnt <= 5:
            call('modis_L1B.py -v -o ' + fname_l1b + ' ' + file_name + ' ' + fname_geo, shell=True)
            cnt+=1
        
    elif hires == 'on':

         fname_l1b =     root_name[0:14] + '.L1B_LAC'
         fname_l1b_hkm = root_name[0:14] + '.L1B_HKM'
         fname_l1b_qkm = root_name[0:14] + '.L1B_QKM'
         cnt = 0
         while (not os.path.exists(fname_l1b)) and cnt <= 5:
             call('modis_L1B.py -v -o ' + fname_l1b + ' ' + file_name + ' ' + fname_geo + ' -k ' + fname_l1b_hkm + ' -q ' + fname_l1b_qkm, shell=True)
             cnt+=1
    
     
      

   
    # -----------------------------------------------------------------------------------------
    # call l2gen c-compiled function for modis
    # -----------------------------------------------------------------------------------------
    
    
    # but first get ancillary data (getanc.py)
    print ('\n >=====>  checking for best ancillary data (Met and Ozone) locally and retrieving from web if needed...')
    call('getanc.py -v ' + file_name, shell=True)
    fname_ancil_list = root_name + '.anc'    
    
    
    if prod_list_sst != 'none':
       
        print ('\n >=====> generating level-2 OC data with SST data using l2gen...')
    
        call(['l2gen',
              'ifile='      + fname_l1b,
              'ofile1='     + color_l2_file_fname,
              'l2prod1='    + prod_list,             
              'ofile2='     + sst_l2_file_fname,
              'l2prod2='    + prod_list_sst,         
              'geofile='    + fname_geo,
              'oformat='    + 'netcdf4',
              'par=' + fname_ancil_list,
              'resolution=' + '1000',
              'proc_sst='   + '1',
              'aer_opt='    + aerosol_corr_type])          
 
    if prod_list_sst == 'none':
        
         print ('\n >=====> generating level-2 OC data without SST using l2gen...')
        
         call(['l2gen',
               'ifile='      + fname_l1b,
               'ofile1='     + color_l2_file_fname,
               'l2prod1='    + prod_list,             
               'geofile='    + fname_geo,
               'oformat='    + 'netcdf4',
               'par=' + fname_ancil_list,
               'resolution=' + '1000',
               'proc_sst='   + '0',
               'aer_opt='    + aerosol_corr_type])                            
                                      
     
    # if high-resolution is specified
    if hires == 'on':
        
        print ('\n >=====>  generating HIRES level-2 data from level-1 data using l2gen...\n')
        
        call(['l2gen',
             'ifile=' + fname_l1b,
             'ofile1=' + color_l2_file_fname + '_HKM',
             'l2prod1=' + 'chl_oc2',
             'geofile=' + fname_geo,
             'oformat='    + 'netcdf4',
             'par=' + fname_ancil_list,
             'resolution=' + '500',
             'aer_opt=' + aerosol_corr_type,
             'ctl_pt_incr=' + '1'])
             
        call(['l2gen',
             'ifile=' + fname_l1b,
             'ofile1=' + color_l2_file_fname + '_QKM',
             'l2prod1=' + 'Rrs_645,Rrs_859',
             'geofile=' + fname_geo,
             'oformat='    + 'netcdf4',
             'par=' + fname_ancil_list,
             'resolution=' + '250',
             'aer_opt=' + aerosol_corr_type,
             'ctl_pt_incr=' + '1'])





# -----------------------------------------------------------------------------------------
#   call l2gen for seawifs/meris/hico PROCESSING STEPS
# -----------------------------------------------------------------------------------------

def seawifs_meris_hico_level12(file_name, root_name, prod_list, color_l2_file_fname):
     
    
    # get ancillary data (getanc.py)
    print ('\n >=====>  checking for best ancillary data (Met and Ozone) locally and retrieving from web if needed...')
   
    if root_name[0] == 'H': 
        os.system('getanc.py -v -s ' + root_name[1:14]) 
        fname_ancil_list = root_name[1:14] + '.anc'
    if root_name[0] != 'H':
        os.system('getanc.py -v ' + file_name)    
        fname_ancil_list = root_name + '.anc'
    
    print ('\n >=====> generating level-2 OC data from level-1 data using l2gen...')
      
    call(['l2gen',
          'ifile='   + file_name,
          'ofile='   + color_l2_file_fname,
          'l2prod1=' + prod_list,
          'par=' + fname_ancil_list])




# -----------------------------------------------------------------------------------------
#   call l2gen for VIIRS PROCESSING STEPS
# -----------------------------------------------------------------------------------------


def viirs_level12(ifile, root_name, prod_list, color_l2_file_fname, prod_list_sst, sst_l2_file_fname): 
   
    
    # located the GEO file in the l1a_dir and if not found, wget the file from the OBPG web site...
    fname_geo_root=  root_name[0:14] + '.GEO-M_SNPP.nc'
    
    # Check to see if GEO file exists in the l1a_dir...
    ifile_dir = os.path.dirname(ifile)
    if os.path.isfile(ifile_dir + '/' + fname_geo_root): 
        fname_geo = ifile_dir + '/' + fname_geo_root
        print ('\n\nUsing VIIRS GEO File from the L1A Directory... ', fname_geo)
    else:                    
        # Get the GEO file and place it in the l1a_dir... 
        url_file= 'http://oceandata.sci.gsfc.nasa.gov/cgi/getfile/' + fname_geo_root    
        cmd= 'wget ' + url_file           
        process = Popen(cmd, shell=True)
        process.wait()       
        cmd= 'mv ' + fname_geo_root + ' ' + ifile_dir
        subprocess.Popen(cmd, shell=True) 
        fname_geo =   ifile_dir + '/' +  fname_geo_root         
        print ('Using VIIRS GEO File DownLoaded from the OPPG Website... ', fname_geo)
    
    exit()
                                                                        
    
    # get ancillary data...
    print ('\n >=====>  checking for best ancillary data (Met and Ozone) locally and retrieving from web if needed...')
    call('getanc.py -v ' + ifile, shell=True)
    fname_ancil_list = root_name + '.anc'    
    
 
    
    if prod_list_sst != 'none':
       
        print ('\n >=====> generating level-2 OC data with SST data using l2gen...')
    
        call(['l2gen',
              'ifile='      + ifile,
              'ofile1='     + color_l2_file_fname,
              'l2prod1='    + prod_list,             
              'ofile2='     + sst_l2_file_fname,
              'l2prod2='    + prod_list_sst,         
              'geofile='    + fname_geo,
              'oformat='    + 'netcdf4',
              'par=' + fname_ancil_list,
              'proc_sst='   + '1'])          
 
    if prod_list_sst == 'none':
        
         print ('\n >=====> generating level-2 OC data without SST using l2gen...')
        
         call(['l2gen',
               'ifile='      + ifile,
               'ofile1='     + color_l2_file_fname,
               'l2prod1='    + prod_list,             
               'geofile='    + fname_geo,
               'oformat='    + 'netcdf4',
               'par=' + fname_ancil_list,
               'proc_sst='   + '0'])                            



# -----------------------------------------------------------------------------------------
#   call l2gen for Landsat 8 OLI PROCESSING STEPS  
# -----------------------------------------------------------------------------------------

def oli_level12(file_name, prod_list, color_l2_file_fname):
   
    
    # Each viirs scene (file_name) comes as a direcotry of several files with the first  
    # band and the geo file being needed for input to l2gen.  The other companion bands 
    # in the directory are found by l2gen automatically...
    # ---
    #PJB delete: first_band_root = file_name.split('/')[-1] # last directory in path
    first_band = file_name #PJB: delete the following; could also just change function call to first_band: glob.glob(file_name + '/' + '*_MTL.txt')[0] # input file for l2gen acting on landsat...
       
    #ancillary list
    os.system('getanc.py -v ' + first_band)
    fname_ancil_list = first_band.split('/')[-1] + '.anc' 
    
       
    print ('\n >=====> generating level-2 OC data from level-1 oli data using l2gen...')
       
    call(['l2gen',
          'ifile='   + first_band,
            'ofile1='  + color_l2_file_fname,
            'l2prod1=' + prod_list,             
            'par=' + fname_ancil_list,
            'resolution=' + '30'])
  
 
# -----------------------------------------------------------------------------------------
#       batch processing
# -----------------------------------------------------------------------------------------

def batch_proc_L12(l1a_dir, l2_dir='not_specified', prod_list='OC_suite', prod_list_sst='none', swir_onoff='off', hires='off'):

    # make sure directories are right (/ and ~)
    l1a_dir = general_utilities.path_reformat(l1a_dir)
    l2_dir = general_utilities.path_reformat(l2_dir)
    
     
    # untar any tar files if necessary.... 
    # Large data orders from OBPR come as tared directories that contain a lot of 
    # individually compressed data files.  This function extracts all of compressed 
    # data files from alll of the tared directories found in l1a_dir
    # ---
    if len(glob.glob(l1a_dir + '/*.tar')) != 0: general_utilities.untar_ocweb(l1a_dir)
        
     
    # decompress files if necessary
    while any([general_utilities.is_compressed(fi) for fi in glob.glob(l1a_dir + '/*')]):
        for fi in glob.glob(l1a_dir + '/*'): general_utilities.decompress_file(fi)
    
    fname_l1a =      glob.glob(l1a_dir + '/*L1*')
   
    fname_landsat =  glob.glob(l1a_dir + '/L*MTL*')      # PJB changed to L*MTL* June 13 2016
    if len(fname_landsat) != 0: fname_l1a=fname_landsat  # if there are landsat directories (equal to level 1 nc files)
                                                         # then rename fname_landsat to fname_l1a and proceed...        
    
    fname_meris =  glob.glob(l1a_dir + '/MER*.N1*')      # if there are meris N1 files directories (equal to level 1 nc files)
    if len(fname_meris) != 0: fname_l1a=fname_meris      # then rename fname_meris to fname_l1a and proceed...  
                                                               
    
    if len(fname_l1a) == 0:
        print ('There are no L1 files in ' + l1a_dir )    # if not L1A and no Landat 8 then there are no usable level 1 data files..
        sys.exit()
   
       
    # if user doesn't specify level 2 directory, make one next to L1 data directory
    if l2_dir == 'not_specified':
        l2_dir = os.path.dirname(l1a_dir) + '/' + 'L2_files'
        
   
    #make L2 directory
    if not os.path.exists(l2_dir):
        os.makedirs(l2_dir)

    root_name = [os.path.basename(i) for i in fname_l1a]  #list of file names split from path = basename
    
                            
    if root_name[0][0] != 'L' and root_name[0][0:3] != 'MER': root_name_trim = [name[0:14] for name in root_name] #list of file names without L1A extension (non landsat) - includes MERIS Extracted L1 files
    if root_name[0][0] == 'L': root_name_trim = [name[0:16] for name in root_name] #list of file root names for landsat
     
    # If meris FRS Level1 (un-extracted N1 files only), then strip out year,mon,day and convert to julian day to from root name for l2 file name generation that is consistent with all L2 files
    if root_name[0][0:3] == 'MER': 
        syear=  [name[14:18] for name in root_name] 
        smon=   [name[18:20] for name in root_name] 
        sday=   [name[20:22] for name in root_name] 
        stime=  [name[23:29] for name in root_name] 
        root_name_trim= []
        for i in range(len(root_name)):            
            sjday= str(datetime.datetime.strptime(syear[i] + '-' + smon[i] +'-'+sday[i], '%Y-%m-%d').timetuple().tm_yday)
            root_name_trim.append('M' + syear[i] + sjday + stime[i])
           
    
    frs_chk = np.asarray(['FRS' in name for name in root_name] )
    
    if frs_chk.all() == True:     
        color_l2_file_fname = [l2_dir + '/' + name + '.L2_FRS_OC'  for name in root_name_trim] #list of level 2 OC files to be created     
    else:
        color_l2_file_fname = [l2_dir + '/' + name + '.L2_OC'  for name in root_name_trim] #list of level 2 OC files to be created 
      
    
    sst_l2_file_fname = [l2_dir + '/' + name + '.L2_SST' for name in root_name_trim] #list of level 2 SST files to be created
    
    
    
    

    # set SWIR option for later use in l2gen call...
    if swir_onoff == 'on':
        aerosol_corr_type = '-9'
    else:
        aerosol_corr_type = '-2'



        
    # identify satellite by first letter in file name and get appropriate
    # satellite name
    sat_id = root_name[0][0]
    
   
    if (sat_id == 'A' or sat_id == 'T'):
        satellite_name = 'modis'
    
    elif sat_id == 'S':
       satellite_name = 'seawifs'  
    
    elif sat_id == 'M':
        satellite_name = 'meris' 
          
    elif sat_id == 'V': 
         satellite_name = 'viirs'  
         
    elif sat_id == 'H': 
         satellite_name = 'hico'       
    
    elif sat_id == 'L': 
         satellite_name = 'oli'  
   
  
    # if "prod_list" was set to "OC_suite" in Batch_Proc.py then get the list
    # of products to be produced in L1 to L2 processing (see very top of this module) 
    # --- 
    if prod_list == 'OC_suite': prod_list= get_oc_suite(satellite_name)
    
    
    
    # MODIS
    if 'modis' in satellite_name:
        for i in range(0,len(fname_l1a)):
            modis_level12(fname_l1a[i], root_name[i], prod_list, prod_list_sst, color_l2_file_fname[i], sst_l2_file_fname[i],  aerosol_corr_type, hires)

    # SEAWIFS, MERIS-RR, HICO
    if ('seawifs' in satellite_name or 'meris' in satellite_name or 'hico' in satellite_name):
        for i in range(0,len(fname_l1a)):
            seawifs_meris_hico_level12(fname_l1a[i], root_name[i], prod_list, color_l2_file_fname[i])

    # VIIRS 
    if 'viirs' in satellite_name:
        for i in range(0,len(fname_l1a)):
            viirs_level12(fname_l1a[i], root_name[i], prod_list, color_l2_file_fname[i], prod_list_sst, sst_l2_file_fname[i]) 
        
    # Landsat 8 
    if 'oli' in satellite_name:        
        for i in range(0,len(fname_l1a)):
            oli_level12(fname_l1a[i], prod_list, color_l2_file_fname[i])



    # clean up... 
    types = ['*.anc','*.atteph','*L1B_LAC','*L1B_HKM*','*L1B_QKM*','*.pcf']
    for t in types:
        for i in glob.glob(l1a_dir + '/' + t): os.remove(i)
        
    types = ['*.GEO*','*.anc','*.atteph','*L1B_LAC','*L1B_HKM*','*L1B_QKM*','*.pcf']  
    for t in types:
        for i in glob.glob(t): os.remove(i)






# -----------------------------------------------------------------------------------------
#  defaults
# -----------------------------------------------------------------------------------------
def main(*args):
    import getopt
    
    arg_options = ['l1a_dir=', 'l2_dir=','prod_list=', 'prod_list_sst=', 'swir_onoff=', 'hires=']
    opts, arg = getopt.getopt(args, '', arg_options)

    if len(args) == 0:
        print ('\nUsage:\n\t batch_L12.py --l1a_dir=<l1a_dir> --l2_dir=<l2_dir> --prod_list=<prod1,prod2,prod2,...> --swir_onoff=<swir_onoff>  --hires=<on/off>\n\n')
        print ('\t--'+arg_options[0][:-1]+' (required) ==> directory containing Level 1 data files')
        print ('\t--'+arg_options[1][:-1]+' (optional, default=L2_files) ==> new directory that Level 2 data files will be written to')
        print ('\t--'+arg_options[2][:-1]+' (optional, default = OC_suite) ==> products to be extracted (chlor_a, poc, cdom_index,...)')
        print ('\t\t-'+'OC_suite defaults contain Kd_490, Rrs_vvv, angstrom, cdom_index, chlor_a, par, pic, and poc')
        print ('\t--'+arg_options[3][:-1]+' (optional, default=none) ==> products to be extracted (sst and/or sst4...)' )
       #print ('\t--'+arg_options[4][:-1]+' (optional, default=off) ==> Nitrogen Dioxide transmittance bitmask selector (on/off)')
        print ('\t--'+arg_options[5][:-1]+' (optional, defualt=off) ==> aerosol mode option (on/off)')
        print ('\t--'+arg_options[6][:-1]+' (optional, default=off) ==> high resolution for modis only (on/off)')
        print ('\n\n')
    
    else:
        for option,value in opts:
            if option == '--' + arg_options[0][:-1]:
                arg1 = value
            if option == '--' + arg_options[1][:-1]:
                arg2 = value           
            if option == '--' + arg_options[2][:-1]:
                arg3 = value
            if option == '--' + arg_options[3][:-1]:
                arg4 = value   
            if option == '--' + arg_options[4][:-1]:
                arg5 = value
            if option == '--' + arg_options[5][:-1]:
                arg6 = value
            if option == '--' + arg_options[6][:-1]:
                arg7 = value

        # make some parameters optional
        if 'arg2' not in locals():
            arg2 = 'not_specified'
        if 'arg3' not in locals(): #in locals() and arg3 = 'OC_suite'
            arg3 = 'OC_suite'      #arg3 = OC_suite
        if 'arg4' not in locals():
            arg4 = 'none'
        if 'arg6' not in locals():
            arg6 = 'off'
        if 'arg7' not in locals():
            arg7 = 'off'

        batch_proc_L12(arg1, arg2, arg3, arg4, arg5, arg6)
        

# -----------------------------------------------------------------------------------------
#   Command Line
# -----------------------------------------------------------------------------------------
if __name__=='__main__':
    main(*sys.argv[1:])
