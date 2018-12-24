#! /usr/bin/env python

from numpy import *
from pylab import *
import sys, os
from subprocess import call
import glob
#from osgeo import gdal, gdal_array
import shutil

sys.path.insert(0, 'utilities')
sys.dont_write_bytecode = True

from my_mapping_utilities import *
from my_hdf_cdf_utilities import *
from my_general_utilities import *

            
#
# make bl2 files (l2bin)
#
def bl2_gen(filelist, l2bin_dir, product, named_flags_2check, space_res):
   
    
    print ('\nloop >>>>>>-----------> ' + product)
    print ('filelist in l2bin sst loop ------> ', filelist)
        
    # take the first 14 characters of each file name and append bin suffix
    l2bin_filelist = [l2bin_dir + '/' + os.path.basename(i)[0:14] + '.bl2bin' for i in filelist]
    
    if not os.path.exists(l2bin_dir):
        os.makedirs(l2bin_dir)
    
    if product != 'sst':
        for j in range(0,len(filelist)):
    
            print ('\n===============> running l2bin binary for color <==================')
            print ('flags to check >>>>>>----------->' +  named_flags_2check + '\n')
            call(['l2bin', 
                  'infile='  + filelist[j], 
                  'ofile='   + l2bin_filelist[j], 
                  'l3bprod=' + product, 
                  'resolve=' + str(space_res).strip(),
                  'prodtype=' + 'regional',
                  'flaguse=' + named_flags_2check])
    
    if product == 'sst':
        
        
        
        for j in range(0,len(filelist)):
    
            print ('\n===============> running l2bin binary for sst <==================' )
            print ('max qual level used = 2    |i.e., exclude -1 = missing and 3, 4 = bad qual...\n' )
            call(['l2bin', 
                  'infile='  + filelist[j], 
                  'ofile='   + l2bin_filelist[j], 
                  'l3bprod=' + product, 
                  'resolve=' + str(space_res).strip(),
                  'prodtype=' + 'regional'                  
                  'qual_max=' + '2'])
 
#'flaguse=' + named_flags_2check,


#
# make ascii file of bl2 files
#
def ascii_gen(bl2dir, filelist):
    
    ascii_file_list = bl2dir + '/' + 'ascii_bl2_list.txt'

    f = open(ascii_file_list, 'w')
    for file in filelist:
        if os.path.exists(bl2dir + '/' + file + '.bl2bin'):
            f.write(bl2dir + '/' + file + '.bl2bin' + '\n')

    f.close()

    return ascii_file_list


#
# make bl3 files (l3bin)
# 
def bl3_gen(product, ascii_file_list, input_coords, bl3_fname):
    print ('\n====> running l3bin binary <=======\n' )
    print ('ascii_list is : ', ascii_file_list )
    
    print ( 'Does file exist?' )
    print ( os.path.exists(ascii_file_list) )
    
    if not os.path.exists(os.path.dirname(bl3_fname)):
        os.makedirs(os.path.dirname(bl3_fname))
    
    call(['l3bin', 
        'in=' + ascii_file_list,
        'out=' + bl3_fname,
        'out_parm=' + product,                       
        'latsouth=' + str(input_coords.south),
        'latnorth=' + str(input_coords.north),
        'lonwest='  + str(input_coords.west), 
        'loneast='  + str(input_coords.east), 
        'noext='    + '1'])
    
    os.remove(ascii_file_list)




# old way to make mapped files (smigen)  -- discontinued jan 2016...
# ---
def smi_gen(product, meas_names, meas_vec, bl3_fname, out_file, smi_proj, space_res, input_coords):
    
    
    print ( '\n====> running smigen binary for color <=======' )
    print ( 'stats = ', meas_names, '\n'  )
    
    if not os.path.exists(os.path.dirname(out_file)):
        os.makedirs(os.path.dirname(out_file))        
         
    call(['smigen',
          'ifile=' + bl3_fname,
          'ofile=' + out_file,
          'oformat=' + 'netcdf4',
          'prod=' + product,
          'meas=' + meas_vec,
          'precision=' + 'F',
          'stype=' + '1',
          'projection=' + smi_proj,
          'resolution=' + space_res + 'km',
          'lonwest=' + str(input_coords.west).strip(),
          'loneast=' + str(input_coords.east).strip(),
          'latnorth=' + str(input_coords.north).strip(),
          'latsouth=' + str(input_coords.south).strip()])
    
    if os.path.exists(out_file):
        print ( '\nwrote file:', out_file )
    else: print ( '\ncould not generate smi file!!!' )





# new way to make mapped files -- began jan 2016...
# ---
def l3map_gen(product_str, meas_names, bl3_fname, out_file, smi_proj, space_res, input_coords):

    print ( '\n====> running l3mapgen binary <=======' )
    print ( 'stats = ', meas_names, '\n'  )
    
    if not os.path.exists(os.path.dirname(out_file)):
        os.makedirs(os.path.dirname(out_file))
                                    
       
    call(['l3mapgen',
          'ifile=' + bl3_fname,
          'ofile=' + out_file,
          'product=' + product_str,
          'deflate=' + '4',          
          'scale_type=' + 'linear',
          'projection=' + smi_proj,
          'resolution=' + space_res + 'km',
          'interp=' + 'nearest',
          'west=' + str(input_coords.west).strip(),
          'east=' + str(input_coords.east).strip(),
          'north=' + str(input_coords.north).strip(),
          'south=' + str(input_coords.south).strip()])
    
    if os.path.exists(out_file):
        print ( '\nwrote file:', out_file    )                       
    else: print ( '\ncould not generate l3mapgen file!!!' )



# create png file with coastlines from an mapped (.nc) file using mapplotlib's basemap 
# ---
def png_gen(ifile, png_dir, product, meas_names, mappping_approach):
   
    print ( '\n====> running png_gen <=======\n'    )
    print ( 'product for which png is being made...  ', product, '\n')
    
    if product == 'chl_ocx':        
               
        # In going from l3bin to l3mapgen output, seadas (internal to l3mapgen' converts the output
        # product name from chl_ocx to either chl_oc3 (modis) or chl_oc4 (seadas). This next chunk of
        # code changes the product name from chl_ocx to the new name that l3mapgen made when outputting
        # this product to the mapped.nc file...
        # ---
        print ( '\nNOTE --> encountered product -- chl_ocx --' )
        print ( '         converting this product name to either chl_oc3 or chl_oc4 for modis or seawifs files...\n' )
         
        l3mapgen_prod_list= get_l3mapgen_prod_list(ifile)
        
        print ( '         list of variables found in the mapped.nc file --> ', l3mapgen_prod_list )
        
        if 'chl_oc3'in l3mapgen_prod_list:
            product= 'chl_oc3'       
        if 'chl_oc4'in l3mapgen_prod_list:
            product= 'chl_oc4'
            
        print ( '         png prod name converted from chl_ocx to --> ', product,  '\n' )
                  
    
    if mappping_approach == 'binmap':
        prod_img = asarray(read_hdf_prod(ifile, product))            #read -smi.nc file
        scale_factor,add_offset= get_l3mapgen_slope_intercept(ifile, product)
        prod_img= scale_factor*prod_img + add_offset
        
    if mappping_approach == 'str_map': prod_img = asarray(read_hdf_prod(ifile, product + '-mean'))  #read -map.nc file
         
   
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
    else: print ( 'could not generate png!!!'    )
        
    
 
 
    
#
# check that preliminary conditions are met
# otherwise exit program


def preliminary_checks(length_filelist, length_uniq_syear, length_uniq_sat_type, smi_proj):

    # first check
    if length_filelist == 0:
        print ( '\n########### PROGRAM batch_binmap.pro HALTED BECAUSE NO L2 FILES FOUND IN THE    #########' )
        print ( '########### SPECIFIED L2DIR. PLEASE GO BACK AND RECHECK THE DIRECTORY PATH THAT #########'  )
        print ( '########### WAS SPECIFIED IN THE BATCH_CMD_BINMAP FILE.                         #########\n' )
        sys.exit()
        
    if length_uniq_syear > 1:
        print ( '##### PROGRAM HALTED BECAUSE L2DIR CONTAINS FILES FROM SEPARATE YEARS ####' )
        print ( '##### PLEASE GO BACK AND PARSE L2 FILE INTO SEPARTE UNIQUE YEARS      ####' )
        sys.exit()

    if length_uniq_sat_type > 1:
        print ( '##### PROGRAM HALTED BECAUSE L2DIR CONTAINS FILES FROM SEPARATE SENSORS        ####' )
        print ( '##### PLEASE GO BACK AND PARSE L2 FILE INTO SEPARTE UNIQUE SENSOR DIRECTORIES  ####' )
        
    #if smi_proj != 'RECT' and smi_proj != 'SIN':
    #    print ( '##### PROGRAM HALTED BECAUSE SMI_PROJ WAS NOT SET TO EITHER: CYL OR SIN ####' )
    #    print ( '##### PLEASE GO BACK AND CORRECT THE PROJECTION NAME...                  ####' )
    #    sys.exit()




#
# loop through and call each processing function (l2bin ==> make ascii ==> l3bin ==> smigen ==> png)
#
def process(filelist, time_period, out_dir, products, named_flags_2check, space_res, input_coords, chk, sat_type, year, stats_yesno, smi_proj, mappping_approach):

    print ( '\n\nMAPPING APPROACH ----> ', mappping_approach )


    if time_period == 'DLY':
        ave_dir = 'daily'
    elif time_period == 'WKY':
        ave_dir = 'weekly'
    elif time_period == 'MON':
        ave_dir = 'monthly'


    # make directories
    temp = os.path.dirname(filelist[0])    
    l2bin_output = temp + '/l2bin'
    l3bin_output = temp + '/l3bin'


    if stats_yesno == 'yes'and time_period != 'DLY':
        meas_vec = ['1', '2', '4']       
        meas_names = ['mean', 'variance', 'npoints']
    else: 
        meas_vec = ['1']
        meas_names = ['mean']

    # returns groupings of file names within prescribe time ranges...
    averages = get_average(filelist, time_period, int(year))

    #print ( 'AVERAGES FILE LIST ---->',averages )
    #print ( 'PRODUCTS LIST ---->',products )
    
    # average 
    for file_group in averages:     #file_group = ([start, end], [file1, file2, file3...])
        
        for prod in products:
            
            
            if mappping_approach == 'binmap':  
                
                print ( 'PROCESSING L2 to L3 USING L2bin, L2bin -> SMI + PNG...\n'  )
                         
                # make bl2 files
                bl2_gen(file_group[1], l2bin_output, prod, named_flags_2check, space_res)
        
                # make ascii file
                bl2_files = [os.path.basename(f)[:14] for f in file_group[1]]
                ascii_file_list = ascii_gen(l2bin_output, bl2_files)
   
                # make bl3 files--- one for each file_group
                l3bin_output_dir = l3bin_output + '/' + ave_dir
                l3bin_output_file = l3bin_output_dir + '/' + sat_type + str(year) + file_group[0][0] + str(year) + file_group[0][1] + '.bl3bin'
                bl3_gen(prod, ascii_file_list, input_coords, l3bin_output_file)
        
                # Now I have averaged files
                for s in range(len(meas_vec)):
    
                    # make smigen files
                    smi_output_dir = out_dir + '/' + ave_dir + '/' + prod + '/' + meas_names[s]
                    smi_output_file = smi_output_dir + '/' + os.path.basename(l3bin_output_file)[:15]+'.'+prod+'-'+meas_names[s]+'.smi.nc' 
                    
                    if stats_yesno == 'yes' and time_period != 'DLY':
                        prod_str = prod + ':avg' + ',' +  prod + ':var' ',' + prod + ':nobs' 
                    else:                        
                        prod_str = prod + ':avg'   
                    
                    if os.path.exists(l3bin_output_file):
                        #smi_gen(prod, meas_names[s], meas_vec[s], l3bin_output_file, smi_output_file, smi_proj, space_res, input_coords)
                        l3map_gen(prod_str, meas_names[s], l3bin_output_file, smi_output_file, smi_proj, space_res, input_coords)
                    # make png's for mean smi only...
                    if meas_names[s] == 'mean':
                        if os.path.exists(smi_output_file):
                            png_gen(smi_output_file, smi_output_dir+'/png', prod, meas_names[s], mappping_approach)
            
              
                 
            if mappping_approach == 'str_map':   
                
                print ( 'PROCESSING L2 to L3 USING pyresample.py TO MAP SWATH DATA + PNG...\n'         )
                 
                map_output_dir = out_dir + '/' + ave_dir + '/' + space_res +'m' + '-' + prod + '/' + 'mean'                
                if not os.path.exists(map_output_dir):
                    os.makedirs(map_output_dir)                 
                
                map_basename = sat_type + str(year) + file_group[0][0] + str(year) + file_group[0][1]
                map_output_file = map_output_dir + '/' + map_basename + '.' + time_period + '.' + prod + '.map.nc'                                                         
                
                str_map_gen(file_group[1], map_output_file, prod, smi_proj, input_coords, space_res, named_flags_2check, stats_yesno)
                    
                if os.path.exists(map_output_file):
                    png_gen(map_output_file, map_output_dir+'/png', prod, meas_names[0], mappping_approach)
    
    
      
    if mappping_approach == 'binmap':  
        shutil.rmtree(l2bin_output)
        shutil.rmtree(l3bin_output) 


 
        



#  
# setup environment variables
#
def setup(l2dir, smi_proj, latlon, stats_yesno, color_flags_to_check, sst_flags_to_check):
 
     
     
    # untar any tar files if necessary.... 
    # Large data orders from OBPR come as tared directories that contain a lot of 
    # individually compressed data files.  This function extracts all of compressed 
    # data files from alll of the tared directories found in l1a_dir
    # ---
    if len(glob.glob(l2dir + '/*.tar')) != 0: untar_ocweb(l2dir) 
     
   
     
    # decompress all individual data files if necessary
    # ---
    while any([is_compressed(fi) for fi in glob.glob(l2dir + '/*')]):
        for fi in glob.glob(l2dir + '/*'): decompress_file(fi)      
    
    
    filelist = asarray(glob.glob(l2dir + '/*L2*'))
       
    if len(filelist) == 0:
        print ( 'There are no L2 files in ' + l2dir )
        sys.exit()
        
       
       
    input_coords = map_coords.map_coords()
    if latlon[0] == 'whole':
        input_coords = get_hdf_latlon(filelist[0])
    else:
        input_coords.north = latlon[2]
        input_coords.south = latlon[0]
        input_coords.east = latlon[3]
        input_coords.west = latlon[1] 
    
 
   
    
    print ( '\nL2 to L3 processing started...' )
    print ( 'total number of files found= ', len(filelist), '\n' )
    
    
    root = [os.path.basename(name) for name in filelist]
    sat_type = [name[0] for name in root]
    
   
    
    
    if sat_type[0] != 'L':
        syear = [name[1:5] for name in root]
        #2to3
        #year = asarray(map(int, syear))[0]    # just defined 'year' as the first year in this list...
        year = asarray(list(map(int, syear)))[0]    # just defined 'year' as the first year in this list...
    
    if sat_type[0] == 'L':
        syear = [name[9:13] for name in root]
        #2to3
        #year = asarray(map(int, syear))[0]    # just defined 'year' as the first year in this list...
        year = asarray(list(map(int, syear)))[0]    # just defined 'year' as the first year in this list...
  
    
    # check standrd resolution "sst" or " OC"
    sst_chk =   ['SST' in name for name in root]
    color_chk = [('OC' in name and not 'QKM' in name and not 'HKM' in name and not 'FRS' in name and not 'H' in name[0] and not 'L' in name[0]) for name in root]
       
    # check for hires OC...
    hkm_chk = ['HKM' in name for name in root]  
    qkm_chk = ['QKM' in name for name in root]  
    frs_chk = ['FRS' in name for name in root]    
    hico_chk= ['H' in name[0] for name in root]  
    oli_chk=  ['L' in name[0] for name in root]  
    
    uniq_syear = list(set(sort(syear))) #sorted list of unique years      
    uniq_sat_type = list(set(sort(sat_type))) #sorted list of unique satellites        
  
                                                            
    # make sure things aren't too screwed up  
    preliminary_checks(len(filelist), len(uniq_syear), len(uniq_sat_type), smi_proj)
    
  

    # get l2_flags to check   
    # ---
    if color_flags_to_check == 'standard':
        color_named_flags_2check = get_sds7_default_l2flags(uniq_sat_type[0], 'color')
    else:
        color_named_flags_2check = color_flags_to_check
    
    if sst_flags_to_check == 'standard' and (uniq_sat_type[0]=='A' or uniq_sat_type[0]=='T' or uniq_sat_type[0]=='V'):   
        sst_named_flags_2check = get_sds7_default_l2flags(uniq_sat_type[0], 'sst')
    else:
        sst_named_flags_2check=sst_flags_to_check

    return color_chk, sst_chk, hkm_chk, qkm_chk, frs_chk, hico_chk, oli_chk, color_named_flags_2check, sst_named_flags_2check, sat_type[0], filelist, input_coords, year

    
    
  
      
# -------------------------------------
# process color and/or sst products
# -------------------------------------
#
def batch_proc_L23(l2dir, output_dir='not_specified', products='all', space_res='9', time_period='DLY', color_flags_to_check='standard', 
                    sst_flags_to_check='standard', latlon='whole', smi_proj='RECT', stats_yesno='no', straight_map='no'):

    products = products.split(',') #split products into list (chlor_a,sst,...)
    time_period = time_period.split(',')
    latlon = latlon.split(',')
    
    
    
    # make sure directories are right (/ and ~)
    l2dir = path_reformat(l2dir)
    output_dir = path_reformat(output_dir)
    
    
    
    #setup variables
    color_chk, sst_chk, hkm_chk, qkm_chk, frs_chk, hico_chk, oli_chk, color_named_flags_2check, sst_named_flags_2check, sat_type, filelist, input_coords, year \
        = setup(l2dir, smi_proj, latlon, stats_yesno, color_flags_to_check, sst_flags_to_check)
          

    # put output data next to input data if not specified by user
    if output_dir == 'not_specified':
        output_dir = os.path.dirname(l2dir) + '/' + 'L3_binmap'

    
    # get rid of bad average values, i.e. only keep 'DLY', 'WKY', 'MON'
    def f(x): return (x == 'DLY' or x == 'WKY' or x == 'MON')
    #2to3
    #time_period = filter(f, time_period) #cut empty groups
    time_period = list(filter(f, time_period)) #cut empty groups
    
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    


    #---------------- setup products ---------------------#
    
    # OC (for standard resolution - (to be l2bin-l3bin first before mapping)
    # color_file_indices = where( ((asarray(color_chk)==True) | (asarray(sat_type)=='M')) & (asarray(hires_chk)==hires) )  
    # ---
    color_file_indices = where( asarray(color_chk)==True )                      
    if len(filelist[color_file_indices]) != 0: 
        color_files = list(asarray(filelist)[color_file_indices])    
        color_prod = get_product_list(color_files[0], products)
        
 
    
    # OC HIRES (for hires modis bands or HICO resolution - (to be mapped from L2-swath to map projection using resample.py)
    # ---
    hkm_color_file_indices = where( asarray(hkm_chk)==True )                       
    if len(filelist[hkm_color_file_indices]) != 0: 
        hkm_color_files = list(asarray(filelist)[hkm_color_file_indices])    
        hkm_color_prod = get_product_list(hkm_color_files[0], ['all'])                                                                                                                                                                                                                           
    
    qkm_color_file_indices = where( asarray(qkm_chk)==True )                       
    if len(filelist[qkm_color_file_indices]) != 0: 
        qkm_color_files = list(asarray(filelist)[qkm_color_file_indices])    
        qkm_color_prod = get_product_list(qkm_color_files[0], ['all'])       
          
    frs_color_file_indices = where( asarray(frs_chk)==True )                       
    if len(filelist[frs_color_file_indices]) != 0: 
        frs_color_files = list(asarray(filelist)[frs_color_file_indices])    
        frs_color_prod = get_product_list(frs_color_files[0], products) 
        
    hico_color_file_indices = where( asarray(hico_chk)==True )                       
    if len(filelist[hico_color_file_indices]) != 0: 
        hico_color_files = list(asarray(filelist)[hico_color_file_indices])    
        hico_color_prod = get_product_list(hico_color_files[0], products) 
        
    oli_color_file_indices = where( asarray(oli_chk)==True )                       
    if len(filelist[oli_color_file_indices]) != 0: 
        oli_color_files = list(asarray(filelist)[oli_color_file_indices])    
        oli_color_prod = get_product_list(oli_color_files[0], products) 

          
    #SST   
    sst_file_indices = where( asarray(sst_chk)==True )
        
    if len(filelist[sst_file_indices]) != 0:         
        sst_files = list(asarray(filelist)[sst_file_indices])                      
        sst_prod = ['sst']    #<---<< FIX LATER TO WORK LIKE COLOR PROD ABOVE...      
      
    #-----------------------------------------------------#



    # create mapped output (and optional stats) for each average 
    # period specified (DLY or WKY or MON)
    #--------------------------------------------------------------------------
    for average in time_period:        
        
        
        # map standad resoluton color products...
        # ---        
         if len(filelist[color_file_indices]) != 0: 
              if straight_map == 'no':
                  mappping_approach= 'binmap' 
                  proc_space_res= space_res    #space_res is the resoluton for standard l2 binning set in Batch_Proc.py
                  process(color_files, average, output_dir, color_prod, color_named_flags_2check, \
                        proc_space_res, input_coords,color_chk, sat_type, year, stats_yesno, smi_proj, mappping_approach)
              
              if straight_map == 'yes':
                  mappping_approach= 'str_map'
                  proc_space_res = str(int(resolution_from_sat_fname(color_files[0])))
                  process(color_files, average, output_dir, color_prod, color_named_flags_2check, \
                        proc_space_res, input_coords, color_chk, sat_type, year, stats_yesno, smi_proj, mappping_approach)
                                                
                        
         
          
         # map high resoluton color products...
         # ---          
         if len(filelist[hkm_color_file_indices]) != 0:
             mappping_approach=  'str_map'
             proc_space_res = str(int(resolution_from_sat_fname(hkm_color_files[0])))
             process(hkm_color_files, average, output_dir, hkm_color_prod, color_named_flags_2check, \
                    proc_space_res, input_coords,hkm_chk, sat_type, year, stats_yesno, smi_proj, mappping_approach)                                                                              
                                                                                                                              
         if len(filelist[qkm_color_file_indices]) != 0:
             mappping_approach=  'str_map'
             proc_space_res = str(int(resolution_from_sat_fname(qkm_color_files[0])))   
             process(qkm_color_files, average, output_dir, qkm_color_prod, color_named_flags_2check, \
                    proc_space_res, input_coords,qkm_chk, sat_type, year, stats_yesno, smi_proj, mappping_approach)                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             
 
         if len(filelist[frs_color_file_indices]) != 0:
             mappping_approach=  'str_map'
             proc_space_res = str(int(resolution_from_sat_fname(frs_color_files[0])))     
             process(frs_color_files, average, output_dir, frs_color_prod, color_named_flags_2check, \
                    proc_space_res, input_coords,frs_chk, sat_type, year, stats_yesno, smi_proj, mappping_approach)                                                                                                                                                                                                                                                                                                                                                                                                                                        
        
         if len(filelist[hico_color_file_indices]) != 0:
             mappping_approach=  'str_map'
             proc_space_res = str(int(resolution_from_sat_fname(hico_color_files[0])))    
             process(hico_color_files, average, output_dir, hico_color_prod, color_named_flags_2check, \
                    proc_space_res, input_coords,hico_chk, sat_type, year, stats_yesno, smi_proj, mappping_approach)                 
         
         if len(filelist[oli_color_file_indices]) != 0:
             mappping_approach=  'str_map'
             proc_space_res = str(int(resolution_from_sat_fname(oli_color_files[0])))    
             process(oli_color_files, average, output_dir, oli_color_prod, color_named_flags_2check, \
                    proc_space_res, input_coords,oli_chk, sat_type, year, stats_yesno, smi_proj, mappping_approach)     
         
        
        
        # map standard resolution sst
        # ---
         if len(filelist[sst_file_indices]) != 0:
             if straight_map == 'no':
                 mappping_approach= 'binmap' 
                 proc_space_res= space_res    #space_res is the resoluton for standard l2 binning set in Batch_Proc.py
                 process(sst_files, average, output_dir, sst_prod, sst_named_flags_2check, \
                        proc_space_res, input_coords, sst_chk, sat_type, year, stats_yesno, smi_proj, mappping_approach)
                                                                                
             if straight_map == 'yes':
                 mappping_approach= 'str_map'
                 proc_space_res = str(int(resolution_from_sat_fname(sst_files[0])))
                 process(sst_files, average, output_dir, sst_prod, sst_named_flags_2check, \
                        proc_space_res, input_coords, sst_chk, sat_type, year, stats_yesno, smi_proj, mappping_approach)

#
# defaults
#   
def main(*args):
    import getopt
    
    arg_options = ['l2_dir=', 'binmap_dir=','products=', 'space_res=', 'time_period=','color_flags=','sst_flags=','latlon=','smi_proj=','stats_yesno=','straight_map']
    opts, arg = getopt.getopt(args, '', arg_options)

    if len(args) == 0:
        print ( '\nUsage:\n\t batch_L23.py --l2_dir=<l2_dir> --binmap_dir=<l2_dir> --products=<prod1,prod2,prod3,...> --space_res=<space_res> --time_period=<time_period> --color_flags=<flag1,flag2,...> --sst_flags=<flag1,flag2,...> --latlon=<S,W,N,E> --smi_proj=<smi_proj> --stats_yesno=<stats_yesno> --straight_map=<straight_map>\n\n')
        print ( '\t--'+arg_options[0][:-1]+' (required) ==> directory containing Level 2 data files' )
        print ( '\t--'+arg_options[1][:-1]+' (optional, default=L3_binmap) ==> new directory that the l2bin and l3bin data files will be written to' )
        print ( '\t--'+arg_options[2][:-1]+' (optional, defualt=all) ==> products contained in the L2 files to be mapped (chlor_a,poc,cdom_index,...)' )
        print ( '\t--'+arg_options[3][:-1]+' (optional, default=9) ==> Spatial Resolution in km (1,4,9,36)' )
        print ( '\t--'+arg_options[4][:-1]+' (optional, defualt=DLY) ==> Time period to be averages (DLY, WKY, MON)' )
        print ( '\t--'+arg_options[5][:-1]+' (optional, defualt=standard) ==> Color flags to check' )
        print ( '\t--'+arg_options[6][:-1]+' (optional, default=standard) ==> SST flags to check' )
        print ( '\t--'+arg_options[7][:-1]+' (optional, default=whole) ==> Lat/LON --S,W,N,E. If not specified, finds latlon from HDF file' )
        print ( '\t--'+arg_options[8][:-1]+' (optional, defualt=RECT) ==> SMI projection (RECT = rectangular, SIN = sinusoidal)' )
        print ( '\t--'+arg_options[9][:-1]+' (optional, defualt=no) ==> Stats Yes/No' )
        print ( '\t--'+arg_options[10][:-1]+' (optional, defualt=no) ==> Straight_Map Yes/No' )
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
            if option == '--' + arg_options[7][:-1]:
                arg8 = value
            if option == '--' + arg_options[8][:-1]:
                arg9 = value
            if option == '--' + arg_options[9][:-1]:
                arg10 = value
            if option == '--' + arg_options[10][:-1]:
                arg11 = value

        if 'arg2' not in locals():
            arg2 = 'not_specified'
        if 'arg3' not in locals():
            arg3 = 'all'
        if 'arg4' not in locals():
            arg4 = '9'
        if 'arg5' not in locals():
            arg5 = 'DLY'
        if 'arg6' not in locals():
            arg6 = 'standard'
        if 'arg7' not in locals():
            arg7 = 'standard'
        if 'arg8' not in locals():
            arg8 = 'whole'
        if 'arg9' not in locals():
            arg9 = 'RECT'
        if 'arg10' not in locals():
            arg10 = 'no'
        if 'arg11' not in locals():
            arg11 = 'no'

        batch_proc_L23(arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10,arg11)
        
        
#
# COMMAND LINE CALL
#
if __name__=='__main__':
    
    main(*sys.argv[1:])
    
            
