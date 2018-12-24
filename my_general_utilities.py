#!/usr/bin/env python
# -*- coding: utf-8 -*-
        
import os, sys
sys.dont_write_bytecode = True
from subprocess import call
import time

import numpy as np
import glob
import zipfile, tarfile
import shutil

from my_hdf_cdf_utilities import *
import map_coords

import subprocess #Added for Sean Bailey's fix in lines 172... to avoid reading hdf file
import bz2


from numpy import *


def rebin(a,new_shape):    
    a_view= a.reshape(new_shape[0],a.shape[0]/new_shape[0],new_shape[1],a.shape[1]/new_shape[1])
    return a_view.mean(axis=3).mean(axis=1)
   

 
def filebreak(full_fname):
    
    if isinstance(full_fname,basestring) == 1:
       full_fname= [full_fname]
  
    result= np.empty([len(full_fname),2], dtype=object)
    
    for i in range(len(full_fname)):
        result[i,0]=os.path.dirname(full_fname[i])
        result[i,1]=os.path.basename(full_fname[i])    
    
    return result
  

def untar_ocweb(indir):
    
    tar_fnames= glob.glob(indir + '/*.tar')
    for name in tar_fnames:
        os.system('tar -C ' + indir + ' -xvf ' + name)                    
    os.system('mv ' + indir + '/requested_files/* ' + indir) 
    os.system('rmdir ' + indir + '/requested_files')
    os.system('rm ' + indir + '/*.tar')                                                                                                                                                                                                                                                                                                                                 
   
             

# input: file ==> data file with name including 'satID,year,jday,hour,min,sec', e.g. 'A2013125220000'
# output: integer julian day
#
def get_jday(file):
    file_base = os.path.basename(file)
    
    if file_base[0] != 'L': return int(file_base[5:8])
    if file_base[0] == 'L': return int(file_base[13:16])
      


  
# unzips and extracts tar and zip to same directory
# puts old compressed files in separate directory
# input:
#   file ==> zipped file
#
def decompress_file(file):
    file_base = os.path.basename(file)
    path = os.path.dirname(file)
    new_dir = os.path.dirname(path) + '/compressed/'
    
    if not os.path.exists(new_dir):
        os.makedirs(new_dir)

    if zipfile.is_zipfile(file):
        zip = zipfile.ZipFile(file)
        zip.extractall(path + '/')
        shutil.move(file, new_dir + file_base)
    elif file[-3:] == 'bz2':
        print ( 'extracting bz2' )
        print ( file )
        os.system('bunzip2 -k ' + file)
        shutil.move(file, new_dir + file_base)
    elif tarfile.is_tarfile(file):
        print ( 'extracting tar' )
        print ( file )
        print ( path )
        os.system('tar -C ' + path + ' -xvf ' + file)
        print ( 'check 1' )
        print ( glob.glob(path + '/*') )
        shutil.move(file, new_dir + file_base)
        print ( 'check 2' )
        print ( glob.glob(path + '/*') )
    elif (file.find('.gz') != -1):
        outfile = file[:file.rindex('.')]
        os.system("gunzip -c " + file + ' > ' + outfile)
        shutil.move(file, new_dir + file_base)
    else:
        print ( "Wrong archive or filename" )
    


#
# returns True is file is a zipped or tarred file
# else returns false
#        
def is_compressed(file):

    if os.path.isfile(file) ==True:
        return zipfile.is_zipfile(file) or tarfile.is_tarfile(file) or (file.find('.gz') != -1) or (file.find('bz2') != -1)
    else: return ('')



# inputs: 
#   data ==> 2D array of values
#   latitudes ==> 1D array of latitudes
#   longitudes ==> 1D array of longitudes
#   xdim ==> number of x pixels in output image
#   ydim ==> number of y pixels in output image
# output: 2D array of size (xdim, ydim)
#
def map_resize(data, latitudes, longitudes, xdim, ydim):
    lon1 = -180.0
    lon2 = 180.0
    lat1 = -90.0
    lat2 = 90.0
    
    mapped = np.empty((ydim,xdim))
    mapped[:,:] = np.nan

    pixels_per_lat = ydim/(lat2 - lat1)
    pixels_per_lon = xdim/(lon2 - lon1)
    
    # map longitudes and latitudes to pixel number
    lon_pix = np.where(longitudes <= 180, (longitudes+180)*pixels_per_lon, (longitudes-179)*pixels_per_lon)
    lat_pix = map( lambda lat: (lat + 90)*pixels_per_lat, latitudes )

    # put each lat/lon in its spot on rectangular grid
    def func(lat,lon,value): mapped[lat,lon] = value    
    map( lambda i: map(lambda j: func(lat_pix[i], lon_pix[j], data[i,j]) ,range(0,len(lon_pix))) ,range(0,len(lat_pix)) )
    
    
    return mapped

    

#
# get rid of ~ at beginning and / at end of file name
#
def path_reformat(path):
    path = path.strip()

    # get rid of troublesome '/' at end
    if path[-1] == '/':
        path = path[:-1]

    # get complete path
    if path[0] == '~':
        path = os.path.expanduser(path)
        
    return path
  
  
  
#
# input: list of files [file1,file2,...], list of averages ['DLY','WKY','MON'] and an integer year
# output: [([start, end], [file1,file2,...]), ([start, end], [file1,file2,...]), ([start, end] ,[file1,file2,...]), ...]
#
def get_average(filelist, time_period, year):

    if time_period == 'DLY':
        start_day = np.arange(1,366)
        end_day = np.arange(1,366)
    elif time_period == 'WKY':
        start_day = np.arange(1,366,8)
        end_day = np.asarray(map(lambda i: i+7, start_day))
    elif time_period == 'MON':
        # check for leap years
        if (year - 1980)%4 != 0:
            start_day = np.array([1,32,60,91,121,152,182,213,244,274,305,335])
            end_day = np.array([31,59,90,120,151,181,212,243,273,304,334,365])
        else:
            start_day = np.array([1,32,61,92,122,153,183,214,245,275,306,336])
            end_day = np.array([31,60,91,121,152,182,213,244,274,305,335,366])
 
    grouping_list = [([],[]) for i in range(0,len(start_day))] #list of tuples
    
    for i in range(0,len(start_day)):
        grouping_list[i][0].append('%03d' % start_day[i])
        grouping_list[i][0].append('%03d' % end_day[i])
        
    # assign each file to a time group
    for file in filelist:

        for i in range(0,len(start_day)):

            if get_jday(file) >= start_day[i] and get_jday(file) <= end_day[i]:
                grouping_list[i][1].append(file)

    def f(x): return len(x[1])!=0
    grouping_list = filter(f, grouping_list) #cut empty groups
    
    
    return grouping_list
    
           



def bulk_extract_meris(l1a_dir, extract_dir, extract_latlon):
#----------------------------------------------------------------------- 

    if not os.path.exists(extract_dir): os.makedirs(extract_dir)
                        
    latlon = extract_latlon.split(',')                   
    south=   latlon[0]  
    west=    latlon[1]
    north=   latlon[2] 
    east=    latlon[3] 
         
    for fi in glob.glob(l1a_dir + '/' + '*FRS*N1*'): decompress_file(fi)    

    fname_l1a = glob.glob(l1a_dir + '/*FRS*N1*')            #re-query now without .gz in fname            
    root_names = [os.path.basename(i) for i in fname_l1a]   #list of file names split from path
  
            
    for i in range(len(root_names)):
        
        root_pieces= root_names[i].split('_')
        
        resoluton= root_pieces[1]
        yyyy=        root_pieces[2][6:10]        
        mm=          root_pieces[2][10:12]
        dd=          root_pieces[2][12:14]
        hhmmss=      root_pieces[3]  
       
                                                                                               
        result= time.strptime(dd+' '+ mm+' '+yyyy, '%d %m %Y')
        jday= result[7]   # the function time.structure returns a full time structure (tuple) 
                          # of which the 7th elsement is the julian day === year_day
    
         
        if jday <= 9: s_jday=  '00' + str(jday)
        if jday >= 10 and jday <= 99: s_jday=  '0' + str(jday)
        if jday >= 100: s_jday=  str(jday)
                                                
        ofname= extract_dir + '/' + 'M' + yyyy + s_jday + hhmmss + '.L1B_' + resoluton                             
                     
        #find start and stop scan lines to coorsponding to prescribed lat/lon bounds        
        #then extract the full meris N1 scene to the smaller extarced scene.
        #note new extraced file has the classic obpg name convention
        #---                             
        p= subprocess.Popen('lonlat2pixline -v ' + fname_l1a[i] + ' ' + west  + ' ' + south  + ' ' + east + ' ' + north, shell=True,stdout=subprocess.PIPE, stderr=subprocess.PIPE)    
        stdout,stderr = p.communicate()
        
       
        for line in stdout.split('\n'):
            if 'sline' in line:
                junk = line.split('=')
                sline= junk[1]
            if 'eline' in line:
                junk = line.split('=')           
                eline= junk[1]    
                            
        print ( '\nextracting meris file ---> ',  fname_l1a[i]  )
        print ( 'new extraced file to be generated:-----> ', ofname )
        print ( '\n'  )
                                                                
        subprocess.Popen( 'l1bextract_meris ' + fname_l1a[i]+' '+ sline +' '+ eline +' '+ ofname,shell=True,stdout=subprocess.PIPE, stderr=subprocess.PIPE) 
       

 
def wget_oceans_gsfc(fname):
#-----------------------------------------------------------------------

    # USAGE:  Termina Prompt>  wget_oceans_gsfc('my_list_of_files_to_get.txt') 
    
    # where my_list_of_files_to_get.txt  is a standard text file containing the list of all
    # files (including "http://...") that you want to retreive.  To get the list of files
    # see below.

    #when at your home institution, first check the wget is installed on your computer by typing
    #the following at the unix prompt: which  wget
    #if wget is not installed then go on the web a search for wget for your type of machine (mac vs linux)
    #Install wget in /usr/local/bin and double check that /usr/local/bin in the the PATH by typing 
    #echo $PATH

    #Go to the oceancolor website's file serach page and specify the files you want: 
        #>>>----->  http://oceandata.sci.gsfc.nasa.gov/cgi/file_search.cgi

        #specify the search criteria for the files you want (no leading blank space before or 
        #after specificed search term).

        # NOTE: >>------> select all three options at the bottom of the page (i.e., check Downlaod result and Add URL)

    #save the file list as a text file (or copy and past from web page list into a new text doc) to a new 
    #directtory where you want to save the downloaded files. Then cd to that directory and:
        # 1. launch python
        # 2. type:  from my_general_utilities import wget_oceans_gsfc
        # 3. Type:  wget_oceans_gsfc('your_fname.txt')
 
    
    f = open(fname)
    all_lines = f.readlines()
    f.close()
    
    for url_name in all_lines:
        if url_name[0:4] == 'http':  
            call(['wget', url_name.strip('\n')])
    
    
 
def wget_aviso(username, password, filelist):
#----------------------------------------------------------------------- 
 
    # NOTE:---> When you get your own user account, replace below 
    #           cornell2:gbp47k3 with your own username:password
 

    # NOTE---> cd to the directory where the text containing the list of
    #          files (which should also be the directory where you want the
    #          new files to downlaod).  
    #          # 1. type: python  
    #          # 2. type: from my_general_utilities import wget_aviso
               # 3. Type:  wget_aviso('username', 'password','your_fname.txt')
        
 
    f = open(filelist)
    all_lines = f.readlines()
    f.close()
   
     
    url_path=  'ftp://' + username + ':' + password + '@ftp.aviso.altimetry.fr/global/delayed-time/grids/madt/all-sat-merged' 
   
    uv_chk = np.asarray(['_uv_' in name for name in all_lines])
    if uv_chk.all() == True: url_path=  url_path + '/uv'  
    
    h_chk=   np.asarray(['_h_'  in name for name in all_lines])
    if h_chk.all() == True: url_path =  url_path + '/h'
    
    if uv_chk.all() == True or h_chk.all() == True:
    
        for aviso_fname in all_lines:
            
            aviso_fname= aviso_fname.strip('\n')
            
            if aviso_fname[0:4] == 'File':                #some recs have the form:   "File:dt_global_allsat_madt_uv_20120107_20140106.nc.gz  2483 KB  3/3/14  3:31:00 AM" 
                split_result= aviso_fname.split(':')      #this piece of code remotes the: "File:"
                aviso_fname=  split_result[1] 
                 
            
            split_result = aviso_fname.split(' ')         #some recs have the form:   "dt_global_allsat_madt_uv_20120107_20140106.nc.gz  2483 KB  3/3/14  3:31:00 AM"  
            aviso_fname= split_result[0]                  #this piece of code removes text trailing info to leave only "dt_global_allsat_madt_uv_20120107_20140106.nc.gz"
                
            split_result= aviso_fname.split('_')
      
            yyyymmdd= split_result[5]
            s_year= yyyymmdd[0:4]
            url_path_year= url_path + '/' + s_year
        
            call(['wget', url_path_year + '/' +  aviso_fname])
   
    
   
 
def wget_ghrsst_l3p_sst(year, start_day, end_day):
#-----------------------------------------------------------------------

    #NOTE:  FOR THIS PROGRAM TO WORK YOU MUST HAVE A TEXT FILE CALLED: .netrc LOCATED IN YOUR HOME DIRECTORY
    #       WHERE WHERE THE .netrc FILE CONTAINS THE FOLLOWING LINE...
    # 
    #      
    #       machine podaac-ftp.jpl.nasa.gov login anonymous password your_email_address 
    # 
    #        
    #  TO USE THIS FTP FUNCTION, YOU SHOULD CD TO THE DATA DIRECTORY YOU WANT THE DATA TO BE
    #  TRANSFERRED TO, THEN LAUNCH PYTHON AND RUN THIS PROGRAM FROM THE PYTHON PROMPT USING THE 
    #  FOLLOWING SYNTAX EXAMPLE
    #  
    #  from my_general_utilities import wget_ghrsst_l3p_sst
    #  wget_ghrsst_l3p_sst(year=2013, start_day=1, end_day=366)   # you can choose any day interval...
    #

    year= str(year)

    jday_vec=      start_day + np.arange(end_day-start_day+1)
    sjday_vec=     str(jday_vec)

    sjday_vec = map(str, jday_vec)
         
    for i in range(len(sjday_vec)):
        if len(sjday_vec[i]) == 1: sjday_vec[i] = '00' + sjday_vec[i]
        if len(sjday_vec[i]) == 2: sjday_vec[i] =  '0' +  sjday_vec[i]
    
    search_string1= '*L3P*120000-v01.7-fv01.0.nc.bz2'
    search_string2= '*L3P*000000-v01.7-fv01.0.nc.bz2'
 
    for i in range(len(sjday_vec)):
        call(['wget', 'ftp://podaac-ftp.jpl.nasa.gov/allData/ghrsst/data/L3P/GLOB/AVHRR_METOP_A/EUR' + '/' + year + '/' + sjday_vec[i] + '/' + search_string1])
        call(['wget', 'ftp://podaac-ftp.jpl.nasa.gov/allData/ghrsst/data/L3P/GLOB/AVHRR_METOP_A/EUR' + '/' + year + '/' + sjday_vec[i] + '/' + search_string2])
   

 


 
def rebin_down_nan(orig_img, new_ydim, new_xdim):
#-----------------------------------------------------------------------

    #THIS FUNCTION ACTS JUST LIKE IDL'S REBIN FUNCTION WHEN REDUCING AN ARRAY BUT WITH THE
    #IMPORTANT EXCEPTION OF ALLOWING NAN IN THE ORIGINAL ARRAY TO BE TREATED AS MISSING
    #VALUES WHEN AVERAGING OVER PIXEL BLOCKS IN THE ORGINAL ARRAY

    #WRITTEN BY BRUCE MONGER, CORNELL UNIVERSITY, FEBRUARY 14, 2015

    orig_ydim, orig_xdim = orig_img.shape
 

    if (np.mod(orig_xdim, new_xdim) != 0) or (np.mod(orig_ydim, new_ydim) != 0): 
        print ( 'dimensions of new array must be integer value of orginal array dimmension' )
        return ['']
     

    new_img=   np.zeros((new_ydim,new_xdim), dtype=float)
    xbox=      orig_xdim/new_xdim
    ybox=      orig_ydim/new_ydim


    start_xi=  np.arange(new_xdim)*xbox
    end_xi=    start_xi + (xbox-1)
    start_yi=  np.arange(new_ydim)*ybox
    end_yi=    start_yi + (ybox-1)
 
 
    for i in range(new_ydim):
        for j in range(new_xdim): 
            box=   orig_img[start_yi[i]:end_yi[i],start_xi[j]:end_xi[j]]
            new_img[i,j]=  np.sum(np.nan_to_num(box))/np.sum(~np.isnan(box).astype(int))
     

    return new_img




