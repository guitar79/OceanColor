#! /usr/bin/env python

from my_mapping_utilities import png_rescale

# This tiny little script can be used to rescale the min and max values for
# the png image output of a given mapped netCDF4 data product file that was 
# processed using the batch_L23 script...


# FIRST, You can use SeaDAS Interactive (GUI) mode on one of the netcdf4 (mapped or smi) files
#        and carefully apply different scales min/max and log/lin settings to find the best settings 
#        for your batch products.
# 
# SECOND, you can go to the following text file to set up your new min/max scale settings that 
#         are specific for your product OR for your specific geographic location or season:
# 
# ===> /Users/bmonger/python_programs/local_procesing_resources/png_min_max_settings/prod_min_max_tab_delimted_txt
# 
# Note: To add a new product min/max to the current list of products to scale, append a line to the 
#       prod_min_max_scales_tbl.txt file  that looks identical to the line above, then then change 
#       the product name to the scale info a appropriate for your product.
# 
# 
#
# THIRD, 1) prescribe (see below) the directory path to the netCDF4 mapped data product files...
#        2) prescribe (see below) the name of the data product to be rescaled for new png output...
#
# -----------------------------------------------------------------------------------------------------------------------------------------




# Directroy path to the mapped netCDF4 data product files (.nc files)
# ---
level3_nc_file_dir= '/rsclass/seadas_examples/test_aqua_13jan13/daily/chlor_a/mean'


# Name of the data product being rescaled for png image output...
# ---
data_product= 'chlor_a'




png_rescale(level3_nc_file_dir,data_product)
