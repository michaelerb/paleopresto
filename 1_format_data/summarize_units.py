#==============================================================================
# Print some basic information about each file.
#    author: Michael P. Erb
#    date  : 5/4/2023
#==============================================================================

import numpy as np
import xarray as xr
import glob

#%% LOAD DATA
data_dir = '/projects/pd_lab/data/paleoclimate_reconstructions/presto_format/'
filenames_all = glob.glob(data_dir+'*.nc')

#%% Loop through the files, printing some basic information
for filename in filenames_all:
    #
    filename_basic = filename.split('/')[-1]
    var_txt = filename_basic.split('_')[-2]
    #
    # Load the data
    handle = xr.open_dataset(filename)
    var_units = handle[var_txt+'_spatial_mean'].units
    handle.close()
    #
    # Print some information
    print(' === '+filename_basic+' | ',var_txt,' ===')
    print('Units:',var_units)

