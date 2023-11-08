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
    #var_txt = filename_basic.split('_')[-2]
    #
    # Load the data
    handle = xr.open_dataset(filename)
    ens_spatial = handle['ens_spatial'].values
    ens_global  = handle['ens_global'].values
    #var_mean   = handle[var_txt+'_mean']
    #var_ens    = handle[var_txt+'_ens']
    #age        = handle['age'].values
    #age_units  = handle['age'].units
    handle.close()
    #
    # Remove the last 1000 years from the reconstruction
    #ind_ref = np.where((age >= 0) & (age <= 1000))[0]
    #var_ens  = var_ens  - np.nanmean(var_mean[:,ind_ref,:,:],axis=1)
    #var_mean = var_mean - np.nanmean(var_mean[:,ind_ref,:,:],axis=1)
    #
    # Print some information
    print(' === '+filename_basic+' ===')
    #print('Range of ages ('+age_units+'):',age[0],age[-1])
    #print('Range of var mean ('+var_units+'):',var_mean.min().values,var_mean.max().values)
    #print('Range of var ens ('+var_units+'): ',var_ens.min().values, var_ens.max().values)
    print('Spatial:',ens_spatial)
    print('Global:',ens_global)

