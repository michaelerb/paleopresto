#==============================================================================
# Make a standardized netcdf file for the LMR.
#    author: Michael P. Erb
#    date  : 5/5/2023
#==============================================================================

import sys
import numpy as np
import xarray as xr
import utils
import matplotlib.pyplot as plt


#%% SETTINGS

version_txt = '2.0'
var_txt     = 'tas'

print('=== Processing LMR, v'+version_txt+' ===')
data_dir = '/projects/pd_lab/data/paleoclimate_reconstructions/LMR_final/'

# Load data
data_xarray = xr.open_dataset(data_dir+'air_MCruns_ensemble_mean_LMRv'+version_txt+'.nc')
var_ens = data_xarray['air']
time    = data_xarray['time'].values
lat     = data_xarray['lat'].values
lon     = data_xarray['lon'].values

data_xarray_stdev = xr.open_dataset(data_dir+'air_MCruns_ensemble_spread_LMRv'+version_txt+'.nc')
var_ens_stdev = data_xarray_stdev['air']
data_xarray_stdev.close()

data_xarray_gmt = xr.open_dataset(data_dir+'gmt_MCruns_ensemble_full_LMRv'+version_txt+'.nc')
var_global = data_xarray_gmt['gmt'].values
data_xarray_gmt.close()
var_global = np.reshape(var_global,(2001,20*100))

# Get the years and ages
year = np.array([time_selected.year for time_selected in time])
age = 1950-year


#%% CALCULATIONS

# Process data
var_mean = np.mean(var_ens,axis=1)

# Since the ensemble members aren't available, approximate the 2*standard deviation values
var_upper_2std = np.mean((var_ens + (var_ens_stdev)),axis=1)  # Mean value for the lower 2*std bound
var_lower_2std = np.mean((var_ens - (var_ens_stdev)),axis=1)  # Mean value for the upper 2*std bound


#%%
# Compute the global mean
lat_weights  = np.cos(np.deg2rad(data_xarray.lat))
#var_global = data_xarray.air.weighted(lat_weights).mean(('lon','lat'))
var_upper_2std_global = var_upper_2std.weighted(lat_weights).mean(('lon','lat'))
var_lower_2std_global = var_lower_2std.weighted(lat_weights).mean(('lon','lat'))


#%%

plt.plot(age,var_upper_2std_global)
plt.plot(age,var_lower_2std_global)
#plt.plot(age,np.percentile(var_global,97.5,axis=1))
#plt.plot(age,np.percentile(var_global,2.5,axis=1))
