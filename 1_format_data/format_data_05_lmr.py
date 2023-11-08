#==============================================================================
# Make a standardized netcdf file for the LMR.
#    author: Michael P. Erb
#    date  : 11/7/2023
#==============================================================================

import sys
import numpy as np
import xarray as xr
import utils


#%% SETTINGS

#version_txt = '2.0'
#version_txt = '2.1'
#var_txt     = 'tas'
#var_txt     = 'precip'
version_txt = sys.argv[1]; var_txt = sys.argv[2]


#%% LOAD DATA
# Note: data can be downloaded at: https://www.ncei.noaa.gov/access/paleo-search/study/27850

print('=== Processing LMR, v'+version_txt+' ===')
data_dir = '/projects/pd_lab/data/paleoclimate_reconstructions/LMR_final/'

# Load data
if var_txt == 'tas':
    #
    data_xarray = xr.open_dataset(data_dir+'air_MCruns_ensemble_mean_LMRv'+version_txt+'.nc')
    var_spatial_mean = data_xarray['air'].values
    time             = data_xarray['time'].values
    lat              = data_xarray['lat'].values
    lon              = data_xarray['lon'].values
    data_xarray.close()
    #
    data_xarray = xr.open_dataset(data_dir+'gmt_MCruns_ensemble_full_LMRv'+version_txt+'.nc')
    var_global_members = data_xarray['gmt'].values
    data_xarray.close()
    #
    # Format the global variables, moving all ensemble members to the ensemble member dimension
    n_time = var_global_members.shape[0]
    n_iter = var_global_members.shape[1]
    n_ens  = var_global_members.shape[2]
    var_global_members = np.reshape(np.swapaxes(var_global_members,0,2),(1,n_ens*n_iter,n_time))
    var_global_mean    = np.mean(var_global_members,axis=1)
    #
    ens_spatial = ['mean']
    ens_global  = np.arange(var_global_members.shape[1])+1
    var_units = 'degrees Celsius'
    #
elif var_txt == 'precip':
    #
    data_xarray = xr.open_dataset(data_dir+'prate_MCruns_ensemble_mean_LMRv'+version_txt+'.nc')
    var_spatial_mean = data_xarray['prate'].values
    time             = data_xarray['time'].values
    lat              = data_xarray['lat'].values
    lon              = data_xarray['lon'].values
    lat_weights = np.cos(np.deg2rad(data_xarray.lat))
    var_global_mean = np.array(data_xarray.prate.weighted(lat_weights).mean(('lon','lat')))  # Compute the global mean
    data_xarray.close()
    #
    # Convert precipitation units from kg/m2/s to mm/day
    var_spatial_mean = var_spatial_mean*60*60*24
    var_global_mean  = var_global_mean*60*60*24
    #
    # Format the global variables
    var_global_mean    = np.mean(var_global_mean,axis=1)[None,:]
    var_global_members = np.expand_dims(var_global_mean,axis=1)
    #
    ens_spatial = ['mean']
    ens_global  = ['mean']
    var_units = 'mm per day'

# Get the years and ages
year = np.array([time_selected.year for time_selected in time])
age = 1950-year

# Compute the mean of the ensemble members
var_spatial_mean = np.mean(var_spatial_mean,axis=1)[:,None,:,:]


#%% CALCULATIONS

# Calculate lat and lon bounds
lat_bounds,lon_bounds = utils.bounding_latlon(lat,lon)

# Format the spatial variables
var_spatial_mean    = np.swapaxes(var_spatial_mean,0,1)
var_spatial_members = np.expand_dims(var_spatial_mean,axis=1)

# Get other metadata
methods = ['LMR']

# If this data can't be reformatted to the standard format, add a note here 
notes = ['']

# Check the shape of the variables
print(var_spatial_members.shape)
print(var_spatial_mean.shape)
print(var_global_members.shape)
print(var_global_mean.shape)


#%% FORMAT DATA

# Create new array
data_xarray_output = xr.Dataset(
    {
        var_txt+'_global_mean':    (['method','age'],                          var_global_mean,    {'units':var_units}),
        var_txt+'_global_members': (['method','ens_global','age'],             var_global_members, {'units':var_units}),
        var_txt+'_spatial_mean':   (['method','age','lat','lon'],              var_spatial_mean,   {'units':var_units}),
        var_txt+'_spatial_members':(['method','ens_spatial','age','lat','lon'],var_spatial_members,{'units':var_units})
    },
    coords={
        'method':     (['method'],methods),
        'notes':      (['notes'],notes),
        'ens_global': (['ens_global'],ens_global,{'description':'ensemble members'}),
        'ens_spatial':(['ens_spatial'],ens_spatial,{'description':'means of the 20 Monte Carlo interations'}),
        'age':        (['age'],age,{'units':'yr BP'}),
        'lat':        (['lat'],lat,{'units':'degrees_north'}),
        'lon':        (['lon'],lon,{'units':'degrees_east'}),
        'lat_bounds': (['lat_bounds'],lat_bounds,{'units':'degrees_north'}),
        'lon_bounds': (['lon_bounds'],lon_bounds,{'units':'degrees_east'}),
    },
    attrs={
        'dataset_name':      'LMR',
        'dataset_source_url':'https://www.ncei.noaa.gov/access/paleo-search/study/27850',
    },
)


#%% SAVE DATA

# Save new array
output_dir = '/projects/pd_lab/data/paleoclimate_reconstructions/presto_format/'
data_xarray_output.to_netcdf(output_dir+'lmr_v'+version_txt.replace('.','_')+'_0_'+var_txt+'_annual.nc')

