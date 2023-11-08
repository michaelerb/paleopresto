#==============================================================================
# Make a standardized netCDF file for Holocene Reconstruction.
#    author: Michael P. Erb
#    date  : 9/12/2023
#==============================================================================

import numpy as np
import xarray as xr
import utils


#%% LOAD DATA
# Note: data can be downloaded at: https://zenodo.org/record/6426332

print('=== Processing Holocene Reconstruction ===')

# Load data
data_dir = '/projects/pd_lab/data/data_assimilation/final_results/'
data_xarray = xr.open_dataset(data_dir+'holocene_reconstruction.nc')
var_global_members  = data_xarray['recon_tas_global_mean'].values
var_spatial_mean    = data_xarray['recon_tas_mean'].values
var_spatial_members = data_xarray['recon_tas_ens'].values
age = data_xarray['ages'].values
lat = data_xarray['lat'].values
lon = data_xarray['lon'].values


#%% CALCULATIONS

# Calculate lat and lon bounds
lat_bounds,lon_bounds = utils.bounding_latlon(lat,lon)

# Format the variables
var_spatial_mean    = np.expand_dims(var_spatial_mean,axis=0)
var_spatial_members = np.expand_dims(np.swapaxes(var_spatial_members,0,1),axis=0)
var_global_members  = np.expand_dims(np.swapaxes(var_global_members,0,1),axis=0)
var_global_mean     = np.mean(var_global_members,axis=1)

# Get other metadata
methods = ['Holocene Reconstruction']
ens_spatial = np.arange(var_spatial_members.shape[1])+1
ens_global  = np.arange(var_global_members.shape[1])+1

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
        'tas_global_mean':    (['method','age'],                          var_global_mean,    {'units':'degrees Celsius'}),
        'tas_global_members': (['method','ens_global','age'],             var_global_members, {'units':'degrees Celsius'}),
        'tas_spatial_mean':   (['method','age','lat','lon'],              var_spatial_mean,   {'units':'degrees Celsius'}),
        'tas_spatial_members':(['method','ens_spatial','age','lat','lon'],var_spatial_members,{'units':'degrees Celsius'})
    },
    coords={
        'method':     (['method'],methods),
        'notes':      (['notes'],notes),
        'ens_global': (['ens_global'],ens_global,{'description':'ensemble members'}),  #TODO: The global and spatial ensemble members won't match. Look into this in all scripts.
        'ens_spatial':(['ens_spatial'],ens_spatial,{'description':'selected ensemble members'}),
        'age':        (['age'],age,{'units':'yr BP'}),
        'lat':        (['lat'],lat,{'units':'degrees_north'}),
        'lon':        (['lon'],lon,{'units':'degrees_east'}),
        'lat_bounds': (['lat_bounds'],lat_bounds,{'units':'degrees_north'}),
        'lon_bounds': (['lon_bounds'],lon_bounds,{'units':'degrees_east'}),
    },
    attrs={
        'dataset_name':      'Holocene Reconstruction',
        'dataset_source_url':'https://zenodo.org/record/6426332',
    },
)


#%% SAVE DATA

# Save new array
output_dir = '/projects/pd_lab/data/paleoclimate_reconstructions/presto_format/'
data_xarray_output.to_netcdf(output_dir+'daholocene_v1_0_0_tas_annual.nc')

