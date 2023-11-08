#==============================================================================
# Make a standardized netCDF file for OWDA.
#    author: Michael P. Erb
#    date  : 9/12/2023
#==============================================================================

import numpy as np
import xarray as xr
import utils


#%% LOAD DATA
# Note: data can be downloaded at: https://www.ncei.noaa.gov/access/paleo-search/study/19419

print('=== Processing OWDA ===')
data_dir = '/projects/pd_lab/data/paleoclimate_reconstructions/drought_atlases/'

# Load data
data_xarray = xr.open_dataset(data_dir+'owda.nc')
var_spatial_mean = data_xarray['pdsi'].values
years            = data_xarray['time'].values
lat              = data_xarray['lat'].values
lon              = data_xarray['lon'].values
data_xarray.close()
age = 1950-years


#%% CALCULATIONS

# Compute the spatial mean
lat_weights = np.cos(np.deg2rad(data_xarray.lat))
var_global_mean = data_xarray['pdsi'].weighted(lat_weights).mean(('lon','lat')).values  #TODO: Check this. Maybe do it manually
data_xarray.close()

# Calculate lat and lon bounds
lat_bounds,lon_bounds = utils.bounding_latlon(lat,lon)

# Format the variables
var_spatial_mean = np.expand_dims(np.swapaxes(var_spatial_mean,0,2),axis=0)
var_global_mean  = np.expand_dims(var_global_mean,axis=0)

# OWDA only has one ensemble member, so the mean and the ensemble member are the same.
var_spatial_members = np.expand_dims(var_spatial_mean,axis=1)
var_global_members  = np.expand_dims(var_global_mean, axis=1)

# Get other metadata
methods = ['OWDA']
ens_spatial = np.array([1])
ens_global  = np.array([1])

# If this data can't be reformatted to the standard format, add a note here 
notes = ['OWDA only has one ensemble member, so the mean and the ensemble member are the same.']

# Check the shape of the variables
print(var_spatial_members.shape)
print(var_spatial_mean.shape)
print(var_global_members.shape)
print(var_global_mean.shape)


#%% FORMAT DATA

# Create new array
data_xarray_output = xr.Dataset(
    {
        'pdsi_global_mean':    (['method','age'],                          var_global_mean,    {'units':'PDSI'}),
        'pdsi_global_members': (['method','ens_global','age'],             var_global_members, {'units':'PDSI'}),
        'pdsi_spatial_mean':   (['method','age','lat','lon'],              var_spatial_mean,   {'units':'PDSI'}),
        'pdsi_spatial_members':(['method','ens_spatial','age','lat','lon'],var_spatial_members,{'units':'PDSI'})
    },
    coords={
        'method':     (['method'],methods),
        'notes':      (['notes'],notes),
        'ens_global': (['ens_global'],ens_global,{'description':'ensemble members'}),
        'ens_spatial':(['ens_spatial'],ens_spatial,{'description':'ensemble members'}),
        'age':        (['age'],age,{'units':'yr BP'}),
        'lat':        (['lat'],lat,{'units':'degrees_north'}),
        'lon':        (['lon'],lon,{'units':'degrees_east'}),
        'lat_bounds': (['lat_bounds'],lat_bounds,{'units':'degrees_north'}),
        'lon_bounds': (['lon_bounds'],lon_bounds,{'units':'degrees_east'}),
    },
    attrs={
        'dataset_name':      'Old World Drought Atlas',
        'dataset_source_url':'https://www.ncei.noaa.gov/access/paleo-search/study/19419',
    },
)


#%% SAVE DATA

# Save new array
output_dir = '/projects/pd_lab/data/paleoclimate_reconstructions/presto_format/'
data_xarray_output.to_netcdf(output_dir+'owda_v1_0_0_pdsi_jja.nc')

