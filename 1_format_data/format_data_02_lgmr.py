#==============================================================================
# Make a standardized netCDF file for the LGMR.
#    author: Michael P. Erb
#    date  : 9/12/2023
#==============================================================================

import sys
import numpy as np
import xarray as xr
import utils


#%% SETTINGS

#var_txt = 'tas'
#var_txt = 'd18Op'
var_txt = sys.argv[1]


#%% LOAD DATA
# Note: data can be downloaded at: https://www.ncei.noaa.gov/access/paleo-search/study/33112

print('=== Processing LGMR ===')
data_dir = '/projects/pd_lab/data/paleoclimate_reconstructions/Holocene_reconstructions/Osman_etal_2021/'

# Load data
if var_txt == 'tas':
    #
    data_xarray = xr.open_dataset(data_dir+'LGMR_SAT_climo.nc')
    var_spatial_mean = data_xarray['sat'].values
    data_xarray.close()
    #
    data_xarray = xr.open_dataset(data_dir+'LGMR_SAT_ens.nc')
    var_spatial_members = data_xarray['sat'].values
    var_units           = data_xarray['sat'].units
    ens_spatial         = data_xarray['nEns'].values
    age = data_xarray['age'].values
    lat = data_xarray['lat'].values
    lon = data_xarray['lon'].values
    data_xarray.close()
    #
    data_xarray = xr.open_dataset(data_dir+'LGMR_GMST_ens.nc')
    var_global_members = data_xarray['gmst'].values
    ens_global         = data_xarray['nEns'].values
    data_xarray.close()
    #
elif var_txt == 'd18Op':
    #
    data_xarray = xr.open_dataset(data_dir+'LGMR_d18Op_climo.nc')
    var_spatial_mean = data_xarray['d18Op'].values
    data_xarray.close()
    #
    data_xarray = xr.open_dataset(data_dir+'LGMR_d18Op_ens.nc')
    var_spatial_members = data_xarray['d18Op'].values
    var_units           = data_xarray['d18Op'].units
    ens_spatial         = data_xarray['nEns'].values
    age = data_xarray['age'].values
    lat = data_xarray['lat'].values
    lon = data_xarray['lon'].values
    #
    # Compute the global mean of d18Op
    lat_weights  = np.cos(np.deg2rad(data_xarray.lat))
    var_global_members = data_xarray.d18Op.weighted(lat_weights).mean(('lon','lat'))
    ens_global         = ens_spatial
    data_xarray.close()


#%% CALCULATIONS

# Calculate lat and lon bounds
lat_bounds,lon_bounds = utils.bounding_latlon(lat,lon)

# Format the variables
var_spatial_members = np.expand_dims(var_spatial_members,   axis=0)
var_spatial_mean    = np.expand_dims(var_spatial_mean,  axis=0)
var_global_members  = np.expand_dims(var_global_members,axis=0)
var_global_mean     = np.mean(var_global_members,axis=1)

# Get other metadata
methods = ['LGMR']

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
        'ens_spatial':(['ens_spatial'],ens_spatial,{'description':'ensemble members'}),
        'age':        (['age'],age,{'units':'yr BP'}),
        'lat':        (['lat'],lat,{'units':'degrees_north'}),
        'lon':        (['lon'],lon,{'units':'degrees_east'}),
        'lat_bounds': (['lat_bounds'],lat_bounds,{'units':'degrees_north'}),
        'lon_bounds': (['lon_bounds'],lon_bounds,{'units':'degrees_east'}),
    },
    attrs={
        'dataset_name':      'LGMR',
        'dataset_source_url':'https://www.ncei.noaa.gov/access/paleo-search/study/33112',
    },
)


#%% SAVE DATA

# Save new array
output_dir = '/projects/pd_lab/data/paleoclimate_reconstructions/presto_format/'
data_xarray_output.to_netcdf(output_dir+'lgmr_v1_0_0_'+var_txt+'_annual.nc')

