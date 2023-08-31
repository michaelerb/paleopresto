#==============================================================================
# Make a standardized netCDF file for Holocene Reconstruction.
#    author: Michael P. Erb
#    date  : 5/5/2023
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
tas_global = data_xarray['recon_tas_global_mean'].values
tas_mean   = data_xarray['recon_tas_mean'].values
tas_ens    = data_xarray['recon_tas_ens'].values
age = data_xarray['ages'].values
lat = data_xarray['lat'].values
lon = data_xarray['lon'].values


#%% CALCULATIONS

# Compute the global mean
lat_weights = np.cos(np.deg2rad(data_xarray.lat))
var_global = data_xarray[load_var_txt].weighted(lat_weights).mean(('longitude','latitude')).values[:,0,:]
data_xarray.close()


# Calculate lat and lon bounds
lat_bounds,lon_bounds = utils.bounding_latlon(lat,lon)

# Format the variables
tas_ens    = np.swapaxes(tas_ens,   0,1)
tas_global = np.swapaxes(tas_global,0,1)
tas_ens    = np.expand_dims(tas_ens,   axis=0)
tas_mean   = np.expand_dims(tas_mean,  axis=0)
tas_global = np.expand_dims(tas_global,axis=0)

# Get other metadata
methods = ['Holocene Reconstruction']
ens_spatial = np.arange(tas_ens.shape[1])+1
ens_global  = np.arange(tas_global.shape[1])+1

# If this data can't be reformatted to the standard format, add a note here 
notes = ['']


#%% FORMAT DATA

# Create new array
data_xarray_output = xr.Dataset(
    {
        'tas_global':(['method','ens_global','age'],             tas_global,{'units':'degrees Celsius'}),
        'tas_mean':  (['method','age','lat','lon'],              tas_mean,  {'units':'degrees Celsius'}),
        'tas_ens':   (['method','ens_spatial','age','lat','lon'],tas_ens,   {'units':'degrees Celsius'})
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
)


#%% SAVE DATA

# Save new array
output_dir = '/projects/pd_lab/data/paleoclimate_reconstructions/presto_format/'
output_name = 'daholocene_v1_0_0_tas_annual.nc'
data_xarray_output.to_netcdf(output_dir+output_name)
