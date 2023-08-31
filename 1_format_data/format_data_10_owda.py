#==============================================================================
# Make a standardized netCDF file for OWDA.
#    author: Michael P. Erb
#    date  : 7/20/2023
#==============================================================================

import numpy as np
import xarray as xr
import utils


#%% LOAD DATA
# Note: data can be downloaded at: https://www.ncei.noaa.gov/access/paleo-search/study/19419

print('=== Processing NADA ===')
data_dir = '/projects/pd_lab/data/paleoclimate_reconstructions/drought_atlases/'

# Load data
data_xarray = xr.open_dataset(data_dir+'owda.nc')
pdsi_mean = data_xarray['pdsi'].values
years     = data_xarray['time'].values
lat       = data_xarray['lat'].values
lon       = data_xarray['lon'].values
data_xarray.close()
age = 1950-years


#%% CALCULATIONS

# Compute the spatial mean
lat_weights = np.cos(np.deg2rad(data_xarray.lat))
pdsi_global = data_xarray['pdsi'].weighted(lat_weights).mean(('lon','lat')).values  #TODO: Check this. Maybe do it manually
data_xarray.close()

# Calculate lat and lon bounds
lat_bounds,lon_bounds = utils.bounding_latlon(lat,lon)

# Format the variables
pdsi_mean   = np.swapaxes(pdsi_mean,0,2)
pdsi_mean   = np.expand_dims(pdsi_mean,axis=0)
pdsi_ens    = np.expand_dims(pdsi_mean,axis=1)
pdsi_global = np.expand_dims(np.expand_dims(pdsi_global,axis=0),axis=0)

# Get other metadata
methods = ['OWDA']
ens_spatial = np.array([1])
ens_global  = ens_spatial
pdsi_units = 'PDSI'

# If this data can't be reformatted to the standard format, add a note here 
notes = ['']


#%% FORMAT DATA

# Create new array
data_xarray_output = xr.Dataset(
    {
        'pdsi_global':(['method','ens_global','age'],             pdsi_global,{'units':pdsi_units}),
        'pdsi_mean':  (['method','age','lat','lon'],              pdsi_mean,  {'units':pdsi_units}),
        'pdsi_ens':   (['method','ens_spatial','age','lat','lon'],pdsi_ens,   {'units':pdsi_units})
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
)

#%% SAVE DATA

# Save new array
output_dir = '/projects/pd_lab/data/paleoclimate_reconstructions/presto_format/'
output_name = 'owda_v1_0_0_pdsi_jja.nc'
data_xarray_output.to_netcdf(output_dir+output_name)
