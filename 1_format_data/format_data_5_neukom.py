#==============================================================================
# Make a standardized netcdf file for the Neukom et al., 2019 reconstructions.
#    author: Michael P. Erb
#    date  : 5/4/2023
#=============================================================================

import numpy as np
import xarray as xr
import utils


#%% LOAD DATA
# Note: data can be downloaded at: https://www.ncei.noaa.gov/access/paleo-search/study/26850

print('=== Processing Neukom et al., 2019 reconstructions ===')

# Load data
data_dir = '/projects/pd_lab/data/paleoclimate_reconstructions/spatial_2k_reconstructions/'
data_xarray_am      = xr.open_dataset(data_dir+'AM.nc')
data_xarray_cca     = xr.open_dataset(data_dir+'CCA.nc')
data_xarray_cps     = xr.open_dataset(data_dir+'CPS.nc')
data_xarray_da      = xr.open_dataset(data_dir+'DA.nc')
data_xarray_graphem = xr.open_dataset(data_dir+'GraphEM.nc')
data_xarray_pcr     = xr.open_dataset(data_dir+'PCR.nc')

# Compute global means for Neukom reconstructions
lat_weights = np.cos(np.deg2rad(data_xarray_am.lat))
tas_global_am      = data_xarray_am.var0.weighted(lat_weights).mean(('lon','lat'))
tas_global_cca     = data_xarray_cca.tas.weighted(lat_weights).mean(('lon','lat'))
tas_global_cps     = data_xarray_cps.tas.weighted(lat_weights).mean(('lon','lat'))
tas_global_da      = data_xarray_da.temp.weighted(lat_weights).mean(('lon','lat'))
tas_global_graphem = data_xarray_graphem.temp.weighted(lat_weights).mean(('lon','lat'))
tas_global_pcr     = data_xarray_pcr.tas.weighted(lat_weights).mean(('lon','lat'))

# Get coordinates
year        = data_xarray_da['time'].values
lat         = data_xarray_da['lat'].values
lon         = data_xarray_da['lon'].values
ens_spatial = data_xarray_da['ens'].values
ens_global  = ens_spatial

# Set the ages
age = 1950-year


#%% CALCULATIONS

# Calculate lat and lon bounds
lat_bounds,lon_bounds = utils.bounding_latlon(lat,lon)

# Get other metadata
methods = ['AM','CCA','CPS','DA','GraphEM','PCR']

# Get dimensions
n_methods = len(methods)
n_ens     = len(ens_spatial)
n_ages    = len(age)
n_lat     = len(lat)
n_lon     = len(lon)

# Create a single spatial variable with the chosen dimensions
tas_ens = np.zeros((n_methods,n_ens,n_ages,n_lat,n_lon)); tas_ens[:] = np.nan
tas_ens[0,:,:,:,:] = np.swapaxes(data_xarray_am['var0'].values,0,1)
tas_ens[1,:,:,:,:] = np.swapaxes(data_xarray_cca['tas'].values,0,1)
tas_ens[2,:,:,:,:] = np.swapaxes(data_xarray_cps['tas'].values,0,1)
tas_ens[3,:,:,:,:] = np.swapaxes(data_xarray_da['temp'].values,0,1)
tas_ens[4,:,:,:,:] = np.swapaxes(data_xarray_graphem['temp'].values,0,1)
tas_ens[5,:,:,:,:] = np.swapaxes(data_xarray_pcr['tas'].values,0,1)

# Create a single spatial variable with the chosen dimensions
tas_global = np.zeros((n_methods,n_ens,n_ages)); tas_global[:] = np.nan
tas_global[0,:,:] = np.swapaxes(tas_global_am.values,0,1)
tas_global[1,:,:] = np.swapaxes(tas_global_cca.values,0,1)
tas_global[2,:,:] = np.swapaxes(tas_global_cps.values,0,1)
tas_global[3,:,:] = np.swapaxes(tas_global_da.values,0,1)
tas_global[4,:,:] = np.swapaxes(tas_global_graphem.values,0,1)
tas_global[5,:,:] = np.swapaxes(tas_global_pcr.values,0,1)

# Create an average across ensemble members
tas_mean = np.mean(tas_ens,axis=1)

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
output_name = 'neukom2019_v1_0_0_tas_annual.nc'
data_xarray_output.to_netcdf(output_dir+output_name)
