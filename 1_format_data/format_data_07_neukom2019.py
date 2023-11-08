#==============================================================================
# Make a standardized netcdf file for the Neukom et al., 2019 reconstructions.
#    author: Michael P. Erb
#    date  : 9/12/2023
#=============================================================================

import numpy as np
import xarray as xr
import utils


#%% LOAD SPATIAL
# Note: data can be downloaded at: https://www.ncei.noaa.gov/access/paleo-search/study/26850

print('=== Processing Neukom et al., 2019 reconstructions ===')

# Load data
data_dir = '/projects/pd_lab/data/paleoclimate_reconstructions/spatial_2k_reconstructions/'
data_xarray_am      = xr.open_dataset(data_dir+'AM.nc',     use_cftime=True)
data_xarray_cca     = xr.open_dataset(data_dir+'CCA.nc',    use_cftime=True)
data_xarray_cps     = xr.open_dataset(data_dir+'CPS.nc',    use_cftime=True)
data_xarray_da      = xr.open_dataset(data_dir+'DA.nc',     use_cftime=True)
data_xarray_graphem = xr.open_dataset(data_dir+'GraphEM.nc',use_cftime=True)
data_xarray_pcr     = xr.open_dataset(data_dir+'PCR.nc',    use_cftime=True)

# Get coordinates
year        = data_xarray_da['time'].values
lat         = data_xarray_da['lat'].values
lon         = data_xarray_da['lon'].values
ens_spatial = data_xarray_da['ens'].values
ens_global  = ens_spatial
age = 1950-year

# Get other metadata
methods = ['AM','CCA','CPS','DA','GraphEM','PCR']

# Get dimensions
n_methods = len(methods)
n_ens     = len(ens_spatial)
n_ages    = len(age)
n_lat     = len(lat)
n_lon     = len(lon)

# Create a single spatial variable with the chosen dimensions
var_spatial_members = np.zeros((n_methods,n_ens,n_ages,n_lat,n_lon)); var_spatial_members[:] = np.nan
var_spatial_members[0,:,:,:,:] = np.swapaxes(data_xarray_am['var0'].values,0,1)
var_spatial_members[1,:,:,:,:] = np.swapaxes(data_xarray_cca['tas'].values,0,1)
var_spatial_members[2,:,:,:,:] = np.swapaxes(data_xarray_cps['tas'].values,0,1)
var_spatial_members[3,:,:,:,:] = np.swapaxes(data_xarray_da['temp'].values,0,1)
var_spatial_members[4,:,:,:,:] = np.swapaxes(data_xarray_graphem['temp'].values,0,1)
var_spatial_members[5,:,:,:,:] = np.swapaxes(data_xarray_pcr['tas'].values,0,1)


#%% LOAD GLOBAL DATA

# Compute global means for Neukom reconstructions
lat_weights = np.cos(np.deg2rad(data_xarray_am.lat))
var_global_am      = data_xarray_am.var0.weighted(lat_weights).mean(('lon','lat'))
var_global_cca     = data_xarray_cca.tas.weighted(lat_weights).mean(('lon','lat'))
var_global_cps     = data_xarray_cps.tas.weighted(lat_weights).mean(('lon','lat'))
var_global_da      = data_xarray_da.temp.weighted(lat_weights).mean(('lon','lat'))
var_global_graphem = data_xarray_graphem.temp.weighted(lat_weights).mean(('lon','lat'))
var_global_pcr     = data_xarray_pcr.tas.weighted(lat_weights).mean(('lon','lat'))

# Create a single spatial variable with the chosen dimensions
var_global_members = np.zeros((n_methods,n_ens,n_ages)); var_global_members[:] = np.nan
var_global_members[0,:,:] = np.swapaxes(var_global_am.values,0,1)
var_global_members[1,:,:] = np.swapaxes(var_global_cca.values,0,1)
var_global_members[2,:,:] = np.swapaxes(var_global_cps.values,0,1)
var_global_members[3,:,:] = np.swapaxes(var_global_da.values,0,1)
var_global_members[4,:,:] = np.swapaxes(var_global_graphem.values,0,1)
var_global_members[5,:,:] = np.swapaxes(var_global_pcr.values,0,1)


#%% CALCULATIONS

# Calculate lat and lon bounds
lat_bounds,lon_bounds = utils.bounding_latlon(lat,lon)

# Create an average across ensemble members
var_spatial_mean = np.mean(var_spatial_members,axis=1)
var_global_mean  = np.mean(var_global_members,axis=1)

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
        'ens_global': (['ens_global'],ens_global,{'description':'ensemble members'}),
        'ens_spatial':(['ens_spatial'],ens_spatial,{'description':'ensemble members'}),
        'age':        (['age'],age,{'units':'yr BP'}),
        'lat':        (['lat'],lat,{'units':'degrees_north'}),
        'lon':        (['lon'],lon,{'units':'degrees_east'}),
        'lat_bounds': (['lat_bounds'],lat_bounds,{'units':'degrees_north'}),
        'lon_bounds': (['lon_bounds'],lon_bounds,{'units':'degrees_east'}),
    },
    attrs={
        'dataset_name':      'Neukom et al., 2019',
        'dataset_source_url':'https://www.ncei.noaa.gov/access/paleo-search/study/26850',
    },
)


#%% SAVE DATA

# Save new array
output_dir = '/projects/pd_lab/data/paleoclimate_reconstructions/presto_format/'
data_xarray_output.to_netcdf(output_dir+'neukom2019_v1_0_0_tas_annual.nc')

