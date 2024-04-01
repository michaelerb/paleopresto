#==============================================================================
# Make a standardized netCDF file for MXDA.
#    author: Michael P. Erb
#    date  : 3/29/2024
#==============================================================================

import numpy as np
import xarray as xr
import utils
import pandas as pd


#%% LOAD DATA
# Note: data can be downloaded at: https://www.ncei.noaa.gov/access/paleo-search/study/20353

print('=== Processing MXDA ===')
data_dir = '/projects/pd_lab/data/paleoclimate_reconstructions/drought_atlases/'

# Load data
data_csv   = pd.read_csv(data_dir+'mxda-recon-jja-scpdsi-1400-2012.txt',delim_whitespace=True)
latlon_csv = pd.read_csv(data_dir+'mxda-grid-points.txt',delim_whitespace=True,header=None)
years = np.array(data_csv['year__#'])
lon_data = np.array(latlon_csv[0])
lat_data = np.array(latlon_csv[1])
age = 1950-years


#%% PUT DATA ON A REGULAR GRID

print(min(lat_data),max(lat_data))
print(min(lon_data),max(lon_data))

lat = np.arange(14.25,35,.5)
lon = np.arange(-120.25,-74.5,.5)
var_spatial_mean = np.zeros((len(years),len(lat),len(lon))); var_spatial_mean[:] = np.nan
for i in range(len(lon_data)):
    print('Gridding the data: '+str(i+1)+'/'+str(len(lon_data)))
    data_column = np.array(data_csv[str(i+1)])
    ind_lat = np.where(lat_data[i] == lat)[0][0]
    ind_lon = np.where(lon_data[i] == lon)[0][0]
    var_spatial_mean[:,ind_lat,ind_lon] = data_column


#%% CALCULATIONS

# Compute the spatial mean
var_global_mean = np.zeros((len(years))); var_global_mean[:] = np.nan
lat_weights = np.cos(np.deg2rad(lat_data))
for i in range(len(years)):
    data_for_year = np.array(data_csv.iloc[i][1:])
    ind_valid = np.isfinite(data_for_year)
    var_global_mean[i] = np.average(data_for_year[ind_valid],weights=lat_weights[ind_valid])

# Calculate lat and lon bounds
lat_bounds,lon_bounds = utils.bounding_latlon(lat,lon)

# Format the variables
var_spatial_mean = np.expand_dims(var_spatial_mean,axis=0)
var_global_mean  = np.expand_dims(var_global_mean, axis=0)

# MXDA only has one ensemble member, so the mean and the ensemble member are the same.
var_spatial_members = np.expand_dims(var_spatial_mean,axis=1)
var_global_members  = np.expand_dims(var_global_mean, axis=1)

# Get other metadata
methods = ['MXDA']
ens_spatial = np.array([1])
ens_global  = np.array([1])

# If this data can't be reformatted to the standard format, add a note here 
notes = ['MXDA only has one ensemble member, so the mean and the ensemble member are the same.']

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
        'dataset_name':      'Mexican Drought Atlas',
        'dataset_source_url':'https://www.ncei.noaa.gov/access/paleo-search/study/20353',
    },
)


#%% SAVE DATA

# Save new array
output_dir = '/projects/pd_lab/data/paleoclimate_reconstructions/presto_format/'
data_xarray_output.to_netcdf(output_dir+'mxda_v1_0_0_pdsi_jja.nc')

