#==============================================================================
# Make a standardized netCDF file for the LGMR.
#    author: Michael P. Erb
#    date  : 8/25/2023
#==============================================================================

import sys
import numpy as np
import xarray as xr
import utils


#%% SETTINGS

season_txt = 'annual'
#season_txt = 'jja'
#season_txt = 'djf'
#season_txt = sys.argv[1]


#%% LOAD DATA
# Note: data can be downloaded at: https://www.ncei.noaa.gov/access/paleo-search/study/33112

print('=== Processing PHYDA ===')
data_dir = '/projects/pd_lab/data/paleoclimate_reconstructions/PHYDA/'
if   season_txt == 'annual': data_file = 'da_hydro_AprMar_r.1-2000_d.05-Jan-2018.nc'
elif season_txt == 'jja':    data_file = 'da_hydro_JunAug_r.1-2000_d.05-Jan-2018.nc'
elif season_txt == 'djf':    data_file = 'da_hydro_DecFeb_r.1-2000_d.05-Jan-2018.nc'

# Load data
data_xarray = xr.open_dataset(data_dir+data_file)
tas_mean   = data_xarray['tas_mn'].values
pdsi_mean  = data_xarray['pdsi_mn'].values
spei_mean  = data_xarray['spei_mn'].values
tas_ens    = data_xarray['tas_pc'].values  # The 5th, 50th, and 95th percentiles
pdsi_ens   = data_xarray['pdsi_pc'].values
spei_ens   = data_xarray['spei_pc'].values
tas_units  = data_xarray['tas_mn'].units
pdsi_units = data_xarray['pdsi_mn'].units
spei_units = data_xarray['spei_mn'].units
years      = data_xarray['time'].values
lat        = data_xarray['lat'].values
lon        = data_xarray['lon'].values

tas_global    = data_xarray['gmt_mn'].values
tas_global_pc = data_xarray['gmt_pc'].values

# Compute the global mean of d18Op
lat_weights = np.cos(np.deg2rad(data_xarray.lat))
pdsi_global = data_xarray.pdsi_mn.weighted(lat_weights).mean(('lon','lat'))
spei_global = data_xarray.spei_mn.weighted(lat_weights).mean(('lon','lat'))
data_xarray.close()

age = 1950 - years



#%% FIGURE
"""
import matplotlib.pyplot as plt

# Global mean temperature
plt.figure(figsize=(20,14))
ax1 = plt.subplot2grid((1,1),(0,0))
#line1, = ax1.plot(years,tas_global,color='k',linewidth=3,label='Mean')
#range1 = ax1.fill_between(years,tas_global_pc[0,:],tas_global_pc[2,:],facecolor='k',alpha=0.15,label='95% range')
#line1, = ax1.plot(years,pdsi_global,color='k',linewidth=3,label='Mean')
line1, = ax1.plot(years,spei_global,color='k',linewidth=3,label='Mean')
ax1.axhline(y=0,color='k',linewidth=1,linestyle='dashed',alpha=0.5)
ax1.set_ylabel('Global mean $\Delta$T ($^\circ$C)',fontsize=20)
ax1.tick_params(axis='both',which='major',labelsize=20)
plt.show()
"""


#%% CALCULATIONS

# Calculate lat and lon bounds
lat_bounds,lon_bounds = utils.bounding_latlon(lat,lon)

# Format the variables
tas_ens     = np.expand_dims(tas_ens,    axis=0)
tas_mean    = np.expand_dims(tas_mean,   axis=0)
tas_global  = np.expand_dims(np.expand_dims(tas_global, axis=0),axis=0)
pdsi_ens    = np.expand_dims(pdsi_ens,   axis=0)
pdsi_mean   = np.expand_dims(pdsi_mean,  axis=0)
pdsi_global = np.expand_dims(np.expand_dims(pdsi_global,axis=0),axis=0)
spei_ens    = np.expand_dims(spei_ens,   axis=0)
spei_mean   = np.expand_dims(spei_mean,  axis=0)
spei_global = np.expand_dims(np.expand_dims(spei_global,axis=0),axis=0)

# Set other metadata
methods = ['PHYDA']
ens_spatial = np.arange(tas_ens.shape[1])+1
ens_global  = np.array([1])

# If this data can't be reformatted to the standard format, add a note here 
notes = ['Ensemble members are actually the 5th, 50th, and 95th percentiles']


#%% FORMAT DATA

# Create new arrays
data_xarray_output_tas = xr.Dataset(
    {
        'tas_global':(['method','ens_global','age'],             tas_global,{'units':tas_units}),
        'tas_mean':  (['method','age','lat','lon'],              tas_mean,  {'units':tas_units}),
        'tas_ens':   (['method','ens_spatial','age','lat','lon'],tas_ens,   {'units':tas_units})
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


data_xarray_output_pdsi = xr.Dataset(
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


data_xarray_output_spei = xr.Dataset(
    {
        'spei_global':(['method','ens_global','age'],             spei_global,{'units':spei_units}),
        'spei_mean':  (['method','age','lat','lon'],              spei_mean,  {'units':spei_units}),
        'spei_ens':   (['method','ens_spatial','age','lat','lon'],spei_ens,   {'units':spei_units})
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
data_xarray_output_tas.to_netcdf(output_dir+'phyda_v1_0_0_tas_'+season_txt+'.nc')
data_xarray_output_pdsi.to_netcdf(output_dir+'phyda_v1_0_0_pdsi_'+season_txt+'.nc')
data_xarray_output_spei.to_netcdf(output_dir+'phyda_v1_0_0_spei_'+season_txt+'.nc')
