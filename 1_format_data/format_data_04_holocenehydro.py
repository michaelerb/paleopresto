#==============================================================================
# Make a standardized netcdf file for the Holocene Hydroclimate composites.
#    author: Michael P. Erb
#    date  : 11/8/2023
#=============================================================================

import sys
import numpy as np
import xarray as xr
import glob
import regionmask
from pyproj import Proj
from shapely.geometry import shape


#%% SETTINGS

#var_txt = 'tas'
#var_txt = 'hydro'
var_txt = sys.argv[1]

print('=== Processing Holocene Hydroclimate ===')


#%% COMPUTE AREAS

# Get all WGI regions
ar6_all = regionmask.defined_regions.ar6.all
ar6_abbreviations = ar6_all.abbrevs
n_regions = len(ar6_abbreviations)

# A function to calculate the area of each region
def calculate_area(region_coords):
    region_lons,region_lats = region_coords
    proj_equalarea = Proj("+proj=aea +lat_1=29.5 +lat_2=42.5")
    region_x,region_y = proj_equalarea(region_lons,region_lats)
    region_new_polygon = {'type':'polygon','coordinates':[zip(region_x,region_y)]}
    area_region = shape(region_new_polygon).area
    area_region = area_region / 1000000  # Convert m2 to km2
    return area_region

# Loop through the regions, calculating the area of each
ar6_area = np.zeros(n_regions); ar6_area[:] = np.nan
for i in range(n_regions):
    region_shape = ar6_all[i]._polygon
    if region_shape.geom_type == 'Polygon':
        region_area = calculate_area(region_shape.exterior.coords.xy)
    elif region_shape.geom_type == 'MultiPolygon':
        region_area = 0
        for subregion_shape in region_shape.geoms: region_area = region_area + calculate_area(subregion_shape.exterior.coords.xy)
    #
    ar6_area[i] = region_area

print('Total area:',sum(ar6_area)/1000000,'million km2')  # This should be approximately 510.1 million km2


#%% SET UP

# Set ages and ens
age = np.arange(0,12001,100)
ens = np.arange(1,501)

# Get dimensions
n_methods = 1
n_ens     = len(ens)
n_ages    = len(age)


#%% LOAD DATA
# Note: data can be downloaded at: https://zenodo.org/record/7939488

data_dir = '/projects/pd_lab/mpe32/HoloceneHydroclimate/Data/RegionComposites/clhancock-HoloceneHydroclimate-6094b85/Data/RegionComposites/'

# Get the names of the available files
if   var_txt == 'tas':   data_subdir = 'T/'
elif var_txt == 'hydro': data_subdir = 'HC/'
filenames_all = glob.glob(data_dir+data_subdir+'*.csv')
filenames_regions = [filename.split('/')[-1].split('.')[0] for filename in filenames_all if filename[-12:] != 'byRegion.csv']

# Load regional composites
var_spatial_members = np.zeros((n_methods,n_ens,n_ages,n_regions,1)); var_spatial_members[:] = np.nan
for i,region_name in enumerate(ar6_abbreviations):
    if region_name not in filenames_regions: print(i,region_name,'NO FILE FOR '+var_txt+' DATA')
    else:
        var_for_region = np.loadtxt(data_dir+data_subdir+region_name+'.csv',dtype=str,delimiter=',',skiprows=1)
        var_for_region[var_for_region == 'NA'] = np.nan
        var_for_region = var_for_region.astype(float)
        var_spatial_members[0,:,:,i,0] = np.swapaxes(var_for_region,0,1)


#%% CALCULATIONS

# Compute mean of ensemble members
var_spatial_mean = np.nanmean(var_spatial_members,axis=1)

# Compute global means
var_global_members = np.zeros((n_methods,n_ens,n_ages)); var_global_members[:] = np.nan
for i in range(n_ages):
    for j in range(n_ens):
        ind_valid_hc = np.isfinite(var_spatial_members[0,j,i,:,0])
        var_global_members[:,j,i] = np.average(var_spatial_members[0,j,i,ind_valid_hc,0],axis=0,weights=ar6_area[ind_valid_hc])

# Compute mean of global ensemble members
var_global_mean = np.nanmean(var_global_members,axis=1)

# Get other metadata
methods = ['Holocene Hydroclimate']
ens_spatial = ens
ens_global  = ens

# If this data can't be reformatted to the standard format, add a note here 
notes = ['lat_is_ar6_region_names']

# Check the shape of the variables
print(var_spatial_members.shape)
print(var_spatial_mean.shape)
print(var_global_members.shape)
print(var_global_mean.shape)


#%% FORMAT DATA

# Create new array
message_txt = np.array(['Use IPCC AR4 regions'])
data_xarray_output = xr.Dataset(
    {
        var_txt+'_global_mean':    (['method','age'],                          var_global_mean,    {'units':'Z-score'}),
        var_txt+'_global_members': (['method','ens_global','age'],             var_global_members, {'units':'Z-score'}),
        var_txt+'_spatial_mean':   (['method','age','lat','lon'],              var_spatial_mean,   {'units':'Z-score'}),
        var_txt+'_spatial_members':(['method','ens_spatial','age','lat','lon'],var_spatial_members,{'units':'Z-score'})
    },
    coords={
        'method':     (['method'],methods),
        'notes':      (['notes'],notes),
        'ens_global': (['ens_global'],ens_global,{'description':'ensemble members'}),
        'ens_spatial':(['ens_spatial'],ens_spatial,{'description':'selected ensemble members'}),
        'age':        (['age'],age,{'units':'yr BP'}),
        'lat':        (['lat'],ar6_abbreviations,{'units':'none'}),
        'lon':        (['lon'],message_txt,{'units':'none'}),
        'lat_bounds': (['lat_bounds'],message_txt,{'units':'none'}),
        'lon_bounds': (['lon_bounds'],message_txt,{'units':'none'}),
    },
    attrs={
        'dataset_name':      'Holocene Hydroclimate',
        'dataset_source_url':'https://zenodo.org/record/7939488',
    },
)


#%% SAVE DATA

# Save new array
output_dir = '/projects/pd_lab/data/paleoclimate_reconstructions/presto_format/'
data_xarray_output.to_netcdf(output_dir+'holocenehydro_v1_0_0_'+var_txt+'_annual.nc')

