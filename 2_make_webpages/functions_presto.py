#==============================================================================
# Functions for Presto
#    author: Michael P. Erb
#    date  : 1/18/2024
#==============================================================================

import numpy as np
import copy


#%% GRID CALCULATIONS

def select_latlons(lat,lon,map_region,dataset_txt):
    #
    # Select the approximate grid to generate time series figures for
    if   map_region == 'global':    lat_grid_desired = 5; lon_grid_desired = 5
    elif map_region == 'n_america': lat_grid_desired = 2; lon_grid_desired = 2
    elif map_region == 'europe':    lat_grid_desired = 1; lon_grid_desired = 1
    #
    # If requested, figure out a reduced set of indices for generating figures
    if len(lat)*len(lon) > 1000:
        #
        # Get the indicies for points along the approximate grid
        lat_grid_mean = np.abs(np.mean(lat[1:]-lat[:-1]))
        lon_grid_mean = np.abs(np.mean(lon[1:]-lon[:-1]))
        lat_skip_factor = max(round(lat_grid_desired/lat_grid_mean),1)
        lon_skip_factor = max(round(lon_grid_desired/lon_grid_mean),1)
        lat_start_value = int(np.floor(np.remainder(len(lat)-1,lat_skip_factor)/2))
        lon_start_value = int(np.floor(np.remainder(len(lon)-1,lon_skip_factor)/2))
        if dataset_txt == 'era5': lon_start_value += 1
        j_for_ts = np.arange(lat_start_value,len(lat),lat_skip_factor)
        i_for_ts = np.arange(lon_start_value,len(lon),lon_skip_factor)
        #
        # Print some stats
        print('--- Skip factors for grid locations ---')
        print('Lat:',lat_skip_factor)
        print('Lon:',lon_skip_factor)
        print('Original number of locations:',len(lat)*len(lon))
        print('New number of locations:     ',len(lat[j_for_ts])*len(lon[i_for_ts]))
        #
    else:
        #
        j_for_ts = np.arange(len(lat))
        i_for_ts = np.arange(len(lon))
    #
    #
    #%% Save grid lats and lons
    #
    if (all(isinstance(x,np.integer) for x in lon)) or (all(isinstance(x,float) for x in lon)) or (all(isinstance(x,np.float32) for x in lon)) or (all(isinstance(x,np.float64) for x in lon)):
        #
        # Make a version of the longitude that goes from -180 to 180
        lon_neg = copy.deepcopy(lon)
        lon_neg[lon_neg > 180] = lon_neg[lon_neg > 180] - 360
        #
        # Make versions for printing
        lat_string = ','.join([str('{:.1f}'.format(value)) for value in lat[j_for_ts]])
        lon_string = ','.join([str('{:.1f}'.format(value)) for value in lon_neg[i_for_ts]])
        #
        lat_string = "          var lat_all = ["+lat_string+"];"
        lon_string = "          var lon_all = ["+lon_string+"];"
        #
    else:
        lat_string = "          var lat_all = [];"
        lon_string = "          var lon_all = [];"
    #
    return lat_string,lon_string,j_for_ts,i_for_ts,lon_neg

