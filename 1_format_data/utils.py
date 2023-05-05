#==============================================================================
# Functions for use in the PReSto project.
#    author: Michael P. Erb
#    date  : 5/4/2023
#==============================================================================

import numpy as np
import copy

#%% Compute the bounding latitudes and longitudes for a regular grid, given central values
def bounding_latlon(lat,lon):
    #
    # Calculate lon bounds
    lon_extended = copy.deepcopy(lon)
    lon_extended = np.insert(lon_extended,0,lon[0]-(lon[1]-lon[0]))
    lon_extended = np.append(lon_extended,lon[-1]+(lon[-1]-lon[-2]))
    lon_bounds = (lon_extended[:-1] + lon_extended[1:])/2
    #
    # Calculate lat bounds
    lat_extended = copy.deepcopy(lat)
    lat_extended = np.insert(lat_extended,0,lat[0]-(lat[1]-lat[0]))
    lat_extended = np.append(lat_extended,lat[-1]+(lat[-1]-lat[-2]))
    lat_bounds = (lat_extended[:-1] + lat_extended[1:])/2
    if max(lat_bounds) > 90:
        print('Maximum latitude bound is '+str(max(lat_bounds))+'. Setting to 90')
        lat_bounds[lat_bounds > 90]  = 90
    if min(lat_bounds) < -90:
        print('Minimum latitude bound is '+str(min(lat_bounds))+'. Setting to -90')
        lat_bounds[lat_bounds < -90] = -90
    #
    return lat_bounds,lon_bounds
