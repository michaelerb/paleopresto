#==============================================================================
# Print the coordinates of the ipcc regions, to use in the html of the
# visualizer webpage.
#    author: Michael P. Erb
#    date  : 5/4/2023
#==============================================================================

import numpy as np
import regionmask


# Settings
coordinates_to_print = 'all'   # 'all' or 'corners'
precision = 'two_decimal'  # 'all_decimal' or 'two_decimal'


#%% LOAD WGI REGIONS

# Get all WGI regions
ar6_all = regionmask.defined_regions.ar6.all
with open('ipcc_region_code_'+coordinates_to_print+'_'+precision+'.txt','w') as f:
    f.write('')


#%% FUNCTIONS

# A function to turn the shapefile coords into html code.
def print_coords(region_lons,region_lats):
    #
    # Output either the full set of coordinates, or just the "corners"
    if coordinates_to_print == 'all':
        region_lats_to_output = region_lats
        region_lons_to_output = region_lons
    elif coordinates_to_print == 'corners':
        #
        # Remove redundant points
        n_points = len(region_lons)
        latlon_angle_last = np.nan
        ind_useful = []
        for j in range(n_points-1):
            if region_lons[j+1]-region_lons[j] == 0: latlon_angle_current = np.nan
            else: latlon_angle_current = (region_lats[j+1]-region_lats[j])/(region_lons[j+1]-region_lons[j])
            if np.isnan(latlon_angle_current) & np.isnan(latlon_angle_last): pass
            elif ~np.isclose(latlon_angle_current,latlon_angle_last): ind_useful.append(j)
            latlon_angle_last = latlon_angle_current
        #
        ind_useful.append(j+1)  # Add the last point
        print(i,n_points,'->',len(ind_useful))
        region_lats_to_output = np.array(region_lats)[ind_useful]
        region_lons_to_output = np.array(region_lons)[ind_useful]
    #
    # Save lat/lon strings, either with full precision or reduced precision (for reducing the size of the webpage)
    latlons_all = []
    for j in range(len(region_lats_to_output)):
        if precision == 'all_decimal':
            lat_txt = str(region_lats_to_output[j]).rstrip('0').rstrip('.')
            lon_txt = str(region_lons_to_output[j]).rstrip('0').rstrip('.')
        elif precision == 'two_decimal':
            lat_txt = str('{:.2f}'.format(region_lats_to_output[j])).rstrip('0').rstrip('.')
            lon_txt = str('{:.2f}'.format(region_lons_to_output[j])).rstrip('0').rstrip('.')
        latlon_txt = '['+lat_txt+','+lon_txt+']'
        latlons_all.append(latlon_txt)
    #
    # Make a file with the lat/lon values
    latlon_all_str = ",".join(latlons_all)
    latlon_all_str = "      ipcc_regions["+str(counter)+"] = ["+latlon_all_str+"];"
    with open('ipcc_region_code_'+coordinates_to_print+'_'+precision+'.txt','a') as f:
        f.write(latlon_all_str+'\n')


#%% LOOP THROUGH REGIONS

# Make website code for clickable region map
n_regions = len(ar6_all)
region_abbrev_all = []
counter = 0
for i in range(n_regions):
    region_number = str(i+1).zfill(2)
    region_abbrev = str(ar6_all.abbrevs[i])
    region_name   = str(ar6_all.names[i])
    region_shape  = ar6_all[i]._polygon
    region_abbrev_all.append(region_abbrev)
    if region_shape.geom_type == 'Polygon':
        region_lons,region_lats = region_shape.exterior.coords.xy
        print_coords(region_lons,region_lats)
        counter += 1
    elif region_shape.geom_type == 'MultiPolygon':
        for subregion_shape in region_shape.geoms:
            region_lons,region_lats = subregion_shape.exterior.coords.xy
            print_coords(region_lons,region_lats)
            counter += 1
    else: print(' === ERROR: Variable is not a Polygon or MultiPolygon ===')

# Format the list of regions
region_abbrev_str = "','".join(region_abbrev_all)
region_abbrev_str = "      var ipcc_region_names = ['"+region_abbrev_str+"'];"

with open('ipcc_region_names.txt','w') as f:
    f.write(region_abbrev_str)

