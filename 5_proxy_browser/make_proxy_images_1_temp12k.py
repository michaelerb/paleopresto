#==============================================================================
# This script loops through all of the Temp12k proxies and makes dashboards.
#    author: Michael Erb (and maybe Chris Hancock?)
#    date  : 12/20/2023
#==============================================================================

import sys
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import numpy  as np
import lipd
import pickle

plt.style.use('ggplot')
save_instead_of_plot = True

# Choose dataset and options
version_txt = '1_0_2'
time_span = sys.argv[1]
#time_span = '12ka'
#time_span = '21ka'
#time_span = 'all'


#%% LOAD DATA

# Set the necessary directories
data_dir   = '/projects/pd_lab/data/data_assimilation/proxies/temp12k/'
output_dir = '/home/mpe32/analysis/4_presto/1_website/4_proxy_browser/'

# Load the Temp12k proxy metadata
file_to_open = open(data_dir+'Temp12k'+version_txt+'.pkl','rb')
proxies_all_12k = pickle.load(file_to_open)['D']
file_to_open.close()

# Extract the time series and use only those which are in Temp12k and in units of degC
all_ts_12k = lipd.extractTs(proxies_all_12k)

# Fix the GISP2 ages - Note: this is a temporary fix, since lipd isn't loading the right ages.
for i in range(len(all_ts_12k)):
    if (all_ts_12k[i]['dataSetName'] == 'Alley.GISP2.2000') and (all_ts_12k[i]['paleoData_variableName'] == 'age'): gisp2_ages = all_ts_12k[i]['paleoData_values']

for i in range(len(all_ts_12k)):
    if (all_ts_12k[i]['dataSetName'] == 'Alley.GISP2.2000') and (all_ts_12k[i]['paleoData_variableName'] == 'temperature') and (np.max(np.array(all_ts_12k[i]['age']).astype(float)) < 50):
        print('Fixing GISP2 ages:',all_ts_12k[i]['paleoData_variableName'],', Index:',i)
        all_ts_12k[i]['age'] = gisp2_ages

filtered_ts = lipd.filterTs(all_ts_12k, 'paleoData_inCompilation == Temp12k')
filtered_ts = lipd.filterTs(filtered_ts,'paleoData_units == degC')

# Set up arrays for lat/lons/TSids
lats_all = []
lons_all = []
tsid_all = []
datasetname_all = []
archive_all     = []
proxy_all       = []

#%% FIGURES
# Loop through all proxies, making dashboards
n_proxies = len(filtered_ts)
i = 0
for i in range(n_proxies):
    #
    print('Making dashboard '+str(i+1)+'/'+str(n_proxies))
    #
    # Get data and metadata
    proxy_ages = np.array(filtered_ts[i]['age']).astype(np.float)
    proxy_data = np.array(filtered_ts[i]['paleoData_values']).astype(np.float)
    #
    # Get proxy metadata
    metadata = {}
    try:    metadata['Dataset name']   = filtered_ts[i]['dataSetName']
    except: metadata['Dataset name']   = ''
    try:    metadata['Archive type']   = filtered_ts[i]['archiveType']
    except: metadata['Archive type']   = ''
    try:    metadata['Proxy type']     = filtered_ts[i]['paleoData_proxy']
    except: metadata['Proxy type']     = ''
    try:    proxytype_general          = filtered_ts[i]['paleoData_proxyGeneral']
    except: proxytype_general          = ''
    try:    metadata['Variable name']  = filtered_ts[i]['paleoData_variableName']
    except: metadata['Variable name']  = ''
    try:    metadata['TSid']           = filtered_ts[i]['paleoData_TSid']
    except: metadata['TSid']           = ''
    try:    metadata['Latitude ($^\circ$N)']  = filtered_ts[i]['geo_meanLat']
    except: metadata['Latitude ($^\circ$N)']  = ''
    try:    metadata['Longitude ($^\circ$E)'] = filtered_ts[i]['geo_meanLon']
    except: metadata['Longitude ($^\circ$E)'] = ''
    try:    metadata['Elevation (m)']  = filtered_ts[i]['geo_meanElev']
    except: metadata['Elevation (m)']  = ''
    try:    metadata['Interp. variable'] = filtered_ts[i]['paleoData_interpretation'][0]['variable']
    except: metadata['Interp. variable'] = ''
    try:    metadata['Interp. detail'] = filtered_ts[i]['paleoData_interpretation'][0]['variableDetail']
    except: metadata['Interp. detail'] = ''
    try:    metadata['Interp. season'] = filtered_ts[i]['paleoData_interpretation'][0]['seasonalityGeneral']
    except: metadata['Interp. season'] = ''
    try:    metadata['Interp. direction'] = filtered_ts[i]['paleoData_interpretation'][0]['direction']
    except: metadata['Interp. direction'] = ''
    try:    metadata['Pub. 1 title']   = filtered_ts[i]['pub1_title']
    except: metadata['Pub. 1 title']   = ''
    try:    metadata['Pub. 1 authors'] = filtered_ts[i]['pub1_author']
    except: metadata['Pub. 1 authors'] = ''
    try:    metadata['Pub. 1 year']    = filtered_ts[i]['pub1_year']
    except: metadata['Pub. 1 year']    = ''
    try:    metadata['Pub. 1 DOI']     = filtered_ts[i]['pub1_doi']
    except: metadata['Pub. 1 DOI']     = ''
    try:    metadata['Original data URL'] = filtered_ts[i]['originalDataUrl']
    except: metadata['Original data URL'] = ''
    try:    data_units = filtered_ts[i]['paleoData_units']
    except: data_units = ''
    try:    age_units  = filtered_ts[i]['ageUnits']
    except: age_units  = ''
    #
    # Set region to plot
    proxy_lat = metadata['Latitude ($^\circ$N)']
    proxy_lon = metadata['Longitude ($^\circ$E)']
    padding = 10
    region_bounds = [proxy_lon-padding,proxy_lon+padding,proxy_lat-padding,proxy_lat+padding]
    #
    #
    # Specify subplots
    plt.figure(figsize=(24,5))
    ax1 = plt.subplot2grid((2,8),(0,0),rowspan=2,colspan=6)
    ax2 = plt.subplot2grid((2,8),(0,7),projection=ccrs.PlateCarree(central_longitude=proxy_lon))
    #
    # Make a time series of the proxy data
    plot_dots, = ax1.plot(proxy_ages,proxy_data,'-ob',color='tab:blue',linewidth=2,markersize=5)
    xlimits = ax1.get_xlim()
    #ax1.axvspan(12000,0,facecolor='gray',alpha=.2)  # Shade the 0-12ka region with a light gray
    if   time_span == '12ka': ax1.set_xlim(0,12000)
    elif time_span == '21ka': ax1.set_xlim(0,21000)
    else:                     ax1.set_xlim(xlimits)
    if len(metadata['TSid']) <= 50: tsid_str = metadata['TSid']
    else:                           tsid_str = metadata['TSid'][0:50]+'...'
    ax1.set_title(metadata['Archive type']+'  |  '+metadata['Proxy type']+'  |  '+metadata['Dataset name']+'  |  '+tsid_str,fontweight='bold',loc='left',fontsize=18)
    if   time_span == '12ka': xloc = 0.42
    elif time_span == '21ka': xloc = 0.47
    else:                     xloc = 0.5
    ax1.set_xlabel('Age ('+age_units+')',fontsize=14,labelpad=-15,x=xloc)
    ax1.set_ylabel(data_units,fontsize=14)
    ax1.tick_params(labelsize=14)
    ax1.invert_xaxis()
    if metadata['Interp. direction'] == 'negative': ax1.invert_yaxis()
    #
    # Make a map of the proxy location
    ax2.set_extent(region_bounds,ccrs.PlateCarree())
    ax2.scatter(metadata['Longitude ($^\circ$E)'],metadata['Latitude ($^\circ$N)'],300,marker='o',facecolor="None",edgecolor='k',linewidth=3,zorder=1,transform=ccrs.PlateCarree())
    ax2.add_feature(cfeature.LAND,facecolor='lightgray')
    ax2.coastlines(zorder=2)
    ax2.add_feature(cfeature.LAKES,facecolor='none',edgecolor='k')
    ax2.set_title('Proxy location',fontsize=14)
    #
    # Print metadata
    metadata_selected = ['Dataset name','Archive type','Proxy type','Variable name','Latitude ($^\circ$N)','Longitude ($^\circ$E)','Elevation (m)','Interp. variable',\
                         'Interp. detail','Interp. season','Interp. direction','Pub. 1 title','Pub. 1 authors','Pub. 1 year','Pub. 1 DOI','Original data URL']
    for offset,key in enumerate(metadata_selected):
        if len(str(metadata[key])) <= 40: metadata_str = str(metadata[key])
        else:                             metadata_str = str(metadata[key])[0:40]+'...'
        plt.text(1.005,.94-(0.06*offset),key+':',     transform=ax1.transAxes,fontsize=8)
        plt.text(1.075,.94-(0.06*offset),metadata_str,transform=ax1.transAxes,fontsize=12)
    #
    #plt.subplots_adjust(left=0.085,right=0.95,top=0.9,bottom=0.05)
    plt.subplots_adjust(left=.04,right=1,top=0.93,bottom=0.07)
    #
    if save_instead_of_plot:
        #filename_txt = 'dashboard_'+metadata['Archive type'].lower()+'_'+str(90-metadata['Latitude ($^\circ$N)'])+'_'+metadata['TSid']+'_'+str(i)
        filename_txt = 'proxy_ts_Temp12k_'+version_txt+'_'+time_span+'_'+metadata['TSid'].replace("'","")
        plt.savefig(output_dir+'figures/'+filename_txt+'.png',dpi=100,format='png')
        plt.close()
    else:
        plt.show()
        #input('Press enter to continue')
    #
    #
    # Save the lats, lons, and TSids
    lats_all.append(filtered_ts[i]['geo_meanLat'])
    lons_all.append(filtered_ts[i]['geo_meanLon'])
    tsid_all.append(filtered_ts[i]['paleoData_TSid'])
    datasetname_all.append(filtered_ts[i]['dataSetName'])
    archive_all.append(filtered_ts[i]['archiveType'])
    proxy_all.append(filtered_ts[i]['paleoData_proxy'])

#%%
# Make versions for printing
lat_string         = ','.join([str('{:.2f}'.format(value)) for value in np.array(lats_all)])
lon_string         = ','.join([str('{:.2f}'.format(value)) for value in np.array(lons_all)])
tsid_string        = ','.join(["'"+value.replace("'","")+"'" for value in np.array(tsid_all)])
datasetname_string = ','.join(["'"+value.replace("'","")+"'" for value in np.array(datasetname_all)])
archive_string     = ','.join(["'"+value.replace("'","")+"'" for value in np.array(archive_all)])
proxy_string       = ','.join(["'"+value.replace("'","")+"'" for value in np.array(proxy_all)])

lat_string         = "      var lat_proxies  = ["+lat_string+"];"
lon_string         = "      var lon_proxies  = ["+lon_string+"];"
tsid_string        = "      var tsid_proxies = ["+tsid_string+"];"
datasetname_string = "      var dataset_proxies = ["+datasetname_string+"];"
archive_string     = "      var archive_proxies = ["+archive_string+"];"
proxy_string       = "      var proxies_proxies = ["+proxy_string+"];"

# Save the latitudes and longitudes to a file
with open(output_dir+'info_temp12k_proxies.txt','w') as f:
    f.write(lat_string+'\n')
    f.write(lon_string+'\n')
    f.write(tsid_string+'\n')
    f.write(datasetname_string+'\n')
    f.write(archive_string+'\n')
    f.write(proxy_string)

print('=== FINISHED ===')

