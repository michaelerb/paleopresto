#==============================================================================
# Make a set of dynamic html files of proxies.
#    author: Michael P. Erb
#    date  : 5/5/2023
#==============================================================================

import numpy as np
import matplotlib.pyplot as plt
import pickle
import lipd
from bokeh.plotting import figure
from bokeh.resources import CDN
from bokeh.embed import file_html
from bokeh.models import HoverTool,Div
from bokeh.layouts import column
from bokeh.models.widgets import Panel,Tabs

plt.style.use('ggplot')

# Choose dataset
version_txt = '1_0_2'


#%% LOAD DATA

# Set the necessary directories
data_dir   = '/projects/pd_lab/data/data_assimilation/proxies/temp12k/'
output_dir = '/projects/pd_lab/mpe32/figures_presto/'

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


#%% Organinze indices by dataset name
n_proxies = len(filtered_ts)
indices_by_dataset = {}
for i in range(n_proxies):
    datasetname = filtered_ts[i]['dataSetName']
    if datasetname not in list(indices_by_dataset.keys()): indices_by_dataset[datasetname] = []
    indices_by_dataset[datasetname].append(i)


#%% For each dataset, put the indices in order of archive type, then 
dataset_names = list(indices_by_dataset.keys())
n_sites = len(dataset_names)
j=4
for j in range(n_sites):
    dataset_name = dataset_names[j]
    indices_selected = indices_by_dataset[dataset_name]
    proxy_and_season = []
    for i in indices_selected:
        proxy_type = filtered_ts[i]['paleoData_proxy']
        season     = filtered_ts[i]['paleoData_interpretation'][0]['seasonalityGeneral']
        proxy_and_season.append(proxy_type+'_'+season)
    #
    ind_sorted = np.argsort(proxy_and_season)
    indices_by_dataset[dataset_name] = np.array(indices_selected)[ind_sorted]
    

#%% Check that the lat and lons are the same for each dataset name
dataset_names = list(indices_by_dataset.keys())
n_sites = len(dataset_names)
for j in range(n_sites):
    dataset_name = dataset_names[j]
    indices_selected = indices_by_dataset[dataset_name]
    for counter,i in enumerate(indices_selected):
        if counter == 0:
            lat_0 = filtered_ts[i]['geo_meanLat']
            lon_0 = filtered_ts[i]['geo_meanLon']
        else:
            if filtered_ts[i]['geo_meanLat'] != lat_0: print('Different lat!',j,i)
            if filtered_ts[i]['geo_meanLon'] != lon_0: print('Different lon!',j,i)


#%% TIME SERIES

# Set up arrays for lat/lons
lats_all = np.zeros(n_sites); lats_all[:] = np.nan
lons_all = np.zeros(n_sites); lons_all[:] = np.nan
latlon_all = []
counters = {}

# Put all proxy records with the same dataset name together and make a figure at every location
j=0
for j in range(n_sites):
    dataset_name = dataset_names[j]
    indices_selected = indices_by_dataset[dataset_name]
    #
    # Get lat and lons
    ind_first = indices_selected[0]
    lats_all[j] = filtered_ts[ind_first]['geo_meanLat']
    lons_all[j] = filtered_ts[ind_first]['geo_meanLon']
    latlon = str(lats_all[j])+'_'+str(lons_all[j])
    if latlon not in latlon_all:
        counters[latlon] = 1
    else:
        # If the lat/lon already appears in the list, shift the new point east by 0.1 degree
        lons_all[j] = lons_all[j] + counters[latlon]*0.1
        counters[latlon] += 1
    #
    latlon_all.append(latlon)
    #
    #print(j,indices_selected)
    tab = {}
    tab_list = []
    for i in indices_selected:
        #
        # Get data and metadata
        proxy_ages = np.array(filtered_ts[i]['age']).astype(float)
        proxy_data = np.array(filtered_ts[i]['paleoData_values']).astype(float)
        #
        # Get proxy metadata
        metadata = {}
        try:    metadata['Dataset name']   = filtered_ts[i]['dataSetName']
        except: metadata['Dataset name']   = ''
        try:    metadata['Archive type']   = filtered_ts[i]['archiveType']
        except: metadata['Archive type']   = ''
        try:    metadata['Proxy type']     = filtered_ts[i]['paleoData_proxy']
        except: metadata['Proxy type']     = ''
        try:    metadata['TSid']           = filtered_ts[i]['paleoData_TSid']
        except: metadata['TSid']           = ''
        try:    metadata['data_url']       = filtered_ts[i]['originalDataUrl']
        except: metadata['data_url']       = ''
        try:    metadata['paper_title']    = filtered_ts[i]['pub1_title']
        except: metadata['paper_title']    = ''
        try:    metadata['seasonality']    = filtered_ts[i]['paleoData_interpretation'][0]['seasonalityGeneral']
        except: metadata['seasonality']    = ''
        try:    data_units = filtered_ts[i]['paleoData_units']
        except: data_units = ''
        try:    age_units  = filtered_ts[i]['ageUnits']
        except: age_units  = ''
        if len(metadata['TSid']) <= 40: tsid_str = metadata['TSid']
        else:                           tsid_str = metadata['TSid'][0:40]+'...'
        lipdverse_url = 'https://lipdverse.org/Temp12k/'+version_txt+'/'+metadata['Dataset name']+'.html'
        html_link = '<a style="text-align: right;" target="_blank" rel="noopener noreferrer" href="'+lipdverse_url+'">Open on LiPDverse</a>'
        #
        title_txt = 'Proxy data: '+metadata['Archive type']+'  |  '+metadata['Proxy type']+'  |  '+metadata['seasonality']+'  |  '+metadata['Dataset name']+'  |  '+tsid_str
        #
        # Make an interactive time series with bokeh
        p1 = figure(plot_width=1200,
                    plot_height=300,
                    title=title_txt,
                    tools='pan,box_zoom,hover,save,reset',
                    active_drag='box_zoom',active_inspect='hover')
        #
        p1.title.text_color = 'purple'
        p1.xaxis.axis_label = 'Age ('+age_units+')'
        p1.yaxis.axis_label = data_units
        p1.x_range.start = 12000
        p1.x_range.end   = 0
        #p1.y_range.start = ts_yrange[0]
        #p1.y_range.end   = ts_yrange[1]
        #
        p1.scatter(proxy_ages,proxy_data,color='purple',size=10,marker='dot')
        p1.line(proxy_ages,proxy_data,color='purple',line_width=1)
        p1.background_fill_color           = 'white'
        p1.grid.grid_line_color            = '#e0e0e0'
        p1.axis.axis_label_text_font_style = 'normal'
        p1.axis.axis_label_text_font_size  = '16px'
        p1.title.text_font_size            = '16px'
        p1.title.align                     = 'center'
        #
        hover = p1.select_one(HoverTool)
        hover.tooltips = [
                #('Age','@x{int} '+age_units),
                ('Age','@x '+age_units),
                ('Data','@y '+data_units),
                ]
        hover.mode='vline'
        #
        season = metadata['seasonality']
        tab[season] = Panel(child=p1,title=metadata['Proxy type']+', '+season)
        tab_list.append(tab[season])
    #
    # Put the tabs together
    all_tabs = Tabs(tabs=tab_list)
    #
    # Create a div with a link to the original data
    link_div = Div(text=html_link,width=1200)
    #
    # Save as html
    full_plot = column(all_tabs,link_div)
    outputfile_txt = 'proxy_temp12k_v'+version_txt+'_'+str(j).zfill(5)+'.html'
    #html = file_html(p1,CDN,outputfile_txt)
    html = file_html(full_plot,CDN,outputfile_txt)
    output_file = open(output_dir+outputfile_txt,'w')
    output_file.write(html)
    output_file.close()


#%% PRINT LAT/LONS FOR WEBSITE

# Make versions for printing
lat_string = ','.join([str('{:.2f}'.format(value)) for value in lats_all])
lon_string = ','.join([str('{:.2f}'.format(value)) for value in lons_all])

lat_string = "              lat_proxies = ["+lat_string+"];"
lon_string = "              lon_proxies = ["+lon_string+"];"

# Save the latitudes and longitudes to a file
with open(output_dir+'latlon_temp12k_proxies.txt','w') as f:
    f.write(lat_string+'\n')
    f.write(lon_string)

print('=== FINISHED ===')

