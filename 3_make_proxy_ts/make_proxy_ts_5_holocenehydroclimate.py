#==============================================================================
# Make a set of dynamic html files of proxies.
#    author: Michael P. Erb
#    date  : 5/5/2023
#==============================================================================

import numpy as np
import matplotlib.pyplot as plt
from bokeh.plotting import figure
from bokeh.resources import CDN
from bokeh.embed import file_html
from bokeh.models import HoverTool
from bokeh.models.widgets import Panel,Tabs

plt.style.use('ggplot')

# Choose dataset
version_txt = '1_0_0'
var_txt = 'hc'
#var_txt = 't'



#%% LOAD METADATA

data_dir     = '/projects/pd_lab/data/proxies/Holocene/HoloceneHydroclimate_github/'
metadata_dir = '/projects/pd_lab/mpe32/HoloceneHydroclimate/Data/Proxy/'
output_dir   = '/projects/pd_lab/mpe32/figures_presto/'

# Load the proxy data. Note: "proxyMetaData_HC_with_pipes.csv" is the same as "proxyMetaData_HC.csv", but commas within the data have been replaced by pipes (|) to make data import easier.
if   var_txt == 'hc': metadata_all = np.loadtxt(metadata_dir+'proxyMetaData_HC_with_pipes.csv',delimiter=',',dtype=str)
elif var_txt == 't':  metadata_all = np.loadtxt(metadata_dir+'proxyMetaData_T.csv',            delimiter=',',dtype=str)
metadata_headers = metadata_all[0,:]
metadata_values  = metadata_all[1:,:]


#%% Organinze indices by dataset name
n_proxies = metadata_values.shape[0]
indices_by_dataset = {}
for i in range(n_proxies):
    datasetname = metadata_values[i,1].replace('"','')
    if datasetname not in list(indices_by_dataset.keys()): indices_by_dataset[datasetname] = []
    indices_by_dataset[datasetname].append(i)


#%% For each dataset, put the indices in order of archive type, then 
dataset_names = list(indices_by_dataset.keys())
n_sites = len(dataset_names)
for j in range(n_sites):
    dataset_name = dataset_names[j]
    indices_selected = indices_by_dataset[dataset_name]
    proxy_and_season = []
    for i in indices_selected:
        proxy_type = metadata_values[i,7].replace('"','')
        season     = metadata_values[i,15].replace('"','')
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
            lat_0 = metadata_values[i,4].astype(float)
            lon_0 = metadata_values[i,3].astype(float)
        else:
            if metadata_values[i,4].astype(float) != lat_0: print('Different lat!',j,i)
            if metadata_values[i,3].astype(float) != lon_0: print('Different lon!',j,i)


#%% TIME SERIES

# Set up arrays for lat/lons
lats_all = np.zeros(n_sites); lats_all[:] = np.nan
lons_all = np.zeros(n_sites); lons_all[:] = np.nan
latlon_all = []
counters = {}

# Put all proxy records with the same dataset name together and make a figure at every location
j,i = 0,0
for j in range(n_sites):
    dataset_name = dataset_names[j]
    indices_selected = indices_by_dataset[dataset_name]
    #
    # Get lat and lons
    ind_first = indices_selected[0]
    lats_all[j] = metadata_values[ind_first,4].astype(float)
    lons_all[j] = metadata_values[ind_first,3].astype(float)
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
        # Get proxy metadata
        metadata = {}
        metadata['Dataset name'] = metadata_values[i,1].replace('"','')
        metadata['Archive type'] = metadata_values[i,6].replace('"','')
        metadata['Proxy type']   = metadata_values[i,7].replace('"','')
        metadata['TSid']         = metadata_values[i,2].replace('"','')
        metadata['seasonality']  = metadata_values[i,15].replace('"','')
        interp                   = metadata_values[i,18].replace('"','')
        if len(metadata['TSid']) <= 40: tsid_str = metadata['TSid']
        else:                           tsid_str = metadata['TSid'][0:40]+'...'
        lipdverse_url = ''
        html_link = '<a style="text-align: right;" target="_blank" rel="noopener noreferrer" href="'+lipdverse_url+'">Open on LiPDverse</a>'
        #
        title_txt = 'Proxy data: '+metadata['Archive type']+'  |  '+metadata['Proxy type']+'  |  '+metadata['seasonality']+'  |  '+metadata['Dataset name']+'  |  '+tsid_str
        #
        # Get data and metadata
        proxy_data_all = np.loadtxt(data_dir+'proxy_'+metadata['TSid']+'.csv',delimiter=',',dtype=str)
        proxy_ages = proxy_data_all[1:,0].astype(float)
        proxy_data = proxy_data_all[1:,1].astype(float)
        age_units  = proxy_data_all[0,0]
        data_units = proxy_data_all[0,1]
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
        if interp == 'negative': p1.y_range.flipped = True
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
    #link_div = Div(text=html_link,width=1200)
    #full_plot = column(all_tabs,link_div)
    #
    # Save as html
    outputfile_txt = 'proxy_holocenehydroclimate_'+var_txt+'_v'+version_txt+'_'+str(j).zfill(5)+'.html'
    #html = file_html(p1,CDN,outputfile_txt)
    html = file_html(all_tabs,CDN,outputfile_txt)
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
with open(output_dir+'latlon_holocenehydroclimate_'+var_txt+'_v'+version_txt+'_proxies.txt','w') as f:
    f.write(lat_string+'\n')
    f.write(lon_string)

print('=== FINISHED ===')

