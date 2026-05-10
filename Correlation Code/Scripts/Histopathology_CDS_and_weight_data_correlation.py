# -*- coding: utf-8 -*-
"""
Created on Wed Aug 27 17:12:23 2025

filename: Histopathology_CDS_and_weight_data_correlation.py

This code is based on data created in EAE_Correlation_Code.py but includes updated code using Pandas DataFrames and more data.

New data includes:
    - weight data for all 32 studied CC strains
    - CDS data for all 32 strains
    - complete B6 control for all data
    - 3 additional CC strains (002, 004, 028) histopathology data
    - An analysis of a binary form of data for periventricular and cuffed perivascular lesions
    
Output goes to a 'Results' folder and contains graphs with pearson and spearman correlations between the various data sets

This code uses the following data sets:
    CDS_Data_All_Strains_8_2025.csv
    Histopathology_Totals_H&E_CC2_4_28_11_18_20_41_43_59_and_B6_8_2025.csv
    Histopathology_Totals_LFB_CC2_4_28_11_18_20_41_43_59_and_B6_8_2025.csv
    Weight_Change_Data_All_Strains_8_2025.csv
    Periventricular-perivascular_cuffing_pathology_H_E_CC11_18_20_41_43_59_B6_10_2025.csv
    Periventricular-perivascular_cuffing_pathology_LFB_CC11_18_20_41_43_59_B6_10_2025.csv
    

@author: Lauren Downs
"""

#%% Modules

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from pathlib import Path as path
import os
import Histopathology_CDS_and_weight_data_correlation_functions_module as fm

#%% Importing Data from Files

# locate input files
cwd = path.cwd()
input_file_path = cwd.parent.joinpath('Input_Data')
result_file_path = cwd.parent.joinpath("Results")

# read input files
CDS_data = pd.read_csv(input_file_path.joinpath('CDS_Data_All_Strains_8_2025.csv'), index_col = 0)
histopath_HE_totals_current_strains = pd.read_csv(input_file_path.joinpath('Histopathology_Totals_H&E_CC2_4_28_11_18_20_41_43_59_and_B6_8_2025.csv'),index_col = 0)
histopath_LFB_totals_current_strains = pd.read_csv(input_file_path.joinpath('Histopathology_Totals_LFB_CC2_4_28_11_18_20_41_43_59_and_B6_8_2025.csv'),index_col = 0)
weight_change_data = pd.read_csv(input_file_path.joinpath('Weight_Change_Data_All_Strains_8_2025.csv'),index_col = 0)
binary_cc43_histopath_HE_current_strains = pd.read_csv(input_file_path.joinpath('Periventricular-perivascular_cuffing_pathology_H_E_CC11_18_20_41_43_59_B6_10_2025.csv'),index_col = 0)
binary_cc43_histopath_LFB_current_strains = pd.read_csv(input_file_path.joinpath('Periventricular-perivascular_cuffing_pathology_LFB_CC11_18_20_41_43_59_B6_10_2025.csv'),index_col = 0)

# making a list for easier sorting 
data_list = [CDS_data,histopath_HE_totals_current_strains,histopath_LFB_totals_current_strains,weight_change_data, binary_cc43_histopath_HE_current_strains, binary_cc43_histopath_LFB_current_strains]
data_list_names = ['CDS',"Histopathology H&E","Histopathology LFB","Weight Change", "Binary Perivascular-Periventricular H&E", "Binary Perivascular-Periventricular LFB"] # labels to follow the indexing

# cleaning blank rows from data
for data_frame in data_list:
    data_frame.dropna(how = 'all', inplace = True)
    data_frame.dropna(how='all', inplace = True, axis = 1)
    
# cleaning up the strain names
for data_frame in data_list:
    data_frame.Strain = data_frame.Strain.str.replace('/','')
    data_frame.Strain = data_frame.Strain.str.strip()
    
# sorting the dataframes row lengths; this way since the IDs are indeces we can use the .index of the smaller data frames to evaluate the larger ones
data_frame_row_len_list = []
for data_frame in data_list:
    data_frame_row_len_list.append(data_frame.shape[0])

row_lengths = pd.DataFrame(data_frame_row_len_list, index=np.arange(len(data_list)))
row_lengths_sort = row_lengths.sort_values(0) # lowest to highest and the index of this will be the who's who of which dataframe gets indexed by the other (smaller indexes larger)

plt.close('all') # closing any previous figures in prep for opening more


#%% Extracting Correlation Data Columns 
# these are the hand selected columns of data that will be plotted and correlated
CDS_corr_data = pd.DataFrame(data_list[0].iloc[:,-4:])
weight_change_corr_data = pd.DataFrame(data_list[3].iloc[:,-3:])
histopath_HE_corr_data = pd.DataFrame(data_list[1].iloc[:,-2:])
histopath_LFB_corr_data = pd.DataFrame(data_list[2].iloc[:,-2:])

corr_data_list = [CDS_corr_data,histopath_HE_corr_data,histopath_LFB_corr_data,weight_change_corr_data]    

# adding a weight loss category that is the negative of average weight change and making Max weight loss positive to match
corr_data_list[3]['Average Weight Loss'] = -1*corr_data_list[3]['Average Weight Change']
corr_data_list[3]['Max Weight Loss'] = -1*corr_data_list[3]['Max Weight Loss']
# Filtering data further for specific correlations

# removing x191
corr_data_no_x191 = []
for data_index, data_frame in enumerate(data_list[0:4]):
    mask = data_frame.EXP != 'X191'
    corr_data_no_x191.append(corr_data_list[data_index][mask])
    



#%% Pearson and Spearman Correlations for all strains across the different data sets

# setting up more organized outputs
all_strain_cor_graph_path = result_file_path.joinpath('All_Strains_Correlation_Graphs')
try: 
    os.mkdir(all_strain_cor_graph_path)
except FileExistsError:
    pass

# for all strains
all_strains_correlation_data = fm.plot_and_correlate_df_list(corr_data_list, all_strain_cor_graph_path, extra_plot_title = 'All Strains' )

# setting up more organized outputs
no_x191_cor_graph_path = result_file_path.joinpath('All_Strains(no_x191)_Correlation_Graphs')
try:
    os.mkdir(no_x191_cor_graph_path)
except FileExistsError:
    pass

# for all strain minus x191 (I think this includes no cc20)
all_strains_no_x191_correlation_data = fm.plot_and_correlate_df_list(corr_data_no_x191, no_x191_cor_graph_path, extra_plot_title = 'All Strains (no x191)' )

#%% Correlation Heat Map

# setting up more organized outputs
correlation_heatmap_path = result_file_path.joinpath('Correlation_Heat_Maps')
try: 
    os.mkdir(correlation_heatmap_path)
except FileExistsError:
    pass
    
# Combing all the data into one dataframe

# getting data columns
corr_data_column_names = []
for dataframe in corr_data_list:
    for column in dataframe.columns:
        corr_data_column_names.append(column)

# creating new combined dataframe - this will contain all the columns from the data sets we are interested in correlating and then as the mouse IDs as the index
correlation_data_combined = pd.DataFrame(columns = corr_data_column_names, index = CDS_data.index)
for dataframe_index, dataframe in enumerate(corr_data_list):
    passed_dataframe, passed_indices = fm.df_index_search(dataframe, CDS_data.index)
    correlation_data_combined.loc[passed_indices,dataframe.columns] = passed_dataframe
    
all_strains_correlation_heat_map_data = fm.map_and_correlate_df(correlation_data_combined,correlation_heatmap_path, extra_plot_title= 'All Strains')

#%% Continued Correlation in individual strains - Heatmaps

# getting all strain names
strain_names = data_list[0].Strain.unique()
strain_dataframe_list = []

# making mini dataframes stored in the list above with all the correlation data for each strain
for each_strain in strain_names:
    strain_ids = CDS_data[CDS_data.Strain == each_strain].index
    strain_dataframe = pd.DataFrame(columns = corr_data_column_names, index = strain_ids)
    for dataframe_index, dataframe in enumerate(corr_data_list):
         passed_dataframe, passed_indices = fm.df_index_search(dataframe, strain_ids)
         strain_dataframe.loc[passed_indices,dataframe.columns] = passed_dataframe
    strain_dataframe_list.append(strain_dataframe)

# heat map with all strains individually 

individual_strain_correlation_heat_map_data = []
for strain_index, strain_data_frame in enumerate(strain_dataframe_list):
    strain_data_frame = strain_data_frame.dropna(axis = 1, how = 'all')
    strain_heat_map_data = fm.map_and_correlate_df(strain_data_frame,correlation_heatmap_path, extra_plot_title=strain_names[strain_index])

#%% Continued Correlation in individual strains - Correlation Graphs
# correlation plots for all strains individually
individual_strain_correlation_data_list = []
for strain_index, each_strain in enumerate(strain_names):
    
    # setting up more organized outputs
    strain_cor_graph_path = all_strain_cor_graph_path.joinpath(f'{each_strain}_Correlation_Graphs')
    try: 
        os.mkdir(strain_cor_graph_path)
    except FileExistsError:
        pass
    
    # breaking up dataframes 
    dissolved_data_frame_list = [strain_dataframe_list[strain_index].iloc[:,1:4], strain_dataframe_list[strain_index].iloc[:,4:8], strain_dataframe_list[strain_index].iloc[:,8:]]
    final_individual_data_frame_list = []
    for index, data_frame in enumerate(dissolved_data_frame_list):
        if data_frame.empty: 
            pass # because pop() wasn't working
        else:
            final_individual_data_frame_list.append(data_frame)

    
    # making plots 
    plot_data = fm.plot_and_correlate_df_list(final_individual_data_frame_list, strain_cor_graph_path, extra_plot_title = each_strain )
    
#%% Assessing Binary Data - Legacy

# THIS AVERAGES BY TISSUE PEICE RATHER THAN TRUE INCIDENCE - DON'T USE THIS
# creating new dataframes to sum the different tissue values and sort them by perivascular vs periventricular
# NOTE: this assumes the data for both H&E and LFB are the same dimesions with the same samples
# no functions used because it is very unlikely that this bare bones scoring system will be used again

# binary_columns = ["HE_SC", "HE_BN","LFB_SC", "LFB_BN"]
# summed_perivascular_data = pd.DataFrame(index=binary_cc43_histopath_HE_current_strains.index, columns = binary_columns )
# summed_periventricular_data = pd.DataFrame(index=binary_cc43_histopath_HE_current_strains.index, columns = binary_columns )

# # looping through and grabbing the H&E data
# for index_entry in binary_cc43_histopath_HE_current_strains.index:
   
#     #spinal cords first
#     perivasc_list = []
#     perivent_list = [] 
    
#     for column in binary_cc43_histopath_HE_current_strains.columns[6:11]: # based on where the SC columns are (iloc will be used rather than names)
#         letter = binary_cc43_histopath_HE_current_strains.loc[index_entry,column] # Becuase the binary data was taken together letter were given for periventricular, perivascular, both or neither
#         if letter == "A":
#             perivasc_list.append(1)
#             perivent_list.append(0)
#         elif letter == "E":
#             perivasc_list.append(0)
#             perivent_list.append(1)
#         elif letter == "B":
#             perivasc_list.append(1)
#             perivent_list.append(1)
#         elif letter != letter: # removes nans
#             pass
#         else:
#             perivasc_list.append(0)
#             perivent_list.append(0)
            
#         # the cummulative score per tissue per stain is the number of tissues that flagged divided by the total number of tissues present for that particular mouse as a way to control for variability in the amount of tissue
#     summed_perivascular_data.loc[index_entry, "HE_SC"] = np.array(perivasc_list).sum() / len(perivasc_list)
#     summed_periventricular_data.loc[index_entry, "HE_SC"] = np.array(perivent_list).sum() / len(perivent_list)
    
#     # brains
#     perivasc_list = []
#     perivent_list = []

#     for column in binary_cc43_histopath_HE_current_strains.columns[11:]: # based on where the SC columns are (iloc will be used rather than names)
#         letter = binary_cc43_histopath_HE_current_strains.loc[index_entry,column] # Becuase the binary data was taken together letter were given for periventricular, perivascular, both or neither
#         if letter == "A":
#             perivasc_list.append(1)
#             perivent_list.append(0)
#         elif letter == "E":
#             perivasc_list.append(0)
#             perivent_list.append(1)
#         elif letter == "B":
#             perivasc_list.append(1)
#             perivent_list.append(1)
#         elif letter != letter:
#             pass
#         else:
#             perivasc_list.append(0)
#             perivent_list.append(0)
            
#         # the cummulative score per tissue per stain is the number of tissues that flagged divided by the total number of tissues present for that particular mouse as a way to control for variability in the amount of tissue
#     summed_perivascular_data.loc[index_entry, "HE_BN"] = np.array(perivasc_list).sum() / len(perivasc_list)
#     summed_periventricular_data.loc[index_entry, "HE_BN"] = np.array(perivent_list).sum() / len(perivent_list)
    
# # looping through and grabbing the LFB data
# for index_entry in binary_cc43_histopath_LFB_current_strains.index:
   
#     #spinal cords first
#     perivasc_list = []
#     perivent_list = [] 
    
#     for column in binary_cc43_histopath_LFB_current_strains.columns[6:11]: # based on where the SC columns are (iloc will be used rather than names)
#         letter = binary_cc43_histopath_LFB_current_strains.loc[index_entry,column] # Becuase the binary data was taken together letter were given for periventricular, perivascular, both or neither
#         if letter == "A":
#             perivasc_list.append(1)
#             perivent_list.append(0)
#         elif letter == "E":
#             perivasc_list.append(0)
#             perivent_list.append(1)
#         elif letter == "B":
#             perivasc_list.append(1)
#             perivent_list.append(1)
#         elif letter != letter:
#             pass
#         else:
#             perivasc_list.append(0)
#             perivent_list.append(0)
            
#         # the cummulative score per tissue per stain is the number of tissues that flagged divided by the total number of tissues present for that particular mouse as a way to control for variability in the amount of tissue
#     summed_perivascular_data.loc[index_entry, "LFB_SC"] = np.array(perivasc_list).sum() / len(perivasc_list) 
#     summed_periventricular_data.loc[index_entry, "LFB_SC"] = np.array(perivent_list).sum() / len(perivent_list)
    
#     # brains
#     perivasc_list = []
#     perivent_list = []

#     for column in binary_cc43_histopath_LFB_current_strains.columns[11:]: # based on where the SC columns are (iloc will be used rather than names)
#         letter = binary_cc43_histopath_LFB_current_strains.loc[index_entry,column] # Becuase the binary data was taken together letter were given for periventricular, perivascular, both or neither
#         if letter == "A":
#             perivasc_list.append(1)
#             perivent_list.append(0)
#         elif letter == "E":
#             perivasc_list.append(0)
#             perivent_list.append(1)
#         elif letter == "B":
#             perivasc_list.append(1)
#             perivent_list.append(1)
#         elif letter != letter:
#             pass
#         else:
#             perivasc_list.append(0)
#             perivent_list.append(0)
            
#         # the cummulative score per tissue per stain is the number of tissues that flagged divided by the total number of tissues present for that particular mouse as a way to control for variability in the amount of tissue
#     summed_perivascular_data.loc[index_entry, "LFB_BN"] = np.array(perivasc_list).sum() / len(perivasc_list)
#     summed_periventricular_data.loc[index_entry, "LFB_BN"] = np.array(perivent_list).sum() / len(perivent_list)
    
# # making a 'Combined' column which sums all the individual ones 
# summed_perivascular_data['Combined'] = summed_perivascular_data.sum(axis=1)
# summed_periventricular_data['Combined'] = summed_periventricular_data.sum(axis=1)

#%% Assessing Binary Data

# Incidence Data Analysis
perivasc_columns = ["HE_SC", "HE_BN","LFB_SC", "LFB_BN"]
perivent_columns = ["HE_BN",'HE_FB',"HE_MB","HE_HB", "LFB_BN", "LFB_FB","LFB_MB","LFB_HB"]
summed_perivascular_data = pd.DataFrame(index=binary_cc43_histopath_HE_current_strains.index, columns = perivasc_columns )
summed_periventricular_data = pd.DataFrame(index=binary_cc43_histopath_HE_current_strains.index, columns = perivent_columns )

# looping through and grabbing the H&E data
for index_entry in binary_cc43_histopath_HE_current_strains.index:
   
    #spinal cords first
    perivasc_list = []
    perivent_list = [] 
    
    for column in binary_cc43_histopath_HE_current_strains.columns[6:11]: # based on where the SC columns are (iloc will be used rather than names)
        letter = binary_cc43_histopath_HE_current_strains.loc[index_entry,column] # Becuase the binary data was taken together letter were given for periventricular, perivascular, both or neither
        if letter == "A":
            perivasc_list.append(1)
            perivent_list.append(0)
        elif letter == "B":
            perivasc_list.append(1)
            perivent_list.append(1)
        elif letter != letter: # removes nans
            pass
        else:
            perivasc_list.append(0)
            perivent_list.append(0)
        # No E because there aren't ventricles in the SC
    
    # if there are any positive flags in the tissue then it gets a 1
    if np.array(perivasc_list).sum() > 0:
        summed_perivascular_data.loc[index_entry, "HE_SC"] = 1
    else:
        summed_perivascular_data.loc[index_entry, "HE_SC"] = 0
    
    
    # brains
    perivasc_list = []
    perivent_list = []

    for column in binary_cc43_histopath_HE_current_strains.columns[11:]: # based on where the SC columns are (iloc will be used rather than names)
        letter = binary_cc43_histopath_HE_current_strains.loc[index_entry,column] # Becuase the binary data was taken together letter were given for periventricular, perivascular, both or neither
        if letter == "A":
            perivasc_list.append(1)
            perivent_list.append(0)
        elif letter == "E":
            perivasc_list.append(0)
            perivent_list.append(1)
        elif letter == "B":
            perivasc_list.append(1)
            perivent_list.append(1)
        elif letter != letter:
            perivent_list.append(np.nan)
        else:
            perivasc_list.append(0)
            perivent_list.append(0)
            
    # if there are any positive flags in the tissue then it gets a 1
    if np.array(perivasc_list).sum() > 0:
        summed_perivascular_data.loc[index_entry, "HE_BN"] = 1
    else:
        summed_perivascular_data.loc[index_entry, "HE_BN"] = 0
        
    if np.array(perivent_list).sum() >0:
        summed_periventricular_data.loc[index_entry, "HE_BN"] = 1
    else:
        summed_periventricular_data.loc[index_entry, "HE_BN"] = 0
   
    # storing the individual brain tissues
    summed_periventricular_data.loc[index_entry, "HE_FB"] = perivent_list[0]
    summed_periventricular_data.loc[index_entry, "HE_MB"] = perivent_list[1]
    summed_periventricular_data.loc[index_entry, "HE_HB"] = perivent_list[2]
        
    
# looping through and grabbing the LFB data
for index_entry in binary_cc43_histopath_LFB_current_strains.index:
   
    #spinal cords first
    perivasc_list = []
    perivent_list = [] 
    
    for column in binary_cc43_histopath_LFB_current_strains.columns[6:11]: # based on where the SC columns are (iloc will be used rather than names)
        letter = binary_cc43_histopath_LFB_current_strains.loc[index_entry,column] # Becuase the binary data was taken together letter were given for periventricular, perivascular, both or neither
        if letter == "A":
            perivasc_list.append(1)
            perivent_list.append(0)
        elif letter == "E":
            perivasc_list.append(0)
            perivent_list.append(1)
        elif letter == "B":
            perivasc_list.append(1)
            perivent_list.append(1)
        elif letter != letter:
            pass
        else:
            perivasc_list.append(0)
            perivent_list.append(0)
            
    # if there are any positive flags in the tissue then it gets a 1
    if np.array(perivasc_list).sum() > 0:
        summed_perivascular_data.loc[index_entry, "LFB_SC"] = 1
    else:
        summed_perivascular_data.loc[index_entry, "LFB_SC"] = 0
    
    # brains
    perivasc_list = []
    perivent_list = []

    for column in binary_cc43_histopath_LFB_current_strains.columns[11:]: # based on where the SC columns are (iloc will be used rather than names)
        letter = binary_cc43_histopath_LFB_current_strains.loc[index_entry,column] # Becuase the binary data was taken together letter were given for periventricular, perivascular, both or neither
        if letter == "A":
            perivasc_list.append(1)
            perivent_list.append(0)
        elif letter == "E":
            perivasc_list.append(0)
            perivent_list.append(1)
        elif letter == "B":
            perivasc_list.append(1)
            perivent_list.append(1)
        elif letter != letter:
            perivent_list.append(np.nan)
        else:
            perivasc_list.append(0)
            perivent_list.append(0)
            
     # if there are any positive flags in the tissue then it gets a 1
    if np.array(perivasc_list).sum() > 0:
        summed_perivascular_data.loc[index_entry, "LFB_BN"] = 1
    else:
        summed_perivascular_data.loc[index_entry, "LFB_BN"] = 0
         
    if np.array(perivent_list).sum() >0:         
        summed_periventricular_data.loc[index_entry, "LFB_BN"] = 1
    else:
        summed_periventricular_data.loc[index_entry, "LFB_BN"] = 0
        
    # storing the individual brain tissues
    summed_periventricular_data.loc[index_entry, "LFB_FB"] = perivent_list[0]
    summed_periventricular_data.loc[index_entry, "LFB_MB"] = perivent_list[1]
    summed_periventricular_data.loc[index_entry, "LFB_HB"] = perivent_list[2]
    
# making combined columns which sums all the individual ones 
summed_perivascular_data['Whole Mouse'] = summed_perivascular_data.sum(axis=1)
summed_periventricular_data['Whole Mouse'] = summed_periventricular_data.sum(axis=1)

summed_perivascular_data['Total_SC'] = summed_perivascular_data[['HE_SC','LFB_SC']].sum(axis=1)
summed_perivascular_data['Total_BN'] = summed_perivascular_data[['HE_BN','LFB_BN']].sum(axis=1)

summed_periventricular_data['Total_FB'] = summed_periventricular_data[['HE_FB','LFB_FB']].sum(axis=1)
summed_periventricular_data['Total_MB'] = summed_periventricular_data[['HE_MB','LFB_MB']].sum(axis=1)
summed_periventricular_data['Total_HB'] = summed_periventricular_data[['HE_HB','LFB_HB']].sum(axis=1)

# converting to binary
summed_perivascular_data['Whole Mouse'][summed_perivascular_data['Whole Mouse'] >0] = 1
summed_periventricular_data['Whole Mouse'][summed_periventricular_data['Whole Mouse'] >0] = 1

summed_perivascular_data['Total_SC'][summed_perivascular_data['Total_SC'] >0] = 1
summed_perivascular_data['Total_BN'][summed_perivascular_data['Total_BN'] >0] = 1

summed_periventricular_data['Total_FB'][summed_periventricular_data['Total_FB'] >0] = 1
summed_periventricular_data['Total_MB'][summed_periventricular_data['Total_MB'] >0] = 1
summed_periventricular_data['Total_HB'][summed_periventricular_data['Total_HB'] >0] = 1

# Plotting the data

# making an index mask by strain
strain_indices_binary = []
strain_labels = np.array(binary_cc43_histopath_HE_current_strains.Strain.unique()) # grabbing strains
strain_labels.sort() # making the listings match the other ordered lists of strains instead of random from the blinding
for strain in strain_labels:
    strain_indices_binary.append(binary_cc43_histopath_HE_current_strains[binary_cc43_histopath_HE_current_strains.Strain == strain].index)

# making results path
binary_data_path = result_file_path.joinpath('Perivascular_Periventricular_Binary_data')
try: 
    os.mkdir(binary_data_path)
except FileExistsError:
    pass

# making bar graphs
perivascular_strain_averages = fm.bar_strain_average(summed_perivascular_data, strain_indices_binary, strain_labels, binary_data_path,title = "Presence of Perivascular Cuffing Lesions all data", x_axis_title= "Strains", y_axis_title="Incidence")
periventricular_strain_averages = fm.bar_strain_average(summed_periventricular_data, strain_indices_binary, strain_labels, binary_data_path,title = "Presence of Periventricular Lesions all data", x_axis_title= "Strains", y_axis_title="Incidence")
# cleaner versions
fm.bar_strain_average(summed_periventricular_data[["HE_BN","LFB_BN","Whole Mouse"]], strain_indices_binary, strain_labels, binary_data_path,title = "Presence of Overall Periventricular Lesions per stain", x_axis_title= "Strains", y_axis_title="Incidence")
fm.bar_strain_average(summed_periventricular_data[["Total_FB","Total_MB","Total_HB","Whole Mouse"]], strain_indices_binary, strain_labels, binary_data_path,title = "Presence of Periventricular Lesions - Tissue Specific", x_axis_title= "Strains", y_axis_title="Incidence")
fm.bar_strain_average(summed_perivascular_data[["Total_SC","Total_BN","Whole Mouse"]], strain_indices_binary, strain_labels, binary_data_path,title = "Presence of Perivascular Cuffing Lesions - Tissue Specific", x_axis_title= "Strains", y_axis_title="Incidence")

# Exporting Data for further Stats testing 

# making new file path
csv_path = result_file_path.joinpath('CSV_files')
try: 
    os.mkdir(csv_path)
except FileExistsError:
    pass

#cleaning up for export
perivascular_path = csv_path.joinpath('perivascular_data.csv')
summed_perivascular_data['Strain'] = binary_cc43_histopath_HE_current_strains['Strain']
summed_perivascular_data['Sex'] = binary_cc43_histopath_HE_current_strains['Sex']
summed_perivascular_data.sort_values('Strain',axis = 0,inplace = True)
summed_perivascular_data.to_csv(perivascular_path)

periventricular_path = csv_path.joinpath('periventricular_data.csv')
summed_periventricular_data['Strain'] = binary_cc43_histopath_HE_current_strains['Strain']
summed_periventricular_data['Sex'] = binary_cc43_histopath_HE_current_strains['Sex']
summed_periventricular_data.sort_values('Strain',axis = 0,inplace = True)
summed_periventricular_data.to_csv(periventricular_path)

#%% Statistics

import Histopathology_CDS_and_weight_data_correlation_functions_module as fm

# Turning the binary data into strain-based outcome data (think contingnecy table for prism)
binary_strain_names = np.array(summed_perivascular_data.Strain.unique())
binary_outcome_list = []
binary_outcome_names = ["Perivascular","Periventricular"]
for data_frame in [summed_perivascular_data,summed_periventricular_data]:
    sub_outcome_list = []
    for strain in binary_strain_names:
        sub_dataframe = data_frame[data_frame.Strain == strain]
        data_columns = sub_dataframe.columns[:-2]
        contingency_data_frame = pd.DataFrame(index = ['yes','no'], columns = data_columns)
        for column in data_columns:
            yes = sub_dataframe[column].sum()
            no = len(sub_dataframe[column]) - yes
            contingency_data_frame[column] = [yes, no]
        sub_outcome_list.append(contingency_data_frame)
    binary_outcome_list.append(sub_outcome_list)

# making new file path
contingency_path = binary_data_path.joinpath('Contingency')
try: 
    os.mkdir(contingency_path)
except FileExistsError:
    pass

# Fisher exact tests between control (B6) and all the other strain
for binary_index, binary_data_set in enumerate(binary_outcome_list):
    contingency_columns = np.array(binary_data_set[0].columns)
    for column in contingency_columns:
        control = binary_data_set[0][column]
        experimental = []
        for strain_index in np.arange(len(binary_strain_names[1:])):
            experimental.append(binary_data_set[strain_index+1][column])
        test = fm.fisher_exact_control_compare(control, experimental, binary_strain_names, contingency_path, extra_plot_title= str(f'for {binary_outcome_names[binary_index]} Lesions in {str(column).replace("_"," ")}'), ylabel = 'Number of Animals', xlabel = 'Strains')
        fm.fisher_exact_control_compare(control, experimental, binary_strain_names, contingency_path, extra_plot_title= str(f'for {binary_outcome_names[binary_index]} Lesions in {str(column).replace("_"," ")} Corrected'), ylabel = 'Number of Animals', xlabel = 'Strains', bon_correction= True)



