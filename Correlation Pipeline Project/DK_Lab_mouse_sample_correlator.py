# -*- coding: utf-8 -*-
"""
Created on Sun Apr 19 14:04:06 2026

filename: DK_Lab_mouse_sample_correlator.py

This pipeline is designed to provide some basic data analysis information, namely correlations, about data collected from
mouse models. Using the provided template (filename), you can add your data to the appropriate columns and have correlation
maps easily produced.

@author: Lauren Downs
"""
#%% Importing Packages

print("Checking for necessary packages...\n")

missing_package_list = []

try:
    import numpy as np
except ImportError:
    missing_package_list.append('Numpy')

try:
    import pandas as pd
except ImportError:
    missing_package_list.append('Pandas')
    
try:
    import seaborn as sns
except ImportError:
    missing_package_list.append('Seaborn')
    
try:
    import scipy.stats as scis
except ImportError:
    missing_package_list.append('SciPy')
    
try:
    from pathlib import Path as path
except ImportError:
    missing_package_list.append('Path')
    
try:
    import os
except ImportError:
    missing_package_list.append('os')
        
try:
    import matplotlib.pyplot as plt
except ImportError:
    missing_package_list.append('Pyplot from Matplotlib')
    
    
if len(missing_package_list) != 0:
    print('The following package(s) need to be installed in order for this script to run:')
    for package in missing_package_list:
        print(f'- {package}')
        
    print('Please make sure all packages are installed and try again.')
else:
    print('All packages successfully imported...\n')
    

#%% Directory and Finding File

# getting the current working directory and setting up a results file path
cwd = path.cwd()
result_file_path = cwd.joinpath("Results")

print(f'Your current working directory is {cwd}\n')

try:
    os.mkdir(result_file_path)
    print('Data from this anallysis will be saved to a folder called "Results" in your current working directory.\n')
except FileExistsError:
    print('Data from this analysis will be saved to the existing "Results" folder.\n')
    
# setting up a sub folders
pairplots_result_file_path = result_file_path.joinpath("Pairplots")
corr_result_file_path = result_file_path.joinpath("Correlation_Heatmaps")
bar_result_file_path = result_file_path.joinpath("Bar_Charts")
box_result_file_path = result_file_path.joinpath("Box_Plots")
data_result_file_path = result_file_path.joinpath("Data")
    
for subfolder in [pairplots_result_file_path,corr_result_file_path,bar_result_file_path, box_result_file_path, data_result_file_path]:
    try:
        os.mkdir(subfolder)
        print(f'Data from this anallysis will be saved to a folder called "{subfolder.name}" in your current working directory.\n')
    except FileExistsError:
        print(f'Data from this analysis will be saved to the existing "{subfolder.name}" folder.\n')

# searching for data file
try:
    data = pd.read_csv('DK_Lab_mouse_sample_correlation_data.csv')
except FileNotFoundError:
    print('"DK_Lab_mouse_sample_correlation_data.csv" not found.')
    bad_file = True
    while bad_file:
        input_filename = input(f'\nPlease enter the name of the .csv file you would like to analyze: ') # checking to see if they saved it as something else
        try:
            data = pd.read_csv(input_filename) # if they include '.csv'
            bad_file = False
            print(f'\n{input_filename} will be used for this analysis')
        except FileNotFoundError:
            try:
                data = pd.read_csv(f'{input_filename}.csv') # if they don't include '.csv'
                bad_file = False
                print(f'\n{input_filename}.csv will be used for this analysis')
            except FileNotFoundError:
                print('File still not found')
    
#%% Functions

def filenamify(name):
    """
    This function takes a string and removes any characters that would not be supported for filenaming. This will also replace spaces with underscores.

    Parameters
    ----------
    name : str
        The input string to be cleaned.

    Returns
    -------
    str
        The output string which can serve as a filename.

    """
    if name.isalnum():
        return name
    else:
        new_name = []
        for character in name:
            if character.isdigit():
                new_name.append(character)
            elif character.isalpha():
                new_name.append(character)
            elif (character == '/') or (character == ':') or (character == '='):
                new_name.append('-')         
            else:
                new_name.append('_')
        return ''.join(new_name)
    
    

#%% Cleaning Data

data.dropna(how = 'all', axis = 1, inplace = True) # getting rid of any columns with no data

columns = np.array(data.columns) # getting this columns


# seeing if the default columns are present
print(f'\ncleaning data...\n')

if np.sort(np.isin(columns,'Strain'))[-1]:
    strain = True
else:
    strain = False
    print('No "Strain" column was found. No analysis based on mouse strain will be performed.\n')

if np.sort(np.isin(columns,'Sex'))[-1]:
    sex = True
else:
    sex = False
    print('No "Sex" column was found. No analysis of sex differences will be performed.\n')
    
if np.sort(np.isin(columns,'EXP'))[-1]:
    exp = True
else:
    exp = False
    print('No "EXP" column was found. No experimental batch effects will be analyzed.\n')
    
# getting some initial data on default columns    

if strain:
    data.Strain = data.Strain.str.replace('/', '')
    data.Strain = data.Strain.str.strip()
    print(f'{len(data.Strain.unique())} unique strain(s) was/were found...\n')

if sex:
    data.Sex = data.Sex.str.strip()
    data.Sex = data.Sex.str.lower()
    print(f'{len(data.Sex.unique())} unique sex(s) was/were found...\n')
    
if exp:
    if data.EXP.dtype != 'int64':
        data.EXP = data.EXP.astype('str')
        data.EXP = data.EXP.str.strip()
    print(f'{len(data.EXP.unique())} unique experiment(s) was/were found...\n')
    
# looking at added columns and sorting them based on whether they are numeric or categorical

added_columns = columns[np.isin(columns, ['Strain','Sex','EXP'], invert = True)]
cat_cols = []
num_cols = []

for col in added_columns:
    # determining if each column's data is numerical....
    col_dtype = data[col].dtype
    if col_dtype == 'float64':
        num_cols.append(col)
    elif col_dtype == 'int64':
        num_cols.append(col)
    elif col_dtype == 'int32':
        num_cols.append(col)
    # ...or categorical 
    else:
        if len(data[col].unique()) > 10: # if there are more than 10 categories (arbitrary number), it will test for is the data is really meant to be numerical
            numeric_datapoints = []
            categorical_datapoints = []
            cat_dp_index = []
            for index, value in enumerate(data[col]):
                try:
                    numeric_datapoints.append(float(value))
                except ValueError:
                    categorical_datapoints.append(value)
                    cat_dp_index.append(index)
            if len(np.unique(numeric_datapoints)) > len(np.unique(categorical_datapoints)):
                num_vs_cat_decision = input(f' For {col} column, {len(np.unique(numeric_datapoints))} unique numeric values were found but only {len(np.unique(categorical_datapoints))} non-numeric value(s) was/were found. Would you like to remove the non-numeric values and treat this column as numeric data? [y] or [n]: ')
                if num_vs_cat_decision.lower() == 'y':
                    data[col] = data[col].drop(cat_dp_index, axis = 0)
                    data[col] = data[col].astype('float')
                    num_cols.append(col)
                else:
                    data[col] = data[col].astype('str')
                    data[col] = data[col].str.strip()
                    cat_cols.append(col)
        else:
            data[col] = data[col].astype('str')
            data[col] = data[col].str.strip()
            cat_cols.append(col)

# adding defaults if they exist
if strain:
    cat_cols = cat_cols + ['Strain']

if sex:
    cat_cols = cat_cols + ['Sex']
    
if exp:
    cat_cols = cat_cols + ['EXP']


# displaying the results
if len(cat_cols) != 0:
    print('\nThe following columns were identified as having categorical data:')
    for col in cat_cols:
        print(f'- {col}')
else:
    print(f'\n No columns were identified as having categorical data.')

if len(num_cols) != 0:
    print('\nThe following columns were identified as having numerical data:')
    for col in num_cols:
        print(f'- {col}')
else:
    print(f'\n No columns were identified as having numerical data.')

#%% Looking at the Data

# Asking for title for plots
if input(f'\nWould you like to add a name for this analysis (applied to all graph titles and file names)? [y] or [n]: ').lower() == 'y':
    extra_graph_title = input(f'Please enter what you would like to add to graph titles: ')
else:
    extra_graph_title = ''
    
print('\nPreparing plots...this may take some time...')

if len(num_cols) != 0:
    
    # looking at the numerical data
    sns.pairplot(data[num_cols],diag_kind = 'hist') #pairplotting
    plot_title = f'{extra_graph_title} Numerical Data'
    plt.suptitle(plot_title, y = 1.0)
    plt.tight_layout()
    plt.savefig(pairplots_result_file_path / f"{filenamify(plot_title)}.png")
    print(f'\nOverall pairplot saved...')
    
    data[num_cols].plot(kind = 'box', subplots = False, figsize = (20,12), rot = 90 ) #box plotting
    plot_title = f'{extra_graph_title} Numerical Data Box Plot'
    plt.title(plot_title)
    plt.savefig(box_result_file_path / f"{filenamify(plot_title)}.png")
    print(f'\nOverall box plot saved...')

# looking at the categorical data

for cat_col in cat_cols:
    plt.figure()
    data[cat_col].value_counts(dropna = False).plot(kind = 'bar')
    plot_title = f'{extra_graph_title} Categorical Data - {cat_col} \n(nan = missing data)'
    plt.title(plot_title)
    plt.savefig(bar_result_file_path / f"{filenamify(plot_title[:-22])}.png", format = 'png')

print(f'\nCategorical bar plots saved...')

#%% Correlations


# correlating numeric data with spearman correlation (all data together)

if len(num_cols) >1:


    spear_corr_coef, spear_corr_p = (scis.spearmanr(data[num_cols].dropna(axis=0))) # doing correlation with testing

    corr_df = pd.DataFrame(spear_corr_coef, # for saving correlation data
                 index = num_cols,
                 columns = num_cols)
    pvalue_df = pd.DataFrame(spear_corr_p, # for saving p-values
                 index = num_cols,
                 columns = num_cols)

    plt.figure(figsize = (15,12)) # plotting
    sns.heatmap(
        corr_df,
        annot= False,
        annot_kws = {'size': 12},
        cmap = 'viridis',
        vmin = -1,
        vmax = 1,
        cbar = False
        )
        
    sns.heatmap( #plotting again for significance
        pd.DataFrame(spear_corr_coef, 
                     index = num_cols,
                     columns = num_cols),
        mask = (spear_corr_p > 0.05), # flipped from regular alpha < # format
        annot= True,
        annot_kws = {'size': 12},
        vmin = -1,
        vmax = 1,
        cmap = 'viridis'
        )
        
    plot_title = f'{extra_graph_title} Spearman Correlations for All Samples \n(shown cofficients = p<0.05)'
    plt.title(plot_title)
    plt.savefig(corr_result_file_path / f"{filenamify(plot_title[:-30])}.png", format = 'png') # saving figure
    corr_df.to_csv(data_result_file_path / f"{filenamify(plot_title[:-30])}_coefficients.csv") # saving data
    pvalue_df.to_csv(data_result_file_path / f"{filenamify(plot_title[:-30])}_pvalues.csv")

    print(f'\nOverall correlation heatmap saved...')
    
    
    
    # correlating numerical data by categorical data
    if len(cat_cols)>0:
        for category in cat_cols:
            
            # showing pairplots
            plot_title = f'{extra_graph_title} Numerical Data by {category}'
            sns.pairplot(data[num_cols + [category]], hue = category, diag_kind = 'hist')            
            plt.suptitle(plot_title, y = 1.00)    
            plt.savefig(pairplots_result_file_path / f"{filenamify(plot_title)}.png")
            print(f'\n{category} pairplot saved...')
            
            # boxplots
            
            # layout formula for box plots
            rows = np.abs((int(len(num_cols)/2))-1)
            if rows == 1:
                cols = 1
            elif len(num_cols)%rows == 0:
                cols = len(num_cols)//rows 
            else:
                cols = (len(num_cols)//rows) +1
                
            data[num_cols + [category]].plot(kind = 'box', by = category, subplots = True, figsize = (20,12), rot = 90, layout = (rows,cols))
            plot_title = f'{extra_graph_title} Numerical Data by {category} Box Plot'
            plt.suptitle(plot_title)
            plt.savefig(box_result_file_path / f"{filenamify(plot_title)}.png")
            print(f'\nOverall box plot saved...')
            
            # correlation heatmaps
            unique_groups = data[category].unique()
            if len(unique_groups) >1:
                for group in unique_groups:
                    
                    #filtering data
                    sub_corr_data = (data[num_cols]
                           .iloc[data.groupby(category)
                           .groups[group]]
                           .dropna(how = 'all', axis = 1)
                           .dropna(axis=0)
                           )
                           
                    sub_corr_data = sub_corr_data.drop(
                        sub_corr_data
                        .nunique()
                        .index[sub_corr_data.nunique() ==1], 
                        axis = 1) # removing cols where values are equal to avoid null correlation results
                    
                           # running and saving correlation
                    spear_corr_coef, spear_corr_p = scis.spearmanr(sub_corr_data)
                                                          
                                                                
                                                     

                    corr_df = pd.DataFrame(spear_corr_coef,
                                 index = sub_corr_data.columns,                                        
                                 columns = sub_corr_data.columns
                                 )
                    pvalue_df = pd.DataFrame(spear_corr_p,
                                 index = sub_corr_data.columns,                                        
                                 columns = sub_corr_data.columns
                                 )

                    plt.figure(figsize = (15,12))
                    sns.heatmap(
                        pd.DataFrame(spear_corr_coef, 
                                     index = sub_corr_data.columns,                                        
                                     columns = sub_corr_data.columns),
                        annot= False,
                        annot_kws = {'size': 12},
                        cmap = 'viridis',
                        vmin = -1,
                        vmax = 1,
                        cbar = False
                        )
                        
                    sns.heatmap(
                        pd.DataFrame(spear_corr_coef, 
                                     index = sub_corr_data.columns,                                        
                                     columns = sub_corr_data.columns),
                        mask = (spear_corr_p > 0.05), # flipped from regular alpha < # format
                        annot= True,
                        annot_kws = {'size': 12},
                        vmin = -1,
                        vmax = 1,
                        cmap = 'viridis'
                        )
                    
                    plot_title = f'{extra_graph_title} Spearman Correlations for {group} \n(shown cofficients = p<0.05)'
                    plt.title(plot_title)
                    plt.savefig(corr_result_file_path / f"{filenamify(plot_title[:-30])}.png", format = 'png') # saving figure
                    corr_df.to_csv(data_result_file_path / f"{filenamify(plot_title[:-30])}_coefficients.csv") # saving data
                    pvalue_df.to_csv(data_result_file_path / f"{filenamify(plot_title[:-30])}_pvalues.csv")
                    print(f'\n{group} correlation heatmap saved...')
                    
                    
else:
    print(f'\nNo correlations made due to no categorical data')


#%% Saving the cleaned data

data.to_csv(data_result_file_path / f"{filenamify(extra_graph_title)}_cleaned_data.csv")
