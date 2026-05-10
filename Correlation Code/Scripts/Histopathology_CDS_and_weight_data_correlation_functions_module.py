# -*- coding: utf-8 -*-
"""
Created on Wed Sep 17 10:44:28 2025

filename: Histopathology_CDS_and_weight_data_correlation_functions_module.py

This file contains functions to be used for Histopathology_CDS_and_weight_data_correlation.py

@author: Lauren Downs
"""

#%% Modules

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import pandas as pd
from pathlib import Path as path
import pandas as pd
import scipy.stats as scis

#%% Functions from Matplotlib

def heatmap(data, row_labels, col_labels, ax=None,
            cbar_kw=None, cbarlabel="", **kwargs):
    """
    Create a heatmap from a numpy array and two lists of labels.

    Parameters
    ----------
    data
        A 2D numpy array of shape (M, N).
    row_labels
        A list or array of length M with the labels for the rows.
    col_labels
        A list or array of length N with the labels for the columns.
    ax
        A `matplotlib.axes.Axes` instance to which the heatmap is plotted.  If
        not provided, use current Axes or create a new one.  Optional.
    cbar_kw
        A dictionary with arguments to `matplotlib.Figure.colorbar`.  Optional.
    cbarlabel
        The label for the colorbar.  Optional.
    **kwargs
        All other arguments are forwarded to `imshow`.
    """

    if ax is None:
        ax = plt.gca()

    if cbar_kw is None:
        cbar_kw = {}

    # Plot the heatmap
    im = ax.imshow(data, **kwargs)

    # Create colorbar
    cbar = ax.figure.colorbar(im, ax=ax, **cbar_kw)
    cbar.ax.set_ylabel(cbarlabel, rotation=-90, va="bottom")

    # Show all ticks and label them with the respective list entries.
    ax.set_xticks(range(data.shape[1]), labels=col_labels,
                  rotation=-40, ha="right", rotation_mode="anchor")
    ax.set_yticks(range(data.shape[0]), labels=row_labels)

    # Let the horizontal axes labeling appear on top.
    ax.tick_params(top=True, bottom=False,
                   labeltop=True, labelbottom=False)

    # Turn spines off and create white grid.
    ax.spines[:].set_visible(False)

    ax.set_xticks(np.arange(data.shape[1]+1)-.5, minor=True)
    ax.set_yticks(np.arange(data.shape[0]+1)-.5, minor=True)
    ax.grid(which="minor", color="w", linestyle='-', linewidth=3)
    ax.tick_params(which="minor", bottom=False, left=False)

    return im, cbar


def annotate_heatmap(im, data=None, valfmt="{x:.2f}",
                     textcolors=("white", "black"),
                     threshold=None, **textkw):
    """
    A function to annotate a heatmap.

    Parameters
    ----------
    im
        The AxesImage to be labeled.
    data
        Data used to annotate.  If None, the image's data is used.  Optional.
    valfmt
        The format of the annotations inside the heatmap.  This should either
        use the string format method, e.g. "$ {x:.2f}", or be a
        `matplotlib.ticker.Formatter`.  Optional.
    textcolors
        A pair of colors.  The first is used for values below a threshold,
        the second for those above.  Optional.
    threshold
        Value in data units according to which the colors from textcolors are
        applied.  If None (the default) uses the middle of the colormap as
        separation.  Optional.
    **kwargs
        All other arguments are forwarded to each call to `text` used to create
        the text labels.
    """

    if not isinstance(data, (list, np.ndarray)):
        data = im.get_array()

    # Normalize the threshold to the images color range.
    if threshold is not None:
        threshold = im.norm(threshold)
    else:
        threshold = im.norm(data.max())/2.

    # Set default alignment to center, but allow it to be
    # overwritten by textkw.
    kw = dict(horizontalalignment="center",
              verticalalignment="center")
    kw.update(textkw)

    # Get the formatter in case a string is supplied
    if isinstance(valfmt, str):
        valfmt = mpl.ticker.StrMethodFormatter(valfmt)

    # Loop over the data and create a `Text` for each "pixel".
    # Change the text's color depending on the data.
    texts = []
    for i in range(data.shape[0]):
        for j in range(data.shape[1]):
            if type(data[i, j]) == str:
                pass
            else:
                kw.update(color=textcolors[int(im.norm(data[i, j]) > threshold)])
                text = im.axes.text(j, i, valfmt(data[i, j], None), **kw)
                texts.append(text)

    return texts



#%% Functions

def plot_and_correlate_df_list(data_frame_list, result_file_path, extra_plot_title = None, significance = 0.05, both_sig_color = (1,.412,.161,1), pear_sig_color = (.129,.859,1,1),spear_sig_color=(.871,.412,1,1), stat_color_key='\nTeal = Pearson significance,\nPink = Spearman,\nOrange = both significance' ):
    """
    This function generates all combinations of correlations between the columns of the different data frames from a list and scatter plots them. Note that the indeces of the rows for the dataframes must be as equivalent as possible. Smaller dataframes should go first and larger ones later in the indexing. If any indeces do not match after this they will be skipped.

    Parameters
    ----------
    data_frame_list : list of pandas DataFrames
        This is a list of dataframes where the columns are the labels of the data and each column needs to be correlated with the columns of the other dataframes
    result_file_path : str or Path object
        This is the name and file path where a .png of the figures will be saved. This will be combined with the title of the plot.
    extra_plot_title : str, optional
        This allows for another title to be added in front of the plot title which is the first column name 'x' second column name. The default is None.
    significance : float, optional
        This is the value the p-value must be less than to have it be colored with the significance color. The default is 0.05.
    both_sig_color : tuple (1x4), optional
        RBGA color for the text when both correlations are significant. The default is (1,.412,.161,1).
    pear_sig_color : tuple (1x4), optional
        RBGA color for the text when the pearson correlation is significant. The default is (.129,.859,1,1).
    spear_sig_color : tuple (1x4), optional
        RBGA color for the text when the spearman correlation is significant. The default is (.871,.412,1,1).
    stat_color_key : str, optional
        Writing describing what the colored text means. The default is '\nTeal = Pearson significance,\nPink = Spearman,\nOrange = both significance'.

    Returns
    -------
    correlation_data: pandas Dataframe [data_frame_list column pairs x 5]
        A dataframe containing the correlation coefficients and p values for both pearson and spearman as well as what indices were used to do the correlation; the indies are the correlated column names separated by " vs "

    """
    
    # set up
    correlation_data_columns = ['Pearson Coefficient', 'Pearson p-value', 'Spearman Coefficient', 'Spearman p-value', 'Indeces Used']
    correlation_data = pd.DataFrame(columns = correlation_data_columns)
    if len(plt.get_fignums()) != 0: # for not rewritting plots that are open
        figure_number = np.array(plt.get_fignums()).max() +1
    else: 
        figure_number = 1
    
    # filtering nans
    for index, data_frame in enumerate(data_frame_list):
        data_frame.dropna(how = 'any', inplace = True)
        if data_frame.empty:
            data_frame_list.pop(index)
        
    # sorting the dataframes row lengths; this way since the IDs are indeces we can use the .index of the smaller data frames to evaluate the larger ones
    data_frame_row_len_list = []
    for data_frame in data_frame_list:
        data_frame_row_len_list.append(data_frame.shape[0])

    row_lengths = pd.DataFrame(data_frame_row_len_list, index=np.arange(len(data_frame_list)))
    row_lengths_sort = row_lengths.sort_values(0) # lowest to highest and the index of this will be the who's who of which dataframe gets indexed by the other (smaller indexes larger)
    list_index = row_lengths_sort.index

    for row_lengths_index, data_frame_index in enumerate(list_index):
        for corr_row_indexes in np.arange(row_lengths_index+1,len(list_index)):
            for i in np.arange(data_frame_list[data_frame_index].shape[1]): # i cycles through the first dataframe's columns
                for j in np.arange(data_frame_list[list_index[corr_row_indexes]].shape[1]):# j cycles through the correlating dataframe's columns
                    
                    # setting up collection for Indeces Used
                    indeces_used = []
                    
                    # these can be used for the cycles through the columns through the data frames (these are pd series)
                    small_temp_series = pd.Series(data_frame_list[data_frame_index].iloc[:,i], dtype = float) # not necessarilly smaller but would be the smaller one if there is a difference
                    
                    # this is to account for any missing entries in the larger data set due to nan filtering 
                    large_temp_series = pd.Series(name = data_frame_list[list_index[corr_row_indexes]].iloc[:,j].name, dtype= float) # it apparently needs the name both here and below to actually show up in the saved data
                    for entry in data_frame_list[data_frame_index].index:
                        try: # the following is unruly because Pandas decided to remove the append attribute :(
                            large_temp_series = pd.concat([large_temp_series, pd.Series(data_frame_list[list_index[corr_row_indexes]].iloc[:,j].loc[entry], index = np.array([entry]), name = data_frame_list[list_index[corr_row_indexes]].iloc[:,j].name)])
                            indeces_used.append(entry)
                        except KeyError:
                            small_temp_series.pop(entry)
                    
                    # doing the correlations                
                    pearson = scis.pearsonr(small_temp_series, large_temp_series)
                    spearman = scis.spearmanr(small_temp_series, large_temp_series)
                    correlation_name = str((f"{extra_plot_title} {small_temp_series.name} vs {large_temp_series.name}"))
                                   
                    # saving the data
                    correlation_data = pd.concat([correlation_data,pd.DataFrame(data = np.array([[pearson[0], pearson[1], spearman[0], spearman[1], str(indeces_used)]]), columns = correlation_data_columns, index = np.array([correlation_name]))])
                                    
                    # plotting
                    plt.figure(figure_number, figsize=(7,5))
                    plt.scatter(small_temp_series, large_temp_series)
                    
                    # plot labels
                    plt.title(correlation_name)
                    plt.xlabel(small_temp_series.name)
                    plt.ylabel(large_temp_series.name)
                    
                    # correlation color coding
                    if (float(pearson[1]) < significance) and (float(spearman[1]) <significance):
                        text_color = both_sig_color
                    elif pearson[1] <significance:
                        text_color=pear_sig_color
                    elif spearman[1] <significance:
                        text_color=spear_sig_color
                    else:
                        text_color = None
                    
                    # adding the correlations
                    plt.figtext(.95,.5,f'Pearson: {pearson[0]:.3f}\np-value: {pearson[1]:.3f}\nSpearman: {spearman[0]:.3f}\np-value: {spearman[1]:.3f}', color = text_color)
                    plt.figtext(.95,.5-0.1, stat_color_key, size = 7)
                    
                    
                    # saving a preparing for next graph
                    filename = result_file_path / f"{str(correlation_name).replace(' ','_')}_correlation.png"
                    plt.savefig(filename, format = 'png',bbox_inches = 'tight')
                    figure_number += 1
                    
    return correlation_data

def df_index_search(dataframe, index_list):
    """
    Searches a pandas data frame for any indices in the given index list; returns np.nans if not found.

    Parameters
    ----------
    dataframe : pandas DataFrame
        A dataframe with indices matching some, none, or all of the ones in the index list.
    index_list : list or similar interative object
        Contains indices to check for in the dataframe.

    Returns
    -------
    return_dataframe : pandas DataFrame
        A dataframe with all values from the matching indices.
    index_pass_list : list
        All the indices that returned a value from the given dataframe.

    """
    
    index_pass_list = []
    return_dataframe = pd.DataFrame(columns = dataframe.columns, index = index_list)
    for index in index_list:
        try:
            return_dataframe.loc[index] = dataframe.loc[index]
            index_pass_list.append(index)
        except KeyError:
            pass
    return return_dataframe, index_pass_list

def map_and_correlate_df(data_frame, result_file_path, extra_plot_title = None, significance = 0.05, font_size = 8, font_weight = 'bold', figure_size = (7,7)):
    """
    This function generates all combinations of correlations between the columns of the different data frames from a list and scatter plots them. Note that the indeces of the rows for the dataframes must be as equivalent as possible. Smaller dataframes should go first and larger ones later in the indexing. If any indeces do not match after this they will be skipped.

    Parameters
    ----------
    data_frame : pandas DataFrame
        TBD
    result_file_path : str or Path object
        This is the name and file path where a .png of the figures will be saved. This will be combined with the title of the plot.
    extra_plot_title : str, optional
        This allows for another title to be added in front of the plot title which is the first column name 'x' second column name. The default is None.
    significance : float, optional
        This is the value the p-value must be less than to have it be colored with the significance color. The default is 0.05.
    
    Returns
    -------
    correlation_data: pandas Dataframe [data_frame_list column pairs x 5]
        A dataframe containing the correlation coefficients and p values for both pearson and spearman as well as what indices were used to do the correlation; the indies are the correlated column names separated by " vs "

    """
    
    # set up
    pear_coefficients = pd.DataFrame(index=data_frame.columns, columns=data_frame.columns)
    spear_coefficients = pd.DataFrame(index=data_frame.columns, columns=data_frame.columns)
    pear_pvalues = pd.DataFrame(index=data_frame.columns, columns=data_frame.columns)
    spear_pvalues = pd.DataFrame(index=data_frame.columns, columns=data_frame.columns)
    if len(plt.get_fignums()) != 0: # for not rewritting plots that are open
        figure_number = np.array(plt.get_fignums()).max() +1
    else: 
        figure_number = 1
    
    # filtering nans
    data_frame.dropna(how = 'any', inplace = True)

    for col1 in data_frame.columns:
        for col2 in data_frame.columns:
            
            # doing the correlations
            pearson = scis.pearsonr(np.array(data_frame[col1], dtype = float), np.array(data_frame[col2], dtype = float))
            spearman = scis.spearmanr(np.array(data_frame[col1], dtype = float), np.array(data_frame[col2], dtype = float))
                                                       
            # saving the data
            pear_coefficients.loc[col1,col2] = pearson[0]
            pear_pvalues.loc[col1,col2] = pearson[1]
            spear_coefficients.loc[col1,col2] = spearman[0]
            spear_pvalues.loc[col1,col2] = spearman[1]
            
    # plotting and making figure
    
    # pearson
    plt.figure(figure_number, figsize = figure_size)
    im, cbar = heatmap(np.array(pear_coefficients, dtype = float),np.array(pear_coefficients.index),np.array(pear_coefficients.columns))
    significant_pearson = np.array(pear_coefficients[pear_pvalues <0.05], dtype = 'object')
    np.putmask(significant_pearson, pd.isna(significant_pearson), 'ns')
    
    annotate_heatmap(im, np.array(significant_pearson), threshold = 0.05, fontsize = font_size, fontweight = font_weight)
    if extra_plot_title != None:
        plot_title = f'{extra_plot_title} Pearson Correlation'
        plt.title(plot_title)
    else:
        plot_title ='Pearson Correlation'
        plt.title(plot_title)
        
    # saving a preparing for next graph
    filename = result_file_path / f"{str(plot_title).replace(' ','_')}_correlation_heat_map.png"
    plt.savefig(filename, format = 'png',bbox_inches = 'tight')
    figure_number += 1
                    
    # spearman
    plt.figure(figure_number, figsize=(7,7))
    im, cbar = heatmap(np.array(spear_coefficients, dtype = float),np.array(spear_coefficients.index),np.array(spear_coefficients.columns))
    significant_spearman = np.array(spear_coefficients[spear_pvalues<0.05], dtype = 'object')
    np.putmask(significant_spearman, pd.isna(significant_spearman),'ns')
    
    annotate_heatmap(im, np.array(significant_spearman), threshold = 0.05, fontsize = font_size, fontweight = font_weight)
    if extra_plot_title != None:
        plot_title = f'{extra_plot_title} Spearman Correlation'
        plt.title(plot_title)
    else:
        plot_title ='Spearman Correlation'
        plt.title(plot_title)
        
    # saving a preparing for next graph
    filename = result_file_path / f"{str(plot_title).replace(' ','_')}_correlation_heat_map.png"
    plt.savefig(filename, format = 'png',bbox_inches = 'tight')
    figure_number += 1
    
    return pear_coefficients, pear_pvalues, significant_pearson, spear_coefficients, spear_pvalues, significant_spearman


def bar_strain_average(data_frame, master_strain_index_list, strain_lables, result_file_path, title = None, x_axis_title = None, y_axis_title = None):
    
    # set up for plotting
    if len(plt.get_fignums()) != 0: # for not rewritting plots that are open
        figure_number = np.array(plt.get_fignums()).max() +1
    else: 
        figure_number = 1
    
    # setting up x axis
    x = np.arange(len(strain_lables))
    
    # getting strain averages of the data for the heights of the bars
    
    # for one data set 
    if len(data_frame.shape) ==1:
        heights = []
        for strain_index in master_strain_index_list:
            heights.append(data_frame[strain_index].sum()/len(strain_index))
    
        # plotting for one data set
        plt.figure(figure_number)
        plt.bar(x,heights)
        plot_title = f"Strain Averages {title}"
        plt.title(plot_title)
        plt.ylabel(y_axis_title)
        plt.xlabel(x_axis_title)
        plt.xticks(x, labels = strain_lables)
        
        # saving the plot
        filename = result_file_path / f"{str(plot_title).replace(' ','_')}_bar_graph.png"
        plt.savefig(filename, format = 'png',bbox_inches = 'tight')
        
        return heights
        
    else:
        plt.figure(figure_number)
        column_number = len(data_frame.columns)
        x_width = .8/column_number
        x_buffer = x_width*column_number/(0-2)
        total_heights = []
        for column in data_frame:
            heights = []
            for strain_index in master_strain_index_list:
                heights.append(data_frame.loc[strain_index, column].sum()/len(strain_index))
        
            # plotting for per column
            plt.bar(x+x_buffer,heights, alpha = 0.5, width = x_width, align = "edge")
            x_buffer +=x_width
        
            # final plot additions
            plot_title = f"Strain Averages {title}"
            plt.title(plot_title)
            plt.ylabel(y_axis_title)
            plt.xlabel(x_axis_title)
            plt.xticks(x, labels = strain_lables)
            plt.legend(data_frame.columns)
        
            # saving the plot
            filename = result_file_path / f"{str(plot_title).replace(' ','_')}_bar_graph.png"
            plt.savefig(filename, format = 'png',bbox_inches = 'tight')
            total_heights.append(heights)
            figure_number+=1
        
        return total_heights
        
        
        
def fisher_exact_control_compare(control_column, experimental_columns, data_names, result_file_path, extra_plot_title = None, ylabel = None, xlabel = None, significance = 0.05, bon_correction = False):
    """
    This function takes contingency data (specifically formatted as 'yes' and 'no') and does a fisher exact test between a given control and n number of experimental data sets as then plots these as stacked bar graphs.

    Parameters
    ----------
    control_column : Pandas DataFrame (2,1)
        Integer number of total 'yes's and 'no's for the control data
    experimental_columns : Pandas DataFrame (2,n)
        Integer number of total 'yes's and 'no's for experimental data sets
    data_names : array-like or length n+1
        Strings labeling the control and experimental data sets
    result_file_path : str or Path object
        This is the name and file path where a .png of the figures will be saved. This will be combined with the title of the plot.
    extra_plot_title : str, optional
        Passed to the plot title for description of what is being compared. The default is None.
    ylabel : str, optional
        Describes y-axis. The default is None.
    xlabel : str, optional
        Describes x-axis. The default is None.
    significance : float, optional
        The p-value below which bars will be labeled as significant. The default is 0.05.

    Returns
    -------
    fisher_data_list : array-like (length n) of SignificanceResult (length 2)
        Coefficients and p-values for every experimental group as compared using a fisher exact test to the control.

    """
    
    # Doing the statistics
    fisher_data_list = []
    for experiment in experimental_columns:
        fisher_dataframe = pd.DataFrame([control_column,experiment])
        fisher_data_list.append(scis.fisher_exact(fisher_dataframe))
        
    # Chi square
    chi_dataframe = pd.DataFrame(control_column)
    for experiment in experimental_columns:
        data_frame = pd.DataFrame(experiment)
        chi_dataframe = pd.concat([chi_dataframe,experiment], axis = 1)
    chi_value = scis.chisquare(chi_dataframe.iloc[1],chi_dataframe.iloc[0], sum_check= False)
    
    # Plotting Set up
    if len(plt.get_fignums()) != 0: # for not rewritting plots that are open
        figure_number = np.array(plt.get_fignums()).max() +1
    else: 
        figure_number = 1
    plotting_dataframe = pd.DataFrame([control_column])
    plotting_dataframe =pd.concat( [plotting_dataframe, pd.DataFrame(experimental_columns)] )

    
    # Plotting the bar graphs
    plt.figure(figure_number)
    plt.bar(data_names, plotting_dataframe.yes+plotting_dataframe.no, label = 'no')
    plt.bar(data_names, plotting_dataframe.yes, label = 'yes')
    
    # labeling the plot
    plt.legend()
    plt.ylabel(ylabel)
    plt.xlabel(xlabel)
    
    # setting up title
    if extra_plot_title != None:
        plot_title = f'Contingency Comparison {extra_plot_title}'
    else:
        plot_title ='Contingency Comparison'
    saved_plot_title = plot_title
    
    # adding significance markings
    for index, fisher_data in enumerate(fisher_data_list):
        if bon_correction:
            if fisher_data[1] <significance/len(experimental_columns):
                plt.text(data_names[index+1], (plotting_dataframe.yes+plotting_dataframe.no).iloc[index+1],'*', ha='center', size = 90/len(data_names))
        elif fisher_data[1] <significance:
            plt.text(data_names[index+1], (plotting_dataframe.yes+plotting_dataframe.no).iloc[index+1],'*', ha='center', size = 90/len(data_names))
    
    if (chi_value[1] < significance) and (chi_value[1] > 0):
        #plt.text(.5*len(experimental_columns),-3, f'Chi Square: {chi_value[0]:.3f}', ha = 'center')
        plot_title = plot_title + str(f'\nChi Square: {chi_value[0]:.3f}')
    plt.title(plot_title)
    
    # Saving the plot
    filename = result_file_path / f"{str(saved_plot_title).replace(' ','_')}_bar_graph.png"
    plt.savefig(filename, format = 'png',bbox_inches = 'tight')
    figure_number+=1

    return fisher_data_list, chi_value

    