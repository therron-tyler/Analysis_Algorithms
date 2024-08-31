#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import matplotlib
import pandas as pd
import numpy as np
import sys
import csv
import os
from scipy.stats import zscore, variation
import time
import seaborn as sns
import matplotlib.pyplot as plt
import argparse # adding this library so this function can be called from the command line

def Corr_Heatmap_VII_40under(_project_name,_df_path,sample_list_path,out_path):
    # Correlation Heatmap edits

    # read the normalized counts in
    #df_path = _df_path
    df_unordered = pd.read_csv(str(_df_path), sep='\t', header='infer', index_col=0)
    print("DataFrame columns:\n", df_unordered)
    # Read in Sample List to reorganize the columns
    file_path = str(sample_list_path)
    sample_list = pd.read_csv(file_path, header=None)
    ordered_sample_names = sample_list[0].tolist()
    print("Ordered sample names:\n", ordered_sample_names)
    # use list to reorganize the CPM counts df for mapping
    # Reorganize the DataFrame using the list of column names
    df = df_unordered[ordered_sample_names]

    # create correlational matrix using Pearson's coefficient
    corr_df = df.corr(method='pearson')
    # creates a list of column values 
    _cols_cpm = list(corr_df.columns.values)


    tick_pos = np.arange(0.5, len(_cols_cpm)+0.5)


    fig = plt.figure(figsize=(26,18)) # increased the figure size
    ax = fig.add_subplot(111)
    colormap = sns.diverging_palette(220, 13, as_cmap=True)
    sns.heatmap(corr_df, cmap=colormap, linecolor='white', annot=False, linewidths=1, ax=ax, cbar_kws={'shrink': 0.7},annot_kws={'size': 80})
    
    # code to swap underscores with spaces '_' -> ' '
    
    name_underscores = str(_project_name)
    title_spaces = name_underscores.replace('_', ' ')
    plt.title(title_spaces, fontsize=16)

    
    
    ax.set_xticks(tick_pos)
    ax.set_xticklabels(_cols_cpm, ha='right', rotation=25, fontsize=50) # 10 larger x axis ticks and labels remove


    ax.set_yticks(tick_pos)
    ax.set_yticklabels(_cols_cpm, ha='right', va='top',rotation=25, minor=False, fontsize=50, y=tick_pos)


    ax.tick_params(axis='x', length=2, color='black', labelsize=5, pad=2, bottom=True, labelbottom=True)
    ax.tick_params(axis='y', length=2, color='black', labelsize=5, pad=2) 
    
    ax.collections[0].colorbar.set_label("Pearson's Correlation Coefficient",fontsize = 14) 
    ax.collections[0].colorbar.ax.tick_params(labelsize=14)

    fig.savefig(str(out_path)+'/'+name_underscores+'_Pearson_Correlation_Heatmap.png', dpi=300) # inlcude '/' before and after the out_path argument
    plt.close()
    return None
if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description="Correlation Heatmap Argument Parser")
    
    parser.add_argument("-name", "--project-name", required=True, 
                        help="[REQUIRED] Name the Title of Plot and File.")
    
    parser.add_argument("-df", "--dataframe-path", type=os.path.abspath, required=True,
                        help="[REQUIRED] Provide the path to the DataFrame CSV file.\n -df [~/DATAFRAME-CSV-FILE-PATH/],\t--dataframe-path [~/DATAFRAME-CSV-FILE-PATH/]\n")

    parser.add_argument("-sl", "--sample-list-path", type=os.path.abspath, required=True,
                        help="[REQUIRED] Provide the path to the Sample List CSV file.\n -sl [~/SAMPLE-LIST-CSV-FILE-PATH/],\t--sample-list-path [~/SAMPLE-LIST-CSV-FILE-PATH/]\n")

    parser.add_argument("-op", "--output-path", type=os.path.abspath, required=True,
                        help="[REQUIRED] Provide the path to save the output plot.\n -op [~/OUTPUT-PLOT-FILE-PATH/],\t--output-path [~/OUTPUT-PLOT-FILE-PATH/]\n")

    args = parser.parse_args()

    _project_name = args.project_name
    _df_path = args.dataframe_path
    sample_list_path = args.sample_list_path
    out_path = args.output_path

    print(f"Name of Plot : {_project_name}")
    print(f"DataFrame Path : {_df_path}")
    print(f"Sample List Path : {sample_list_path}")
    print(f"Output Path : {out_path}")

    
    Corr_Heatmap_VII_40under(_project_name,_df_path,sample_list_path,out_path)

