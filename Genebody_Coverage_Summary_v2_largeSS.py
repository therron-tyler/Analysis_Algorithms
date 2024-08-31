import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
from sklearn.preprocessing import normalize
import argparse
from sklearn.preprocessing import MinMaxScaler
from matplotlib.backends.backend_pdf import PdfPages
import seaborn

def geneBody_summaryTable(_input_txt_dir, _out_dir, _project_name):
# define the directory path
    my_path = _input_txt_dir
 
    # get the list of files in the directory
    files_list = os.listdir(my_path) 

    # create output directory for figure
   #  - no need to make new directory
    # definer scaler as object
    scaler = MinMaxScaler()   
# create an empty list 
    df_list = []

    df_list2 = []
#empty pandas dataframe to append to
    genebody_df = pd.DataFrame()

# iterate through each file in the directory
    for file in files_list:
    # read the data of each file
        with open(my_path + "/" + file, 'r') as f: # 2/6/23 - added slash because of file location error
            data = pd.read_csv(f, sep="\t")
        # read in the tab delimited data into the 'data' object and make a dataframe from it
            pd.DataFrame(data)
        
        # append data to list
            df_list.append(data)
        
        

    result = pd.concat(df_list, axis=0, join='outer', ignore_index=True)

# Rearranging index
    result.index = np.arange(1, len(result) + 1)

# bam files have the '_chr' cut from the name
    result['Percentile'] = result['Percentile'].str.replace('_chr', '')
    
    # this is where the changes begin with the new normalization methods after talking with Debbie
    result_transform = result
    result_transform = result_transform.drop('Percentile',axis=1)
    


    result_scaled=pd.DataFrame(scaler.fit_transform(result_transform.T).T,columns=result_transform.columns)
    
    # MaxMin scale transformation should occur here - then the scaled results can be loaded into the other dataframe 


    # index and insert column names to normalized data
    result_scaled.index = np.arange(1, len(result) + 1)
    result_scaled.insert(0, 'Percentile', result['Percentile'])
    result_scaled = result_scaled.sort_values("Percentile")

# 2/15/2023 - insert chunk for batch genebody coverage line plot ------------------------------------

    # transpose the columns so they can be group as an x-axis and y-axis   
    result_scaled_noindex = result_scaled.reset_index(drop=True)
    result_scaled_noindex.index = np.arange(1, len(result_scaled_noindex) + 1)
    result_scaled_transposed = result_scaled_noindex.transpose()
    result_scaled_transposed = result_scaled_transposed.reset_index()

# Use the rename() method to make the first row the column headers:
    result_scaled_transposed = result_scaled_transposed.rename(columns=result_scaled_transposed.iloc[0])
    result_scaled_transposed = result_scaled_transposed.drop(result_scaled_transposed.index[0])

# make everything in the dataframe numeric type
    result_scaled_transposed = result_scaled_transposed.apply(pd.to_numeric)

# convert to long (tidy) form
    dfm = result_scaled_transposed.melt('Percentile', var_name='Samples', value_name='Normalized Counts')

# percentile is a column now - so it can be uzed as an x-variable on the plot
# plot size
    plt.figure(figsize=(20,16))

    #colormap_L = sns.cubehelix_palette(start=1, rot=.5, gamma=1.1,hue=1, as_cmap=False,light=1, n_colors=68)
    print("about to make line plot")
    #palette = sns.color_palette("husl", 83)
    # instead of hardcoded number of samples - dynamically count them for palette object creation

    num_samples = dfm["Samples"].nunique()
    palette = sns.color_palette("husl", num_samples)
    seaborn.lineplot(x = "Percentile", y = "Normalized Counts", data=dfm, hue = "Samples", palette = palette)
    

# Move the legend to the upper left
    ax = plt.gca()
    ax.legend(loc='upper left', fontsize = "xx-small")
    ax.tick_params(axis='both', which='major', labelsize=20)
    ax.set_xlabel("Genebody Percentile (5' -> 3')",fontsize = 20)
    ax.set_ylabel("Normalized Counts",fontsize = 20)
    ax.set_title(str(_project_name)+" Genebody Coverage",fontsize = 28)
# save the figure to a pdf
    plt.savefig(str(_out_dir)+"/"+str(_project_name)+"_GeneBody_Coverage_LinePlot.pdf", format='pdf')

# write datatable to the csv file ---------------------------------------------------
    print("about to make CSV")
    result_scaled.to_csv(str(_out_dir)+"/"+str(_project_name)+"_Genebody_Coverage_summary.csv", index=True) # the output is the specified out-directory and desired filename
#  make the heatmap for genebody coverage ---------------------------------------------
    result_scaled = result_scaled.set_index(result_scaled.columns[0])

    fig = plt.figure(figsize=(35,44))
    ax = fig.add_subplot(111)
    colormap = sns.diverging_palette(220, 13, as_cmap=True)
    print("about to make heatmap")
#colormap = sns.cubehelix_palette(start=1, rot=1, gamma=0.8, as_cmap=True, n_colors=80)
    sns.heatmap(result_scaled, cmap=colormap, annot=False, ax=ax, cbar_kws={'shrink': 0.7},square=False)
    plt.title('Percentile Genebody Coverage Heatmap for '+str(_project_name), fontsize=28, pad=15)

    ax.tick_params(axis='x', length=8, color='black', labelsize=14, pad=5, bottom=True, labelbottom=True,rotation=70)
    ax.tick_params(axis='y', length=8, color='black', labelsize=9, pad=10, rotation=0)
    ax.set_xlabel("Genebody Percentile (5' -> 3')",fontsize = 28)
    ax.set_ylabel("Sample",fontsize = 28)
    ax.collections[0].colorbar.set_label("Min-Max Normalized Counts",fontsize = 24)
    ax.collections[0].colorbar.ax.tick_params(labelsize=20)
    fig.savefig(str(_out_dir)+"/"+str(_project_name)+"_Genebody_Coverage_Heatmap.pdf", bbox_inches='tight', format='pdf',dpi=300)

    
    return None

if __name__ == "__main__":
### Run Genebody_Coverage_summary_v2.py() to take input using commandline arguments (USER Input)   
    
    parser = argparse.ArgumentParser(description="Genebody Summary script for tab-demlimited un-normalized counts in a sample batch directory -- Argument Parser")

    parser.add_argument("-txtdir", "--txt-dir", type=os.path.abspath, required=True,
                        help="[REQUIRED] Provide the path to the USER input dirctory of tab-delimited txt files for genebodycoverage.\n -txtdir [~/absolute/path/to/tab-delimited-txtdata/directory],\t--txt-dir [~/absolute/path/to/tab-delimited-txtdata/directory]\n")

    parser.add_argument("-out", "--output-directory", required=True, type=os.path.abspath,
                        help="[REQUIRED] Provide the desire output path for the directory containing the plot.\n -out [~/OUTPUT-FILE-PATH/],\t--output-directory [~/OUTPUT-FILE-PATH/]\n")

    parser.add_argument("-name", "--project-name", required=True, 
                        help="[REQUIRED] Name your project for whiich is table is being generated for.")

    args = parser.parse_args()

    _input_txt_dir     = args.txt_dir
    _out_dir         = args.output_directory
    _project_name    = args.project_name
    
    geneBody_summaryTable(_input_txt_dir, _out_dir, _project_name)
