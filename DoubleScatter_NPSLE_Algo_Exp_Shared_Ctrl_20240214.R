# Tyler Therron, MS - Winter Lab, Macrophage Genomics
# Rheumatology Department, Feinberg School of Medicine, Northwestern University

#LIBRARY PATH
libs <- .libPaths("/home/ttm3567/63_tylert/DEseq_libpath_dir/")

# load libraries
library(purrr)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(tidyr)
library(readr)
library(pheatmap)
library(dendextend)
library(MatrixGenerics)
library(plotly)
library(htmlwidgets)

# INPUT for script
args <- commandArgs(TRUE) # allows script to run in CmdLine

DEseq_output_CSV <- as.character(args[1]) # data used for the heatmap
FPKM_table <- as.character(args[2]) # genes in the same order as the kmeans heatmap outputs
output_directory <- as.character(args[3]) # where to put the output heatmap figs
sample_list <- as.character(args[4]) # sets the x-axis order on the PDF
# =============================================
# functions for script

# DEseq outputs contain the ref and experimental conditions in the column name - this fxn prevents this from being an issue by renaming
# the Symbol column is read in as ...1 so that needed to be renamed as well
column_renamer <- function(DEG_data) {
  DEG_data <- DEG_data %>%
    rename_with(~"LFC", starts_with("log2FoldChange"))
  
  DEG_data <- DEG_data %>%
    rename_with(~"Symbol", starts_with("...1"))
  return(DEG_data)
}

# only take genes that have a cutoff above LFC of 1 and below -1
filter_LFC <- function(DEG_data) {
  DEG_data <- DEG_data %>% dplyr::filter(LFC > 1 | LFC < -1)
  return(DEG_data)
}

# only take genes with a p-adjusted value of less than 0.05
# p-adjusted is used because there are so many genes that false positives will occur at the high rate otherwise
filter_padj <- function(DEG_filtered) {
  DEG_data <- DEG_filtered %>% dplyr::filter(padj < 0.05)
  return(DEG_data)
}

annotate_genes <- function(master_df, control_only, experimental_only, shared_genes) {
  # Create a new column for the annotation
  master_df$Annotation <- NA
  
  # Annotate 'Control Only'
  master_df$Annotation[master_df$Genes %in% control_only] <- 'Control Only'
  
  # Annotate 'Experimental Only'
  master_df$Annotation[master_df$Genes %in% experimental_only] <- 'Experimental Only'
  
  # Annotate 'Shared'
  master_df$Annotation[master_df$Genes %in% shared_genes] <- 'Shared'
  
  # Remove genes that are not in any of the lists (not annotated)
  master_df <- master_df[!is.na(master_df$Annotation), ]
  
  return(master_df)
}

# Min-max normalization function - improves clarity on heatmap
min_max_normalization <- function(x, min_value = -1, max_value = 1) {
  min_x <- min(x)
  max_x <- max(x)
  normalized_x <- (x - min_x) / (max_x - min_x) * (max_value - min_value) + min_value
  return(normalized_x)
}
# ===========. file read in  ==================================
# this is the deseq outputs to be used for gene selection of tha DEGs
DEseq_output_CSV <- "/Users/ttm3567/Documents/February2024/NPSLE_DoubleScat_QuantCut_BinCut_20240213/double_scatter/Input_CSVs_DESeq_Info/NZBW_ctrl_ifn_adv_hilo_20240215.csv"

# this is the FPKM table with the samples and data for visualization later
FPKM_table <- readr::read_csv("/Users/ttm3567/Documents/February2024/NPSLE_DoubleScat_QuantCut_BinCut_20240213/NZBW_HQ_FPKM_Data.csv")

# this is the samples to include and order by which they should appear - sample_list
samples <- readr::read_csv("/Users/ttm3567/Documents/February2024/nzbw_HQ_samles.csv")
head(FPKM_table)
print(samples)


# location of the output file
output_directory <- "/Users/ttm3567/Documents/February2024/NPSLE_DoubleScat_QuantCut_BinCut_20240213/double_scatter/NZBW/"
print(paste0("output Directory: ", output_directory))
## =============================================

# Assuming the first column of 'samples' dataframe is named 'hq_samples'
sample_names <- samples$hq_samples

# Filter and reorder the columns in FPKM_table based on sample_names
print("Filter and reorder the columns in FPKM_table based on sample_names")
FPKM_table <- FPKM_table %>%
  select(c(Symbol,Gene_Symbol), all_of(sample_names))

print(FPKM_table)
# Read the file containing the paths to the DEseq outputs since there should be multiple for any given one experiment
paths_data <- readr::read_csv(DEseq_output_CSV)
print(paths_data[1:10,])
# Use lapply to read in each CSV specified in the Deseq_Output_Path column
print("Use lapply to read in each CSV specified in the Deseq_Output_Path column")
DEG_data_list <- lapply(1:nrow(paths_data), function(i) {
  readr::read_csv(paths_data$Deseq_Output_Path[i])
})
print(DEG_data_list[[1]][1:10,])

# Set names of the list items based on the Comparison_Name column
names(DEG_data_list) <- paths_data$Comparison_Name
print(DEG_data_list[[1]])
## =============================================

# DESeq Filtering ------------------------------
# Apply the functions to each item in DEG_data_list
print("Check that the genes are the same between both DESeq outputs (they should be)")
all(DEG_data_list[[1]]$Gene_Symbol %in% DEG_data_list[[2]]$Gene_Symbol)
all(DEG_data_list[[2]]$Gene_Symbol %in% DEG_data_list[[1]]$Gene_Symbol)



DEG_Control_data_union_of_genes_clean <- DEG_data_list[[1]] %>% 
  filter(Gene_Symbol != "NoGeneName") %>% 
  distinct(Gene_Symbol, .keep_all = TRUE) %>%
  rename_with(~"LFC_control", starts_with("log2FoldChange"))

DEG_Experimental_data_union_of_genes_clean <- DEG_data_list[[2]] %>% 
  filter(Gene_Symbol != "NoGeneName") %>% 
  distinct(Gene_Symbol, .keep_all = TRUE) %>%
  rename_with(~"LFC_experimental", starts_with("log2FoldChange"))

all(DEG_Control_data_union_of_genes_clean$Gene_Symbol %in% DEG_Experimental_data_union_of_genes_clean$Gene_Symbol)
all(DEG_Experimental_data_union_of_genes_clean$Gene_Symbol %in% DEG_Control_data_union_of_genes_clean$Gene_Symbol)

master_DEG_dataframe <- data.frame(Genes = DEG_Control_data_union_of_genes_clean$Gene_Symbol, 
                                   Control_LFC = DEG_Control_data_union_of_genes_clean$LFC_control,
                                   Experimental_LFC = DEG_Experimental_data_union_of_genes_clean$LFC_experimental)

print(master_DEG_dataframe[1:10,])
print("processing DEseq datasets")
processed_DEG_data_list <- lapply(DEG_data_list, function(data) {
  data <- column_renamer(data)
  data <- filter_LFC(data)
  data <- filter_padj(data)
  return(data)
})

print(processed_DEG_data_list)

# Filter out rows where Gene_Symbol is "NoGeneName"
print(" Filter out rows where Gene_Symbol is NoGeneName")
filtered_DEG_data_list <- lapply(processed_DEG_data_list, function(df) {
  df %>% filter(Gene_Symbol != "NoGeneName")
})

print(filtered_DEG_data_list)

# Apply distinct to each data frame in the list
print("removing duplicated gene symbols")
distinct_DEG_data_list <- lapply(filtered_DEG_data_list, function(df) {
  distinct(df, Gene_Symbol, .keep_all = TRUE)
})
print(distinct_DEG_data_list)

# DESeq Filtering ------------------------------

# DESeq Gene Group Creation ------------------------------

control_hilo <- data.frame(distinct_DEG_data_list[[1]])
experimental_hilo <- data.frame(distinct_DEG_data_list[[2]])

print(control_hilo[1:10,])
print(experimental_hilo[1:10,])

# Significant DEGs in control but not in experimental
genes_in_control_not_in_experimental <- setdiff(control_hilo$Gene_Symbol, experimental_hilo$Gene_Symbol)

# Significant DEGs in experimental but not in control
genes_in_experimental_not_in_control <- setdiff(experimental_hilo$Gene_Symbol, control_hilo$Gene_Symbol)

# Significant DEGs in control AND Experimental
common_values <- intersect(control_hilo$Gene_Symbol, experimental_hilo$Gene_Symbol)

annotated_master_DEG_df <- annotate_genes(master_DEG_dataframe, 
                                      genes_in_control_not_in_experimental, 
                                      genes_in_experimental_not_in_control, 
                                      common_values)

# View the updated dataframe
head(annotated_master_DEG_df)

# get the LFC into one dataframe

# DESeq Gene Group Creation ------------------------------

# Double Scatter Plotting ------------------------------

plot <- plot_ly(data = annotated_master_DEG_df, 
                x = ~Control_LFC, 
                y = ~Experimental_LFC, 
                type = 'scatter', 
                mode = 'markers',
                text = ~paste("Gene:", Genes), # Add gene name to hover text
                hoverinfo = 'text+x+y',
                color = ~Annotation,  # Color points based on Classification
                colors = c("Control Only" = "red", "Shared" = "blue", "Experimental Only" = "green")) # Define custom colors

plot <- plot %>% layout(title = paste0("Scatter plot of Significant DEG LFC values: ",paths_data$Title_Info[1], " ",paths_data$Title_Info[2]),
                        xaxis = list(title = "Control_LFC"),
                        yaxis = list(title = "Experimental_LFC"))
plot

# Double Scatter Plotting ------------------------------

# Save the plot as an interactive HTML file -----------------
htmlwidgets::saveWidget(plot, paste0(output_directory,"HTML_",paths_data$Title_Info[1], "_",paths_data$Title_Info[2],".html"))
write_csv(annotated_master_DEG_df,paste0(output_directory,"CSV_",paths_data$Title_Info[1], "_",paths_data$Title_Info[2],".csv"))



