#!/usr/bin/env Rscript

# USAGE: 1. module load R/4.2.0
# 2. Rscript GeneExpressionMatrix_CmdLine.R [1]fpkm_counts_file.txt [2]genes_of_interest.csv [3]HQ_samples_in_desired_order.csv [4]/output_directory/ [5]model_name [6]regex_pattern_for_splitting_CD11lo_CD11hi_in_sample_HM [7]regex_pattern_for_splitting_CD11lo_CD11hi_in_grp_HM

# Example of REGEX expression to split heatmap: pattern2 <- "^CD11chi_NZBW$"
# ^[string]$ = exact matches only

#LIBRARY PATH
libs <- .libPaths("/Users/ttm3567/Library/R/arm64/4.2/library")

#install.packages("rlang")

#install.packages("purrr")

library(purrr)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(tidyr)
library(readr)
library(pheatmap)
library(dendextend)
library(MatrixGenerics)


# INPUT for script
args <- commandArgs(TRUE) # allows script to run in CmdLine

# Argument inputs
fpkm_counts_file <- as.character(args[1])
genes_cluster_table <- as.character(args[2])
HQ_samples_in_desired_order <- as.character(args[3])
output_directory <- as.character(args[4])
model_name <- as.character(args[5])
#regex_pattern_for_splitting_CD11lo_CD11hi_in_sample_HM <- as.character(args[6])

# FPKM counts that will be min-max normalized
FPKM <- read_tsv(fpkm_counts_file)
print("FPKM counts file")
head(FPKM)
# Read the gene cluster table
gene_cluster_table <- read.csv(genes_cluster_table)
print("Experimental group genes and their respectice cluster assignments")
head(gene_cluster_table)
# HQ sample list
HQsamples <- read.csv(HQ_samples_in_desired_order, header = TRUE)
print("the samples to be analyzed in desired order")
print(HQsamples)
# Get the indices of the matching rows
print("matching FPKM genes to the ones in the gene cluster table")
indices <- which(FPKM$Gene_Symbol %in% gene_cluster_table$X)

# Subset the data frame using the indices
FPKM_GoI <- FPKM[indices, ]

print("filtering and reordering columns in FPKM counts matrix based on desired sample order CSV")
# Filter the counts based on HQ samples
FPKM_HQsamples_GoI <- FPKM_GoI[colnames(FPKM_GoI) %in% HQsamples$HQ_samples]

# Create a vector of the new column order based on the row order in HQsamples
new_order <- HQsamples$HQ_samples %>% match(colnames(FPKM_HQsamples_GoI))

# Reorder the columns
FPKM_HQsamples_GoI_ordered <- FPKM_HQsamples_GoI %>% dplyr::select(all_of(new_order))
head(FPKM_HQsamples_GoI_ordered)

print("normalization across rows (min-max)")
# Min-max normalization function
min_max_normalization <- function(x, min_value = -1, max_value = 1) {
  min_x <- min(x)
  max_x <- max(x)
  normalized_x <- (x - min_x) / (max_x - min_x) * (max_value - min_value) + min_value
  return(normalized_x)
}

# Apply the min-max normalization function to each row
normalized_data <- apply(FPKM_HQsamples_GoI_ordered, 1, min_max_normalization)

# Convert back to a data frame
normalized_FPKM_HQsamples_GoI <- t(as.data.frame(normalized_data))
rownames(normalized_FPKM_HQsamples_GoI) <- FPKM_GoI$Gene_Symbol
head(normalized_FPKM_HQsamples_GoI)

print("order the genes in data to match the gene cluster table")
gene_cluster_table$order <- seq_len(nrow(gene_cluster_table))

# Order the genes according to the gene_cluster_table
merged_data <- merge(gene_cluster_table, normalized_FPKM_HQsamples_GoI, by.x = "X", by.y = "row.names")
print("Debug: merged_data")
print(head(merged_data))
#sorted_data <- merged_data[order(merged_data$Cluster),]
sorted_data <- merged_data[order(merged_data$order), ]

print("Debug: sorted_data")
print(head(sorted_data))

# sorted_data$Cluster <- NULL

# normalized_FPKM_HQsamples_GoI_ordered <- as.data.frame(t(sorted_data))
# colnames(normalized_FPKM_HQsamples_GoI_ordered) <- sorted_data$Gene
# print(head(sorted_data))
# colnames(normalized_FPKM_HQsamples_GoI_ordered)
# head(normalized_FPKM_HQsamples_GoI_ordered)
# normalized_FPKM_HQsamples_GoI_ordered <- as.data.frame(t(normalized_FPKM_HQsamples_GoI_ordered))
# head(normalized_FPKM_HQsamples_GoI_ordered)

# print("regex to split the heatmap for visual ease")
# # Regex to split CD11clo from CD11chi
# pattern <- regex_pattern_for_splitting_CD11lo_CD11hi_in_sample_HM
# matching_indices <- grep(pattern, colnames(sorted_data))
# gaps_indices <- matching_indices - 1

print("making heatmap")
# Create the heatmap
# annots_HM_FPKM_HQsamples_8clusters <- pheatmap(
#   sorted_data,
#   cluster_rows = FALSE, 
#   cluster_cols = FALSE,
#   main = paste0("Gene Expression Matrix - ", model_name),
#   fontsize_row = 0.5,
#   treeheight_row = 30,
#   angle_col = 315,
#   fontsize_col = 7,
#   width = 10,
#   height = 12,
#   gaps_col = gaps_indices
# )

# Exclude the non-sample columns
data_to_plot <- sorted_data[, !(names(sorted_data) %in% c("X", "Cluster", "order"))]

# Set the row names to be the "X" column (the gene names)
rownames(data_to_plot) <- sorted_data$X
head(data_to_plot)

annotation_df <- data.frame(Cluster = sorted_data$Cluster)
rownames(annotation_df) <- sorted_data$X
head(annotation_df)

# Create the heatmap
heatmap_gene_cluster_assignments <- pheatmap(
  data_to_plot,
  cluster_rows = FALSE,  # No reordering of rows
  cluster_cols = FALSE,  # No reordering of columns
  annotation_row = annotation_df, # If you want to use gene names as annotation
  main = paste0("Gene Expression Matrix - ", model_name),
  fontsize_row = 0.5,  # Font size for row labels
  fontsize_col = 5,  # Font size for column labels
  treeheight_row = 40,
  angle_col = 315,
  width = 14,
  height = 10
)

print("saving heatmap")
# Save the heatmap
ggsave(filename = paste0(output_directory, "8clusters_Sample_Gene_Expression_Matrix_Group_GeneOrder_", model_name, ".pdf"), plot = heatmap_gene_cluster_assignments, width = 10, height = 12)


# old version
# INPUT for script
# args <- commandArgs(TRUE) # allows script to run in CmdLine
# fpkm_counts_file <- as.character(args[1]) # data used for the matrix
# genes_cluster_table <- as.character(args[2]) # object of interest for matrix
# HQ_samples_in_desired_order <- as.character(args[3]) # needed for figure prep, data cleaning
# output_directory <- as.character(args[4]) # include '/' at beginning and end of path. E.g. "/path/to/dir/"
# model_name <- as.character(args[5]) # used for titles, and file names
# regex_pattern_for_splitting_CD11lo_CD11hi_in_sample_HM <- as.character(args[6]) # splits columns by cell type
# regex_pattern_for_splitting_CD11lo_CD11hi_in_grp_HM <- as.character(args[7]) # splits columns by cell type
# #  ================================= ================================= ================================= ================================= ================================= ================================= ================================= ================================= ================================= 
# output_directory <- output_directory
# # FPKM counts that will be min-max normalized
# FPKM <- read_tsv(fpkm_counts_file) # arg1
# 
# # GoI list
# genes_cluster_table <- read.csv(genes_cluster_table) # arg2
# 
# # HQ sample list
# HQsamples <- read.csv(HQ_samples_in_desired_order, header = TRUE) # arg3
# 
# # HQ samples should also contain the experimental group, so I can easily average the values in prior to plotting the data in HEATMAP 3
# 
# # filter the FPKM counts based on the genes of interest before selecting the HQ samples
# # Get the indices of the matching rows
# indices <- which(FPKM$Gene_Symbol %in% GoI_Cuda$GoI)
# 
# # Subset the data frame using the indices
# FPKM_GoI <- FPKM[indices, ]
# # filter the counts based on HQ samples
# FPKM_HQsamples_GoI <- FPKM_GoI[colnames(FPKM_GoI) %in% HQsamples$HQ_samples] # count file needs to be order CD11cLOW C8, CD11cLOW fingo, etc, CD11cHI C8, etc.
# 
# # order column names to fit the specified organization as put in the HQsamples df
# # Extract column names from df1
# FPKM_HQsamples_GoI_colnames <- colnames(FPKM_HQsamples_GoI)
# 
# # Create a vector of the new column order based on the row order of the same names in df2
# new_order <- HQsamples$HQ_samples %>% match(FPKM_HQsamples_GoI_colnames)
# # Reorder the columns in df1 based on the new column order vector
# FPKM_HQsamples_GoI_ordered <- FPKM_HQsamples_GoI %>% dplyr::select(all_of(new_order))
# 
# # normalize the selected data before adding the gene symbol column from 'FPKM_GoI'
# # Perform min-max normalization to range from -1 to 1 across each row ========================================
# 
# min_max_normalization <- function(x, min_value = -1, max_value = 1) {
#   min_x <- min(x)
#   max_x <- max(x)
#   normalized_x <- (x - min_x) / (max_x - min_x) * (max_value - min_value) + min_value
#   return(normalized_x)
# }
# 
# # Apply the min-max normalization function to each row of the data frame
# normalized_data <- apply(FPKM_HQsamples_GoI_ordered, 1, min_max_normalization, min_value = -1, max_value = 1)
# 
# # Convert the matrix back to a data frame with the same column names
# normalized_FPKM_HQsamples_GoI <- t(as.data.frame(normalized_data))
# rownames(normalized_FPKM_HQsamples_GoI) <- FPKM_GoI$Gene_Symbol
# 
# # Omit any NA values - 5/17/2023 b/c certain gene lists may not contain the genes listed in the bulk pipeline output if it is a GoI list found from the literature
# #normalized_FPKM_HQsamples_GoI <- na.omit(normalized_FPKM_HQsamples_GoI)  # This removes rows with any NA values
# 
# # The data table has the HQ samples selected, is filtered by the desired gene list, and has had every row min-max normalized from -1 to 1 +++++++++++++++++++++++++++++======================++++++++++++++++++============================
# # Next: create the correlation figure - absolute scale -1 to 1. Pearson's. Blue[low] to red[high].
# 
# 
# # REGEX1 to split CD11clo from CD11chi for visual ease
# pattern <-  regex_pattern_for_splitting_CD11lo_CD11hi_in_sample_HM # Replace with your desired regex pattern
# matching_indices <- grep(pattern, colnames(normalized_FPKM_HQsamples_GoI))
# gaps_indices <- matching_indices - 1
# gaps_indices
# 
# 
# # HEATMAP 1 - clustered samples, columns should remain in specified order
# 
# annots_HM_FPKM_HQsamples_8clusters <- pheatmap(
#   normalized_FPKM_HQsamples_GoI_ordered,
#   cluster_rows = FALSE,  # Setting this to FALSE ensures that rows are not reordered
#   cluster_cols = FALSE,  # Setting this to FALSE ensures that columns are not reordered
#   main = paste0("Gene Expression Matrix - ", model_name),
#   fontsize_row = 0.5,
#   treeheight_row = 30,
#   angle_col = 315,
#   fontsize_col = 7,
#   width = 10,
#   height = 12,
#   gaps_col = gaps_indices
# )
# 
# ggsave(filename = paste0(output_directory,"kmeans_8clusters_Sample_Gene_Expression_Matrix_",model_name,".pdf"),plot = annots_HM_FPKM_HQsamples_8clusters, width = 10, height = 12)
# 
# 
