#LIBRARY PATH
libs <- .libPaths("/home/ttm3567/63_tylert/DEseq_libpath_dir/")

library(purrr)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(tidyr)
library(readr)
library(pheatmap)
library(dendextend)
library(MatrixGenerics)
library(stats)


# INPUT for script
args <- commandArgs(TRUE) # allows script to run in CmdLine
TSV_counts_file <- as.character(args[1]) # data used for the matrix - TSV/txt
HQ_samples_in_desired_order_with_ExpGrp <- as.character(args[2]) # needed for figure prep, data cleaning - CSV
genes <- as.character(args[3]) # object of interest for matrix - CSV
k_means_centers <- as.numeric(args[4])
output_directory <- as.character(args[5]) # include '/' at beginning and end of path. E.g. "/path/to/dir/"

#  ================================= ================================= ================================= ================================= ================================= ================================= ================================= ================================= ================================= 
output_directory <- output_directory

counts <- read_tsv(TSV_counts_file)
# 1st column - Symbol
# 2nd column - Gene_Symbol
# column 3: whatever range - ####'s

print(counts[1:10,])

# HQ sample list
HQsamples <- read.csv(HQ_samples_in_desired_order_with_ExpGrp, header = TRUE)
print(HQsamples)
# 1st column - HQ_samples
# 2nd column - Exp_Grp

# curated Gene list or variable gene list
variable_genes <- read.csv(genes, header = TRUE)
# 1 column = Variable_Genes

# filter the FPKM counts based on the genes of interest before selecting the HQ samples - Get the indices of the matching rows
indices <- which(counts$Gene_Symbol %in% variable_genes$Variable_Genes)
print(length(indices))

# Subset the data frame using the indices
counts_variableGenes <- data.frame(counts[indices, ])
print(counts_variableGenes[1:10,])

# Set the row names of counts_variableGenes to the Gene_Symbol column - the table only contains genes of interest at this point
rownames(counts_variableGenes) <- counts_variableGenes$Gene_Symbol

# filter the columns in the counts file based on HQ samples - used as input for hierarchical clustering
counts_variableGenes_HQsamples <- counts_variableGenes[colnames(counts_variableGenes) %in% HQsamples$HQ_samples] 
print(counts_variableGenes_HQsamples[1:10,])

# scale and center the counts so they are k-means clustered on the z-scores
# Perform min-max normalization to range from -1 to 1 across each row ========================================

counts_variableGenes_HQsamples_scaled <- t(scale(t(counts_variableGenes_HQsamples)))  # Convert data to z-scores
print(counts_variableGenes_HQsamples_scaled[1:10,])

min_max_normalization <- function(x, min_value = -1, max_value = 1) {
  min_x <- min(x)
  max_x <- max(x)
  normalized_x <- (x - min_x) / (max_x - min_x) * (max_value - min_value) + min_value
  return(normalized_x)
}

# Apply the min-max normalization function to each row of the data frame
normalized_data <- apply(counts_variableGenes_HQsamples_scaled, 1, min_max_normalization, min_value = -1, max_value = 1)

# Convert the matrix back to a data frame with the same column names
normalized_counts_variableGenes <- t(as.data.frame(normalized_data))
print(normalized_counts_variableGenes[1:10,])

# Genes from variable_genes not found in counts
genes_not_found <- data.frame(Genes_not_in_data = variable_genes$Variable_Genes[!variable_genes$Variable_Genes %in% counts$Gene_Symbol])

# Print the genes not found
print(genes_not_found)
write_csv(genes_not_found, paste0(output_directory,"genes_not_found.csv"))

# k-means CLUSTERS =================================
# Perform k-means clustering and create a dataframe with rownames and their corresponding cluster
set.seed(40)
kmeans_result <- kmeans(normalized_counts_variableGenes, centers=k_means_centers, iter.max=20, nstart=200)
#View(kmeans_result)

# function for cluster assignments 
assignClusterNames <- function(kmeans_result) {
  # Extract the cluster assignments from the kmeans result
  clusters <- data.frame(gene_name = names(kmeans_result$cluster), cluster = paste("Gene Cluster", kmeans_result$cluster))
  
  return(clusters)
}

my_gene_column_kmeans_clusters <- assignClusterNames(kmeans_result)
print(my_gene_column_kmeans_clusters)

# Order the data frame rows according to k-means cluster assignments - 'my_gene_col_<#>clusters' object used 8 times 
my_gene_column_kmeans_clusters <- my_gene_column_kmeans_clusters[order(my_gene_column_kmeans_clusters$cluster), ]

# order the count variable gene rows to match the k-means result
normalized_counts_variableGenesOrder <- normalized_counts_variableGenes[my_gene_column_kmeans_clusters$gene_name, ]

# set the genes for each cluster as the row names and remove the column
rownames(my_gene_column_kmeans_clusters) <- my_gene_column_kmeans_clusters$gene_name
my_gene_column_kmeans_clusters$gene_name <- NULL

#  ================================= ================================= ================================= ================================= ================================= ================================= ================================= ================================= ================================= 

# Create the heatmap
print("about to make heatmap")
HM_normalized_counts_variableGenesOrder <- pheatmap(normalized_counts_variableGenesOrder, annotation_row = my_gene_column_kmeans_clusters,
                                                    cluster_rows = FALSE, 
                                                    cluster_cols = FALSE, 
                                                    fontsize_row = .10,
                                                    main = paste0("K-means Gene Expression Heatmap - Correlation Distance"),
                                                    treeheight_row = 30, 
                                                    angle_col = 315, 
                                                    fontsize_col = 7,
                                                    width = 10, 
                                                    height = 12)

ggsave(filename = paste0(output_directory,"kmeans_corr_dist_HM_normalized_counts_variableGenesOrder",".pdf"),plot = HM_normalized_counts_variableGenesOrder, width = 10, height = 12)
write.csv(my_gene_column_kmeans_clusters,
          file = paste0(output_directory, "kmeans_corr_dist_HM_normalized_counts_variableGenesOrder_CSV",".csv"),
          row.names = TRUE)


# group average heatmap version
print(HQsamples)
colnames(normalized_counts_variableGenesOrder)

calculateGroupMeans <- function(normalized_counts, sample_info) {
  # Create an empty dataframe to store the results
  averaged_counts <- data.frame(row.names = rownames(normalized_counts))
  
  # Loop through each unique experimental group in the sample information
  for(exp_grp in unique(sample_info$Exp_Grp)) {
    # Find the samples belonging to this experimental group
    samples_in_grp <- sample_info$HQ_samples[sample_info$Exp_Grp == exp_grp]
    
    # Subset the normalized counts matrix to only these samples
    counts_subset <- normalized_counts[, colnames(normalized_counts) %in% samples_in_grp, drop = FALSE]
    
    # Calculate the mean across the columns (samples) for each gene
    group_means <- rowMeans(counts_subset, na.rm = TRUE)
    
    # Add the results to the dataframe, naming the column after the experimental group
    averaged_counts[[exp_grp]] <- group_means
  }
  
  return(averaged_counts)
}

averaged_counts_result <- calculateGroupMeans(normalized_counts_variableGenesOrder, HQsamples)

print("about to make heatmap")
HM_normalized_counts_variableGenesOrder_grp <- pheatmap(averaged_counts_result, annotation_row = my_gene_column_kmeans_clusters,
                                                    cluster_rows = FALSE, 
                                                    cluster_cols = FALSE, 
                                                    fontsize_row = .10,
                                                    main = paste0("K-means Group Gene Expression Heatmap - Correlation Distance"),
                                                    treeheight_row = 30, 
                                                    angle_col = 315, 
                                                    fontsize_col = 7,
                                                    width = 10, 
                                                    height = 12)

ggsave(filename = paste0(output_directory,"kmeans_HM_normalized_counts_variableGenesOrder_Grp",".pdf"),plot = HM_normalized_counts_variableGenesOrder_grp, width = 10, height = 12)
# write.csv(my_gene_column_kmeans_clusters,
#           file = paste0(output_directory, "kmeans_HM_normalized_counts_variableGenesOrder_Grp_CSV", ".csv"),
#           row.names = TRUE)

#  ================================= ================================= ================================= ================================= ================================= ================================= ================================= ================================= ================================= 
# Samples Hierarchical Correlation-based distance Heatmap

# call the min-max fxn again
Hier_Clust_normalized_data <- apply(counts_variableGenes_HQsamples, 1, min_max_normalization, min_value = -1, max_value = 1)
Hier_Clust_normalized_data <- t(as.data.frame(Hier_Clust_normalized_data))
print(Hier_Clust_normalized_data[1:10,])


HM_HierClustering_variableGenesOrder <- pheatmap(Hier_Clust_normalized_data,
                                                 cluster_rows = TRUE,
                                                 clustering_distance_rows = "correlation",
                                                 cluster_cols = FALSE, 
                                                 fontsize_row = .10,
                                                 main = paste0("Hierarchical Clustering Gene Expression Heatmap - Correlation Distance"),
                                                 treeheight_row = 30, 
                                                 angle_col = 315, 
                                                 fontsize_col = 7,
                                                 width = 10, 
                                                 height = 12)

ggsave(filename = paste0(output_directory,"Hierarchical_Clustering_HM",".pdf"),plot = HM_HierClustering_variableGenesOrder, width = 10, height = 12)


# Experimental Group Hierarchical Correlation-based distance Heatmap
HierClustering_averaged_counts_result <- calculateGroupMeans(Hier_Clust_normalized_data, HQsamples)

HM_HierClustering_variableGenesOrder_Grp <- pheatmap(HierClustering_averaged_counts_result,
                                                     cluster_rows = TRUE,
                                                     clustering_distance_rows = "correlation",
                                                     cluster_cols = FALSE, 
                                                     fontsize_row = .10,
                                                     main = paste0("Hierarchical Clustering Group Gene Expression Heatmap - Correlation Distance"),
                                                     treeheight_row = 30, 
                                                     angle_col = 315, 
                                                     fontsize_col = 7,
                                                     width = 10, 
                                                     height = 12)

ggsave(filename = paste0(output_directory,"Hierarchical_Clustering_HM_Grp",".pdf"),plot = HM_HierClustering_variableGenesOrder_Grp, width = 10, height = 12)

