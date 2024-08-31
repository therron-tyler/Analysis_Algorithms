# USAGE: 1. module load R/4.2.0
# 2. Rscript GeneExpressionMatrix_CmdLine.R [1]fpkm_counts_file.txt [2]genes_of_interest.csv [3]HQ_samples_in_desired_order.csv [4]/output_directory/ [5]model_name [6]regex_pattern_for_splitting_CD11lo_CD11hi_in_sample_HM [7]regex_pattern_for_splitting_CD11lo_CD11hi_in_grp_HM

# Example of REGEX expression to split heatmap: pattern2 <- "^CD11chi_NZBW$"
# ^[string]$ = exact matches only

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


# INPUT for script
args <- commandArgs(TRUE) # allows script to run in CmdLine
fpkm_counts_file <- as.character(args[1]) # data used for the matrix
genes_of_interest <- as.character(args[2]) # object of interest for matrix
HQ_samples_in_desired_order <- as.character(args[3]) # needed for figure prep, data cleaning
output_directory <- as.character(args[4]) # include '/' at beginning and end of path. E.g. "/path/to/dir/"
model_name <- as.character(args[5]) # used for titles, and file names
regex_pattern_for_splitting_CD11lo_CD11hi_in_sample_HM <- as.character(args[6]) # splits columns by cell type
regex_pattern_for_splitting_CD11lo_CD11hi_in_grp_HM <- as.character(args[7]) # splits columns by cell type
#  ================================= ================================= ================================= ================================= ================================= ================================= ================================= ================================= ================================= 
output_directory <- output_directory
# FPKM counts that will be min-max normalized
FPKM <- read_tsv(fpkm_counts_file)

# GoI list
GoI_Cuda <- unique(read.csv(genes_of_interest))

# HQ sample list
HQsamples <- read.csv(HQ_samples_in_desired_order, header = TRUE)

# HQ samples should also contain the experimental group, so I can easily average the values in prior to plotting the data in HEATMAP 3

# filter the FPKM counts based on the genes of interest before selecting the HQ samples
# Get the indices of the matching rows
indices <- which(FPKM$Gene_Symbol %in% GoI_Cuda$GoI)

# Subset the data frame using the indices
FPKM_GoI <- FPKM[indices, ]
# filter the counts based on HQ samples
FPKM_HQsamples_GoI <- FPKM_GoI[colnames(FPKM_GoI) %in% HQsamples$HQ_samples] # count file needs to be order CD11cLOW C8, CD11cLOW fingo, etc, CD11cHI C8, etc.

# order column names to fit the specified organization as put in the HQsamples df
# Extract column names from FPKM_HQsamples_GoI
FPKM_HQsamples_GoI_colnames <- colnames(FPKM_HQsamples_GoI)

# Create a vector of the new column order based on the row order of the same names in FPKM_HQsamples_GoI_colnames
new_order <- HQsamples$HQ_samples %>% match(FPKM_HQsamples_GoI_colnames)
# Reorder the columns in FPKM_HQsamples_GoI based on the new column order vector - 'new_order'
FPKM_HQsamples_GoI_ordered <- FPKM_HQsamples_GoI %>% dplyr::select(all_of(new_order))

# normalize the selected data before adding the gene symbol column from 'FPKM_GoI'
# Perform min-max normalization to range from -1 to 1 across each row ========================================

min_max_normalization <- function(x, min_value = -1, max_value = 1) {
  min_x <- min(x)
  max_x <- max(x)
  normalized_x <- (x - min_x) / (max_x - min_x) * (max_value - min_value) + min_value
  return(normalized_x)
}

# Apply the min-max normalization function to each row of the data frame
normalized_data <- apply(FPKM_HQsamples_GoI_ordered, 1, min_max_normalization, min_value = -1, max_value = 1)

# Convert the matrix back to a data frame with the same column names
normalized_FPKM_HQsamples_GoI <- t(as.data.frame(normalized_data))
rownames(normalized_FPKM_HQsamples_GoI) <- FPKM_GoI$Gene_Symbol

# The data table has the HQ samples selected, is filtered by the desired gene list, and has had every row min-max normalized from -1 to 1 +++++++++++++++++++++++++++++======================++++++++++++++++++============================
# Next: create the correlation figure - absolute scale -1 to 1. Pearson's. Blue[low] to red[high].

# by sample - 2 clusters
set.seed(40)
#kmeans_gene_normalized_FPKM_HQsamples_2clusters <- kmeans(normalized_FPKM_HQsamples_GoI, centers = 2,iter.max = 10, nstart = 200)
# by sample - 4 clusters
#kmeans_gene_normalized_FPKM_HQsamples_4clusters <- kmeans(normalized_FPKM_HQsamples_GoI, centers = 4,iter.max = 10, nstart = 200)

# annotation 2 clusters
# my_gene_col_2clusters <- data.frame(Cluster = ifelse(test = kmeans_gene_normalized_FPKM_HQsamples_2clusters$cluster == 1, yes = "Gene Cluster 1", no = "Gene Cluster 2"))
# head(my_gene_col_2clusters)

# annotation 4 clusters
#my_gene_col_4clusters <- data.frame(Cluster = ifelse(test = kmeans_gene_normalized_FPKM_HQsamples_4clusters$cluster == 1,
#                                                     yes = "Gene Cluster 1",
#                                                     no = ifelse(test = kmeans_gene_normalized_FPKM_HQsamples_4clusters$cluster == 2,
#                                                                 yes = "Gene Cluster 2",
#                                                                 no = ifelse(test = kmeans_gene_normalized_FPKM_HQsamples_4clusters$cluster == 3,
#                                                                             yes = "Gene Cluster 3",
#                                                                             no = "Gene Cluster 4"))))
#head(my_gene_col_4clusters)

# write.csv(my_gene_col_2clusters, 
#           file = paste0(output_directory, "Sample_Kmeans_2_Gene_Clusters_",model_name, ".csv"), 
#           row.names = TRUE)

#write.csv(my_gene_col_4clusters, 
#          file = paste0(output_directory, "Sample_Kmeans_4_Gene_Clusters_",model_name, ".csv"), 
#          row.names = TRUE)

# REGEX1 to split CD11clo from CD11chi for visual ease
pattern <-  regex_pattern_for_splitting_CD11lo_CD11hi_in_sample_HM # Replace with your desired regex pattern
matching_indices <- grep(pattern, colnames(normalized_FPKM_HQsamples_GoI))
gaps_indices <- matching_indices - 1
gaps_indices


# HEATMAP 1 - clustered samples, columns should remain in specified order
# annots_HM_FPKM_HQsamples_2clusters <- pheatmap(normalized_FPKM_HQsamples_GoI, annotation_row = my_gene_col_2clusters, cluster_cols = FALSE,cluster_rows = FALSE, cutree_rows = 2, fontsize_row = .5, main = paste0("K-means Gene Expression Matrix (2 clusters) - ",model_name), treeheight_row = 30, angle_col = 315, fontsize_col = 7, width = 10, height = 12, gaps_col = gaps_indices)
# ggsave(filename = paste0(output_directory,"kmeans_2clusters_Sample_Gene_Expression_Matrix_",model_name,".pdf"),plot = annots_HM_FPKM_HQsamples_2clusters, width = 10, height = 12)


# NEW HEATMAP CODE *********** 7/12/2023 *******

# 2 CLUSTERS  =================================
# Perform k-means clustering and create a dataframe with rownames and their corresponding cluster
kmeans_2_result <- kmeans(normalized_FPKM_HQsamples_GoI, centers=2, iter.max=10, nstart=200)
my_gene_col_2clusters <- data.frame(RowName = rownames(normalized_FPKM_HQsamples_GoI), 
                                    Cluster = ifelse(test = kmeans_2_result$cluster == 1, yes = "Gene Cluster 1",
                                                     no = "Gene Cluster 2"))

head(my_gene_col_2clusters)


# Order the data frame rows according to k-means cluster assignments
my_gene_col_2clusters <- my_gene_col_2clusters[order(my_gene_col_2clusters$Cluster), ]
normalized_FPKM_HQsamples_GoI_ordered <- normalized_FPKM_HQsamples_GoI[my_gene_col_2clusters$RowName, ]
rownames(my_gene_col_2clusters) <- my_gene_col_2clusters$RowName
my_gene_col_2clusters$RowName <- NULL
# Create the heatmap
print("about to make heatmap")
annots_HM_FPKM_HQsamples_2clusters <- pheatmap(normalized_FPKM_HQsamples_GoI_ordered, annotation_row = my_gene_col_2clusters, 
                                               cluster_rows = FALSE, cluster_cols = FALSE, 
                                               cutree_rows = 2, fontsize_row = .5, 
                                               main = paste0("K-means Gene Expression Matrix (2 clusters) - ",model_name), 
                                               treeheight_row = 30, angle_col = 315, fontsize_col = 7, 
                                               width = 10, height = 12, gaps_col = gaps_indices)
ggsave(filename = paste0(output_directory,"kmeans_2clusters_Sample_Gene_Expression_Matrix_",model_name,".pdf"),plot = annots_HM_FPKM_HQsamples_2clusters, width = 10, height = 12)
write.csv(my_gene_col_2clusters,
          file = paste0(output_directory, "Sample_Kmeans_2_Gene_Clusters_",model_name, ".csv"),
          row.names = TRUE)

# 4 CLUSTERS =================================
# Perform k-means clustering and create a dataframe with rownames and their corresponding cluster
kmeans_result <- kmeans(normalized_FPKM_HQsamples_GoI, centers=4, iter.max=10, nstart=200)
my_gene_col_4clusters <- data.frame(RowName = rownames(normalized_FPKM_HQsamples_GoI),
					Cluster = ifelse(test = kmeans_result$cluster == 1, yes = "Gene Cluster 1",
                                                     no = ifelse(test = kmeans_result$cluster == 2, yes = "Gene Cluster 2",
                                                                 no = ifelse(test = kmeans_result$cluster == 3, yes = "Gene Cluster 3",
                                                                             no = "Gene Cluster 4"))))

head(my_gene_col_4clusters)


# Order the data frame rows according to k-means cluster assignments
my_gene_col_4clusters <- my_gene_col_4clusters[order(my_gene_col_4clusters$Cluster), ]
normalized_FPKM_HQsamples_GoI_ordered <- normalized_FPKM_HQsamples_GoI[my_gene_col_4clusters$RowName, ]
rownames(my_gene_col_4clusters) <- my_gene_col_4clusters$RowName
my_gene_col_4clusters$RowName <- NULL
# Create the heatmap
print("about to make heatmap")
annots_HM_FPKM_HQsamples_4clusters <- pheatmap(normalized_FPKM_HQsamples_GoI_ordered, annotation_row = my_gene_col_4clusters, 
                                               cluster_rows = FALSE, cluster_cols = FALSE, 
                                               cutree_rows = 4, fontsize_row = .5, 
                                               main = paste0("K-means Gene Expression Matrix (4 clusters) - ",model_name), 
                                               treeheight_row = 30, angle_col = 315, fontsize_col = 7, 
                                               width = 10, height = 12, gaps_col = gaps_indices)
ggsave(filename = paste0(output_directory,"kmeans_4clusters_Sample_Gene_Expression_Matrix_",model_name,".pdf"),plot = annots_HM_FPKM_HQsamples_4clusters, width = 10, height = 12)
write.csv(my_gene_col_4clusters,
          file = paste0(output_directory, "Sample_Kmeans_4_Gene_Clusters_",model_name, ".csv"),
          row.names = TRUE)

# 6 CLUSTERS =================================
# Perform k-means clustering and create a dataframe with rownames and their corresponding cluster
kmeans_6_result <- kmeans(normalized_FPKM_HQsamples_GoI, centers=6, iter.max=10, nstart=200)
# assign cluster names from the kmeans() function
my_gene_col_6clusters <- data.frame(RowName = rownames(normalized_FPKM_HQsamples_GoI),
                                    Cluster = ifelse(test = kmeans_6_result$cluster == 1, yes = "Gene Cluster 1",
                                                     no = ifelse(test = kmeans_6_result$cluster == 2, yes = "Gene Cluster 2",
                                                                 no = ifelse(test = kmeans_6_result$cluster == 3, yes = "Gene Cluster 3",
                                                                             no = ifelse(test = kmeans_6_result$cluster == 4, yes = "Gene Cluster 4",
                                                                                         no = ifelse(test = kmeans_6_result$cluster == 5, yes = "Gene Cluster 5",
                                                                                                     no = "Gene Cluster 6"))))))

head(my_gene_col_6clusters)


# Order the data frame rows according to k-means cluster assignments - 'my_gene_col_<#>clusters' object used 8 times 
my_gene_col_6clusters <- my_gene_col_6clusters[order(my_gene_col_6clusters$Cluster), ]
normalized_FPKM_HQsamples_GoI_ordered <- normalized_FPKM_HQsamples_GoI[my_gene_col_6clusters$RowName, ]
rownames(my_gene_col_6clusters) <- my_gene_col_6clusters$RowName
my_gene_col_6clusters$RowName <- NULL
# Create the heatmap
print("about to make heatmap")
annots_HM_FPKM_HQsamples_6clusters <- pheatmap(normalized_FPKM_HQsamples_GoI_ordered, annotation_row = my_gene_col_6clusters, 
                                               cluster_rows = FALSE, cluster_cols = FALSE, 
                                               cutree_rows = 6, fontsize_row = .5, 
                                               main = paste0("K-means Gene Expression Matrix (6 clusters) - ",model_name), 
                                               treeheight_row = 30, angle_col = 315, fontsize_col = 7, 
                                               width = 10, height = 12, gaps_col = gaps_indices)
ggsave(filename = paste0(output_directory,"kmeans_6clusters_Sample_Gene_Expression_Matrix_",model_name,".pdf"),plot = annots_HM_FPKM_HQsamples_6clusters, width = 10, height = 12)
write.csv(my_gene_col_6clusters,
          file = paste0(output_directory, "Sample_Kmeans_6_Gene_Clusters_",model_name, ".csv"),
          row.names = TRUE)

# 8 CLUSTERS =================================
# Perform k-means clustering and create a dataframe with rownames and their corresponding cluster
kmeans_8_result <- kmeans(normalized_FPKM_HQsamples_GoI, centers=8, iter.max=10, nstart=200)
# assign cluster names from the kmeans() function
my_gene_col_8clusters <- data.frame(RowName = rownames(normalized_FPKM_HQsamples_GoI),
                                    Cluster = ifelse(test = kmeans_8_result$cluster == 1, yes = "Gene Cluster 1",
                                                     no = ifelse(test = kmeans_8_result$cluster == 2, yes = "Gene Cluster 2",
                                                                 no = ifelse(test = kmeans_8_result$cluster == 3, yes = "Gene Cluster 3",
                                                                             no = ifelse(test = kmeans_8_result$cluster == 4, yes = "Gene Cluster 4",
                                                                                         no = ifelse(test = kmeans_8_result$cluster == 5, yes = "Gene Cluster 5",
                                                                                                     no = ifelse(test = kmeans_8_result$cluster == 6, yes = "Gene Cluster 6",
                                                                                                                 no = ifelse(test = kmeans_8_result$cluster == 7, yes = "Gene Cluster 7",
                                                                                                                             no = "Gene Cluster 8"))))))))

head(my_gene_col_8clusters)


# Order the data frame rows according to k-means cluster assignments - 'my_gene_col_<#>clusters' object used 8 times 
my_gene_col_8clusters <- my_gene_col_8clusters[order(my_gene_col_8clusters$Cluster), ]
normalized_FPKM_HQsamples_GoI_ordered <- normalized_FPKM_HQsamples_GoI[my_gene_col_8clusters$RowName, ]
rownames(my_gene_col_8clusters) <- my_gene_col_8clusters$RowName
my_gene_col_8clusters$RowName <- NULL
# Create the heatmap
print("about to make heatmap")
annots_HM_FPKM_HQsamples_8clusters <- pheatmap(normalized_FPKM_HQsamples_GoI_ordered, annotation_row = my_gene_col_8clusters, 
                                               cluster_rows = FALSE, cluster_cols = FALSE, 
                                               cutree_rows = 8, fontsize_row = .5, 
                                               main = paste0("K-means Gene Expression Matrix (8 clusters) - ",model_name), 
                                               treeheight_row = 30, angle_col = 315, fontsize_col = 7, 
                                               width = 10, height = 12, gaps_col = gaps_indices)
ggsave(filename = paste0(output_directory,"kmeans_8clusters_Sample_Gene_Expression_Matrix_",model_name,".pdf"),plot = annots_HM_FPKM_HQsamples_8clusters, width = 10, height = 12)
write.csv(my_gene_col_8clusters,
          file = paste0(output_directory, "Sample_Kmeans_8_Gene_Clusters_",model_name, ".csv"),
          row.names = TRUE)



# HEATMAP EXPERIMENTAL GROUP - sample group averages, column order same as specified, samples are clustered ++++++++++++++++++++++++++++++++++++++++++++++++++


# use the selected genes FPKM counts df, and the HQsamples df

# match column names to row entries in HQsamples df, then rename the dataframe by the experimental group instead. Make a copy of FPKM counts df first as to not mix up the cols.
# use normalized counts 
HEATMAP3_fpkm <- normalized_FPKM_HQsamples_GoI

# CHECKS to make sure the columns are in the same order from the FPKM counts file to the rows in the HQsamples - which is where the experimental group info is contained

all(colnames(HEATMAP3_fpkm) %in% HQsamples$HQ_samples)
all(colnames(HEATMAP3_fpkm) == HQsamples$HQ_samples)

ifelse(colnames(HEATMAP3_fpkm) == HQsamples$HQ_samples, colnames(HEATMAP3_fpkm) <- HQsamples$Exp_Grp, colnames(HEATMAP3_fpkm))

# average the values of all the columns based on their column name
# created list to subset the larger dataframe
Experimental_grps <- unique(HQsamples$Exp_Grp)
Experimental_grps
subset_list <- list()
for (group in Experimental_grps) {
  # Get the indices of the rows in df2 that belong to the current group
  group_indices <- which(HQsamples$Exp_Grp == group)
  
  # Subset df1 using the group_indices
  subset <- HEATMAP3_fpkm[, group_indices]
  
  # Add the subset to the list with a key corresponding to the group name
  subset_list[[group]] <- subset
}

# the experimental groups each have their own matrix of values - next step is to convert to use rowMeans2()
avg_exp_groups <- lapply(subset_list, function(x) as.data.frame(rowMeans2(x)))

labeled_genes_avg_exp_grps <- lapply(seq_along(avg_exp_groups), function(x) {
  df <- avg_exp_groups[[x]]
  rownames(df) <- rownames(normalized_FPKM_HQsamples_GoI)
  return(df)
})
names(labeled_genes_avg_exp_grps) <- Experimental_grps 

# Merge the data frames into one larger data frame
merged_exp_grps_df <- bind_cols(labeled_genes_avg_exp_grps)
colnames(merged_exp_grps_df) <- Experimental_grps
head(merged_exp_grps_df)

# Generate the Gene Clusters for the Experimental Group averaged heatmap
set.seed(40)
# kmeans_gene_merged_exp_grps_df_4clusters <- kmeans(merged_exp_grps_df, centers = 4,iter.max = 10, nstart = 200)
# 
# # make two types of gene expression output, one with 2 clusters and one with 4 clusters
# 
# merged_exp_grps_gene_col_4clusters <- data.frame(Cluster = ifelse(test = kmeans_gene_merged_exp_grps_df_4clusters$cluster == 1,
#                                                                   yes = "Gene Cluster 1",
#                                                                   no = ifelse(test = kmeans_gene_merged_exp_grps_df_4clusters$cluster == 2,
#                                                                               yes = "Gene Cluster 2",
#                                                                               no = ifelse(test = kmeans_gene_merged_exp_grps_df_4clusters$cluster == 3,
#                                                                                           yes = "Gene Cluster 3",
#                                                                                           no = "Gene Cluster 4"))))
# 
# # 2 clusters
# kmeans_gene_merged_exp_grps_df_2clusters <- kmeans(merged_exp_grps_df, centers = 2,iter.max = 10, nstart = 200)
# merged_exp_grps_gene_col_2clusters <- data.frame(Cluster = ifelse(test = kmeans_gene_merged_exp_grps_df_2clusters$cluster == 1, yes = "Gene Cluster 1", no = "Gene Cluster 2"))
# 
# write.csv(merged_exp_grps_gene_col_4clusters, 
#           file = paste0(output_directory,"ExpGrps_Kmeans_4_Gene_Clusters_",model_name, ".csv"), 
#           row.names = TRUE)
# 
# write.csv(merged_exp_grps_gene_col_2clusters, 
#           file = paste0(output_directory,"ExpGrps_Kmeans_2_Gene_Clusters_",model_name, ".csv"), 
#           row.names = TRUE)

# REGEX2 to split CD11clo from CD11chi for visual ease
pattern2 <- regex_pattern_for_splitting_CD11lo_CD11hi_in_grp_HM # Replace with your desired regex pattern
matching_indices2 <- grep(pattern2, colnames(merged_exp_grps_df))
gaps_indices2 <- matching_indices2 - 1
gaps_indices2


# plot the averaged, normalized, and merged gene expression data
# 4 clusters
# ExpGrp_HeatMap_4clusters <- pheatmap(merged_exp_grps_df, annotation_row = merged_exp_grps_gene_col_4clusters, cluster_cols = FALSE, cutree_rows = 4, fontsize_row = .5, main = paste0("K-means Group Gene Expression Matrix (4 clusters) - ",model_name), treeheight_row = 30, angle_col = 315, fontsize_col = 7, width = 10, height = 12, gaps_col = gaps_indices2)
# ggsave(filename = paste0(output_directory,"kmeans_4clusters_Group_Gene_Expression_Matrix_",model_name,".pdf"),plot = ExpGrp_HeatMap_4clusters, width = 8, height = 12)
# 
# 
# # 2 clusters
# ExpGrp_HeatMap_2clusters <- pheatmap(merged_exp_grps_df, annotation_row = merged_exp_grps_gene_col_2clusters, cluster_cols = FALSE, cutree_rows = 2, fontsize_row = .5, main = paste0("K-means Group Gene Expression Matrix (2 clusters) - ",model_name), treeheight_row = 30, angle_col = 315, fontsize_col = 7, width = 10, height = 12, gaps_col = gaps_indices2)
# ggsave(filename = paste0(output_directory,"kmeans_2clusters_Group_Gene_Expression_Matrix_",model_name,".pdf"),plot = ExpGrp_HeatMap_2clusters, width = 8, height = 12)


# Generating Experimental Group Clusters ***************************************

# 2 CLUSTERS  =================================
# Perform k-means clustering and create a dataframe with rownames and their corresponding cluster
ExpGrps_kmeans_2_result <- kmeans(merged_exp_grps_df, centers=2, iter.max=10, nstart=200)
ExpGrps_my_gene_col_2clusters <- data.frame(RowName = rownames(merged_exp_grps_df), 
                                    Cluster = ifelse(test = ExpGrps_kmeans_2_result$cluster == 1, yes = "Gene Cluster 1",
                                                     no = "Gene Cluster 2"))

head(ExpGrps_my_gene_col_2clusters)


# Order the data frame rows according to k-means cluster assignments
ExpGrps_my_gene_col_2clusters <- ExpGrps_my_gene_col_2clusters[order(ExpGrps_my_gene_col_2clusters$Cluster), ]
merged_exp_grps_df_ordered <- merged_exp_grps_df[ExpGrps_my_gene_col_2clusters$RowName, ]
rownames(ExpGrps_my_gene_col_2clusters) <- ExpGrps_my_gene_col_2clusters$RowName
ExpGrps_my_gene_col_2clusters$RowName <- NULL
# Create the heatmap
print("about to make heatmap")
HM_merged_exp_grps_df_ordered_2clusters <- pheatmap(merged_exp_grps_df_ordered, annotation_row = ExpGrps_my_gene_col_2clusters, 
                                               cluster_rows = FALSE, cluster_cols = FALSE, 
                                               cutree_rows = 2, fontsize_row = .5, 
                                               main = paste0("K-means Group Gene Expression Matrix (2 clusters) - ",model_name), 
                                               treeheight_row = 30, angle_col = 315, fontsize_col = 7, 
                                               width = 10, height = 12, gaps_col = gaps_indices2)
ggsave(filename = paste0(output_directory,"kmeans_2clusters_Group_Gene_Expression_Matrix_",model_name,".pdf"),plot = HM_merged_exp_grps_df_ordered_2clusters, width = 10, height = 12)
write.csv(ExpGrps_my_gene_col_2clusters,
          file = paste0(output_directory, "Group_Kmeans_2_Gene_Clusters_",model_name, ".csv"),
          row.names = TRUE)

# 4 CLUSTERS  =================================
# Perform k-means clustering and create a dataframe with rownames and their corresponding cluster
ExpGrps_kmeans_4_result <- kmeans(merged_exp_grps_df, centers=4, iter.max=10, nstart=200)
ExpGrps_my_gene_col_4clusters <- data.frame(RowName = rownames(merged_exp_grps_df), 
                                            Cluster = ifelse(test = ExpGrps_kmeans_4_result$cluster == 1, yes = "Gene Cluster 1",
                                                             no = ifelse(test = ExpGrps_kmeans_4_result$cluster == 2, yes = "Gene Cluster 2",
                                                                         no = ifelse(test = ExpGrps_kmeans_4_result$cluster == 3, yes = "Gene Cluster 3",
                                                                                     no = "Gene Cluster 4"))))

head(ExpGrps_my_gene_col_4clusters)


# Order the data frame rows according to k-means cluster assignments
ExpGrps_my_gene_col_4clusters <- ExpGrps_my_gene_col_4clusters[order(ExpGrps_my_gene_col_4clusters$Cluster), ]
merged_exp_grps_df_ordered <- merged_exp_grps_df[ExpGrps_my_gene_col_4clusters$RowName, ]
rownames(ExpGrps_my_gene_col_4clusters) <- ExpGrps_my_gene_col_4clusters$RowName
ExpGrps_my_gene_col_4clusters$RowName <- NULL
# Create the heatmap
print("about to make heatmap")
HM_merged_exp_grps_df_ordered_4clusters <- pheatmap(merged_exp_grps_df_ordered, annotation_row = ExpGrps_my_gene_col_4clusters, 
                                                    cluster_rows = FALSE, cluster_cols = FALSE, 
                                                    cutree_rows = 4, fontsize_row = .5, 
                                                    main = paste0("K-means Group Gene Expression Matrix (4 clusters) - ",model_name), 
                                                    treeheight_row = 30, angle_col = 315, fontsize_col = 7, 
                                                    width = 10, height = 12, gaps_col = gaps_indices2)
ggsave(filename = paste0(output_directory,"kmeans_4clusters_Group_Gene_Expression_Matrix_",model_name,".pdf"),plot = HM_merged_exp_grps_df_ordered_4clusters, width = 10, height = 12)
write.csv(ExpGrps_my_gene_col_4clusters,
          file = paste0(output_directory, "Group_Kmeans_4_Gene_Clusters_",model_name, ".csv"),
          row.names = TRUE)

# 6 CLUSTERS  =================================
# Perform k-means clustering and create a dataframe with rownames and their corresponding cluster
ExpGrps_kmeans_6_result <- kmeans(merged_exp_grps_df, centers=6, iter.max=10, nstart=200)
ExpGrps_my_gene_col_6clusters <- data.frame(RowName = rownames(merged_exp_grps_df), 
                                            Cluster = ifelse(test = ExpGrps_kmeans_6_result$cluster == 1, yes = "Gene Cluster 1",
                                                             no = ifelse(test = ExpGrps_kmeans_6_result$cluster == 2, yes = "Gene Cluster 2",
                                                                         no = ifelse(test = ExpGrps_kmeans_6_result$cluster == 3, yes = "Gene Cluster 3",
                                                                                     no = ifelse(test = ExpGrps_kmeans_6_result$cluster == 4, yes = "Gene Cluster 4",
                                                                                                 no = ifelse(test = ExpGrps_kmeans_6_result$cluster == 5, yes = "Gene Cluster 5",
                                                                                                             no = "Gene Cluster 6"))))))

head(ExpGrps_my_gene_col_6clusters)


# Order the data frame rows according to k-means cluster assignments
ExpGrps_my_gene_col_6clusters <- ExpGrps_my_gene_col_6clusters[order(ExpGrps_my_gene_col_6clusters$Cluster), ]
merged_exp_grps_df_ordered <- merged_exp_grps_df[ExpGrps_my_gene_col_6clusters$RowName, ]
rownames(ExpGrps_my_gene_col_6clusters) <- ExpGrps_my_gene_col_6clusters$RowName
ExpGrps_my_gene_col_6clusters$RowName <- NULL
# Create the heatmap
print("about to make heatmap")
HM_merged_exp_grps_df_ordered_6clusters <- pheatmap(merged_exp_grps_df_ordered, annotation_row = ExpGrps_my_gene_col_6clusters, 
                                                    cluster_rows = FALSE, cluster_cols = FALSE, 
                                                    cutree_rows = 6, fontsize_row = .5, 
                                                    main = paste0("K-means Group Gene Expression Matrix (6 clusters) - ",model_name), 
                                                    treeheight_row = 30, angle_col = 315, fontsize_col = 7, 
                                                    width = 10, height = 12, gaps_col = gaps_indices2)
ggsave(filename = paste0(output_directory,"kmeans_6clusters_Group_Gene_Expression_Matrix_",model_name,".pdf"),plot = HM_merged_exp_grps_df_ordered_6clusters, width = 10, height = 12)
write.csv(ExpGrps_my_gene_col_6clusters,
          file = paste0(output_directory, "Group_Kmeans_6_Gene_Clusters_",model_name, ".csv"),
          row.names = TRUE)

# 8 CLUSTERS  =================================
# Perform k-means clustering and create a dataframe with rownames and their corresponding cluster
ExpGrps_kmeans_8_result <- kmeans(merged_exp_grps_df, centers=8, iter.max=10, nstart=200)
ExpGrps_my_gene_col_8clusters <- data.frame(RowName = rownames(merged_exp_grps_df), 
                                            Cluster = ifelse(test = ExpGrps_kmeans_8_result$cluster == 1, yes = "Gene Cluster 1",
                                                             no = ifelse(test = ExpGrps_kmeans_8_result$cluster == 2, yes = "Gene Cluster 2",
                                                                         no = ifelse(test = ExpGrps_kmeans_8_result$cluster == 3, yes = "Gene Cluster 3",
                                                                                     no = ifelse(test = ExpGrps_kmeans_8_result$cluster == 4, yes = "Gene Cluster 4",
                                                                                                 no = ifelse(test = ExpGrps_kmeans_8_result$cluster == 5, yes = "Gene Cluster 5",
                                                                                                             no = ifelse(test = ExpGrps_kmeans_8_result$cluster == 6, yes = "Gene Cluster 6",
                                                                                                                         no = ifelse(test = ExpGrps_kmeans_8_result$cluster == 7, yes = "Gene Cluster 7",
                                                                                                                                     no = "Gene Cluster 8"))))))))

head(ExpGrps_my_gene_col_8clusters)


# Order the data frame rows according to k-means cluster assignments
ExpGrps_my_gene_col_8clusters <- ExpGrps_my_gene_col_8clusters[order(ExpGrps_my_gene_col_8clusters$Cluster), ]
merged_exp_grps_df_ordered <- merged_exp_grps_df[ExpGrps_my_gene_col_8clusters$RowName, ]
rownames(ExpGrps_my_gene_col_8clusters) <- ExpGrps_my_gene_col_8clusters$RowName
ExpGrps_my_gene_col_8clusters$RowName <- NULL
# Create the heatmap
print("about to make heatmap")
HM_merged_exp_grps_df_ordered_8clusters <- pheatmap(merged_exp_grps_df_ordered, annotation_row = ExpGrps_my_gene_col_8clusters, 
                                                    cluster_rows = FALSE, cluster_cols = FALSE, 
                                                    cutree_rows = 8, fontsize_row = .5, 
                                                    main = paste0("K-means Group Gene Expression Matrix (8 clusters) - ",model_name), 
                                                    treeheight_row = 30, angle_col = 315, fontsize_col = 7, 
                                                    width = 10, height = 12, gaps_col = gaps_indices2)
ggsave(filename = paste0(output_directory,"kmeans_8clusters_Group_Gene_Expression_Matrix_",model_name,".pdf"),plot = HM_merged_exp_grps_df_ordered_8clusters, width = 10, height = 12)
write.csv(ExpGrps_my_gene_col_8clusters,
          file = paste0(output_directory, "Group_Kmeans_8_Gene_Clusters_",model_name, ".csv"),
          row.names = TRUE)
