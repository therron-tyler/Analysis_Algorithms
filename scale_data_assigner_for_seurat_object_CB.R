
# This is only for Seurat objects that have RNA and ADT assay slots filled

#LIBRARY PATH
libs <- .libPaths("/home/ttm3567/63_tylert/DEseq_libpath_dir/")
#install.packages("Seurat")
library(Seurat)
library(dplyr)
library(tidyr)
library(readr)
# INPUT for script
args <- commandArgs(TRUE) # allows script to run in CmdLine

seurat_object_final <- as.character(args[1]) # INPUT Seurat S4 .rds object locale
desired_output_name <- as.character(args[2]) # NAME of Seurat object that is output
default_ident <- as.character(args[3]) # HAS to be a FACTOR, will be loaded identity CB
output_directory <- as.character(args[4]) # OUTPUT locale - include the '/' at the end because of how code writes the output 
ADT_or_RNA <- as.character(args[5])


# ===================================================================================================
# STEPS that sould be COMPLETED before using this script

# 1. remove desired metadata columns
#p13_merged_20230722@meta.data[["seurat_clusters"]] <- NULL

# 2.  rename any metadata columns as desired
#colnames(tmrc2@meta.data)[colnames(tmrc2@meta.data) == "predicted.id"] <- "Published Annotations"

# ===================================================================================================

# 10/31/23 - Combine ADT and RNA assay features with NCM data ----------------------------------------------

seurat_object_readIN <- readRDS(seurat_object_final)
print("read in data")


seurat_obj <- seurat_object_readIN

# Convert the counts sparse matrix to a dense matrix
dense_counts <- as.matrix(seurat_obj@assays[["RNA"]]@counts)

# Assign the dense matrix to the scale.data slot
seurat_obj@assays[["RNA"]]@scale.data <- dense_counts

# Verify the assignment
dim(seurat_obj@assays[["RNA"]]@scale.data)

# 12. Save the seurat object as RDS
saveRDS(seurat_obj, paste0(output_directory,desired_output_name,".rds"))
print("successfully saved")
