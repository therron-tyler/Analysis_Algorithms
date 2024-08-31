
#LIBRARY PATH
libs <- .libPaths("/home/ttm3567/63_tylert/Analysis_Algorithms/Seurat_CellBrowser_Renv_ONLY/")
#install.packages("Seurat")
library(Seurat)
library(dplyr)
library(readr)
# INPUT for script
args <- commandArgs(TRUE) # allows script to run in CmdLine

# the filtered BC matrix folder within the cell ranger output folder
Ten_x_directory_1 <- as.character(args[1]) # "/path/gto/Sample1"
Ten_x_directory_2 <- as.character(args[2]) # "/path/gto/Sample2"
Project_name_1 <- as.character(args[3]) # "Sample1"
Project_name_2 <- as.character(args[4]) # "Sample2"
Merged_Project_name <- as.character(args[5]) # "Merged_Sample1_Sample2"
Output_directory <- as.character(args[6]) # "/path/gto/save/"

# read the 10x directory
data.sample1 <- Read10X(data.dir = Ten_x_directory_1)
data.sample2 <- Read10X(data.dir = Ten_x_directory_2)
print("read in 10x directories")

# create object from the 10x directory
seurat.sample1 <- CreateSeuratObject(counts = data.sample1, project = Project_name_1)
seurat.sample2 <- CreateSeuratObject(counts = data.sample2, project = Project_name_2)
print("create seurat object")

# what to name the metadata columns
seurat.sample1$sample <- Project_name_1
seurat.sample2$sample <- Project_name_2
print("annotated cell origins before merging")

# merge the two objects
seurat.merged <- merge(seurat.sample1, y = list(seurat.sample2), add.cell.ids = c(Project_name_1, Project_name_2), project = Merged_Project_name)
print("merged objects")

# save the raw RDS object
saveRDS(seurat.merged, file = paste0(Output_directory,Merged_Project_name,"_raw.rds"))
print("save raw RDS")

# process the RDS object ==================

# Normalize data
seurat.merged <- NormalizeData(seurat.merged, normalization.method = "LogNormalize", scale.factor = 10000)
print("log transform - factor is 10k")
# Find variable features
seurat.merged <- FindVariableFeatures(seurat.merged, selection.method = "vst", nfeatures = 2000)
print("var features - 2000")

# Scale data
seurat.merged <- ScaleData(seurat.merged)
print("scale")

# Run PCA
seurat.merged <- RunPCA(seurat.merged, features = VariableFeatures(object = seurat.merged))
print("pca")

# Find neighbors and clusters
seurat.merged <- FindNeighbors(seurat.merged, dims = 1:10)
print("neighbors")

seurat.merged <- FindClusters(seurat.merged, resolution = 0.5)
print("find clusters")

# Run UMAP
seurat.merged <- RunUMAP(seurat.merged, dims = 1:10)
print("UMAP")

# save the raw RDS object
saveRDS(seurat.merged, file = paste0(Output_directory,Merged_Project_name,"_default_processing.rds"))
print("save processed RDS")

