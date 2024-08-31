# USAGE: 1. module load R/4.2.0
# 2. Rscript Automated_CellBrowser_Cmds_SeuratObj_prep.R </this/is/directory/path/random_seurat_object.rds> <NAME> <seurat_default_metadata_column> </output/directory/> <"ADT"/"RNA"> 
#library(Seurat)
#library(dplyr)
#library(tidyr)
#library(readr)

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
#  ================================= ================================= ================================= ================================= ================================= ================================= ================================= ================================= =================================

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

seurat_object_readIN@active.assay <- 'RNA'
print("that active assay, is forsure RNA now")

# this is to make sure the markers are set for the published annotation
Idents(seurat_object_readIN) <- default_ident
print("default Ident set")

# 3. Attach markers manually to the Seurat object before saving it
seurat_object_readIN@misc$markers <- FindAllMarkers(seurat_object_readIN)
print("attached markers")


Seurat_Obj <- seurat_object_readIN
# --- adding in ADT protein features ---
# need to normalize the data first

Seurat_Obj <- NormalizeData(Seurat_Obj, normalization.method = "LogNormalize", scale.factor = 10000)
  
Seurat_Obj@assays$RNA@counts <- Seurat_Obj@assays[["RNA"]]@data

# 5.  Grab ADT assay
if (ADT_or_RNA == "ADT") {
 adt_assay <- GetAssayData(Seurat_Obj, assay = "ADT")
 print("grabbed ADT assay")

# 6. Pull name of ADT features
 adt_features <- rownames(adt_assay)
 print("grabbed ADT names")

# 7. Append "adt_" to the beginning of all the ADT feature names
 adt_features <- paste0("ADT_", adt_features)
 print("appended _ADT to features")

# 8. Replace ADT assay feature names with edited feature names
 rownames(adt_assay) <- adt_features
 print("replaced feature ADT names")

# 9. Combine ADT assay with RNA assay
 rna_assay <- GetAssayData(Seurat_Obj, assay = "RNA", slot = "counts")
 combined_assay <- rbind(rna_assay, adt_assay)
 print("Combined RNA and ADT assay")
# 10. Replace RNA assay --> counts
 Seurat_Obj@assays$RNA@counts <- combined_assay
 print("replaced counts in RNA assay")
} else {
 print("only RNA assay so not combining any ADT counts")
}
# 11. Inspect results
tail(rownames(data.frame(Seurat_Obj@assays$RNA@counts)), n = 200)

# 4. specify the default Cell Identities -> Idents(object) <- "metadata_column_name" - HAS to be last step and HAS to be a FACTOR
Idents(Seurat_Obj) <- default_ident 
print("default Ident set")
# 12. Save the seurat object as RDS
saveRDS(Seurat_Obj, paste0(output_directory,desired_output_name,".rds"))
print("successfully saved")
