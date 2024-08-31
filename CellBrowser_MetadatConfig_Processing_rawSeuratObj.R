#LIBRARY PATH
libs <- .libPaths("/home/ttm3567/63_tylert/Analysis_Algorithms/Seurat_CellBrowser_Renv_ONLY/")

# libraries to call on in script
library(Seurat)
library(SeuratObject)
library(dplyr)
library(readr)
# allows for command line arguments
args <- commandArgs(TRUE)

# this is where input file path goes and is casted as a string 
rds_path <- as.character(args[1]) # INPUT Seurat S4 .rds object locale
out_path <- as.character(args[2]) # include ending '/'
updated_sObj_name <- as.character(args[3])
meta_column_csv <- read_csv(as.character(args[4]))
ADT_or_RNA <- as.character(args[5])
Idents_MetaColumn <- as.character(args[6])
# testing inputs
# rds_path <- "/Users/ttm3567/Documents/December2023/stromalData_10Mt_subset_cc_regressed.rds"
# meta_column_csv <- read_csv(as.character("/Users/ttm3567/Documents/January2024/CellBrowser_sObj_meta_CSV.csv"))
# 
# updated_sObj_name <- "testingSTROMdata"
# out_path <- "/Users/ttm3567/Documents/January2024/"
# View(meta_column_csv)
# sObj
sObj <- readRDS(rds_path)
print("loaded the Seurat Object")

# create a list of columns to delete
columns_to_remove <- meta_column_csv$delete
for (col in columns_to_remove) {
  sObj@meta.data[[col]] <- NULL
}

replace_old_with_new_Metacol <- function(sObj, meta_column_csv) {
  for (i in 1:nrow(meta_column_csv)) {
    old_name <- as.character(meta_column_csv$old_name[i])
    new_name <- as.character(meta_column_csv$new_name[i])
    
    if (old_name %in% colnames(sObj@meta.data)) {
      colnames(sObj@meta.data)[colnames(sObj@meta.data) == old_name] <- new_name
    } else {
      warning(paste("Column", old_name, "not found in Seurat object's metadata, or is an NA value."))
    }
  }
  return(sObj)
}

# Apply the function to the Seurat object
sObj <- replace_old_with_new_Metacol(sObj, meta_column_csv)

print("old names have been swapped with new")

# Add function that automatically adds a color metadata column for the published datasets if specified in args
# hexidecimal color codes and metadata column as input
#Idents_metadata_column <- meta_column_csv$Colored_Metadata_Column
#Idents_name <- as.character(Idents_metadata_column[1])
Idents_name <- Idents_MetaColumn
Idents_name

# Verify Idents_name is in the metadata
if (!Idents_name %in% colnames(sObj@meta.data)) {
  stop(paste("The Idents column", Idents_name, "is not present in Seurat object metadata."))
}

# Creation of a named color vector which will contain the values which match the manuscript figure
Color_vectors <- meta_column_csv$Color_Vector
Colors <- meta_column_csv$Color

color_df <- na.omit(data.frame(Color_Vector = Color_vectors, Color = Colors))

# Convert to a named vector
color_vector <- setNames(color_df$Color, color_df$Color_Vector)

print("color vector is created")

# take the Idents Name and assign it in the seurat object
Idents(sObj) <- Idents_name
Idents(sObj)
#Idents(sObj) <- sObj[[as.character(Idents_name)]]

print("Specified Idents name is now the active identity")
#print(table(Idents(sObj)))
# Assign colors to cells based on their metadata

missing_colors <- setdiff(Idents(sObj), names(color_vector))
if(length(missing_colors) > 0) {
  warning("The following Idents are missing in color_vector: ", paste(missing_colors, collapse = ", "))
}

sObj[['Published_Annotation_Colors']] <- sapply(sObj@meta.data[[Idents_name]], function(x) color_vector[as.character(x)])
#View(sObj)
print("Colors added to published metadata column")
# check presence of SCT and if so then remove it because it slows everything down
sObj@active.assay <- "RNA"
print("RNA specified as the active assay")

if("SCT" %in% names(sObj@assays)) {
  sObj[["SCT"]] <- NULL
  print("SCT assay removed.")
} else {
  print("No SCT assay present in the object.")
}
### if SCT transform was done - remove it ####
# Not the purpose of this process, the RNA slot is already having a log transformation performed

# copy seurat object in case things go wrong
Seurat_Obj <- sObj

print("data about to be log transformed and hace raw counts replaced")

# perform log transformation on data and merge ADT protein expression with the RNA data 
logT_ADTmerger <- function(seurat_object, ADT_or_RNA) {
  # Normalize RNA assay
  Seurat_Obj <- NormalizeData(Seurat_Obj, assay = "RNA", normalization.method = "LogNormalize", scale.factor = 10000)
  
  # Get normalized RNA data
  rna_normalized <- GetAssayData(Seurat_Obj, assay = "RNA", slot = "data")
  
  if (ADT_or_RNA == "ADT") {
  # 5.  Grab ADT assay
    adt_assay <- GetAssayData(Seurat_Obj, assay = "ADT") # keep
  
  # 6. Pull name of ADT features
    adt_features <- rownames(adt_assay) # keep
  
  # 7. Append "adt_" to the beginning of all the ADT feature names
    adt_features <- paste0("ADT_", adt_features) # keep
  
  # 8. Replace ADT assay feature names with edited feature names
    rownames(adt_assay) <- adt_features # keep
  
  # 9. Combine ADT assay with RNA assay
 # rna_assay <- GetAssayData(Seurat_Obj, assay = "RNA")
    rna_assay <- rna_normalized
    combined_assay <- rbind(rna_assay, adt_assay)
  
  # 10. Replace RNA assay --> counts
    Seurat_Obj@assays$RNA@counts <- combined_assay
  } else {
    print("only RNA assay so not combining any ADT counts, but still replacing raw counts")
    Seurat_Obj@assays$RNA@counts <- rna_normalized
    
  }
  # 11. Inspect results
  tail(rownames(data.frame(Seurat_Obj@assays$RNA@counts)), n = 200)
}

# run function on the seurat object
sObj_merged_transformed <- logT_ADTmerger(Seurat_Obj, ADT_or_RNA)

print("data has been log transformed and if ADT was present then it is merged and replaced the raw counts matrix")

# save the updated object
saveRDS(sObj_merged_transformed,paste0(out_path,updated_sObj_name,".rds"))

