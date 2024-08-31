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
# INPUT for script # =============================================
args <- commandArgs(TRUE) # allows script to run in CmdLine

DEseq_BinaryScatter_output_CSV <- readr::read_csv(as.character(args[1])) # data used for the heatmap
FPKM <- readr::read_csv(as.character(args[2])) # genes in the same order as the kmeans heatmap outputs
OUT <- as.character(args[3]) # where to put the output heatmap figs
SAMPLES <- readr::read_csv(as.character(args[4])) # sets the x-axis order on the PDF
Experiment_name <- as.character(args[5])
# =============================================
# functions for script

min_max_normalization <- function(x, min_value = -1, max_value = 1) {
  min_x <- min(x)
  max_x <- max(x)
  normalized_x <- (x - min_x) / (max_x - min_x) * (max_value - min_value) + min_value
  return(normalized_x)
}

update_annotation <- function(df) {
  df <- df %>%
    mutate(Annotation = case_when(
      Annotation == 'Control Only' & Control_LFC > 0 ~ 'Up-regulated Control Only',
      Annotation == 'Control Only' & Control_LFC < 0 ~ 'Down-regulated Control Only',
      Annotation == 'Experimental Only' & Experimental_LFC > 0 ~ 'Up-regulated Experimental Only',
      Annotation == 'Experimental Only' & Experimental_LFC < 0 ~ 'Down-regulated Experimental Only',
      Annotation == 'Shared' & Control_LFC > 0 ~ 'Up-regulated Shared',
      Annotation == 'Shared' & Control_LFC < 0 ~ 'Down-regulated Shared',
      TRUE ~ Annotation
    ))
  return(df)
}
# =============================================

# ============================================= read in code # =============================================
#annotated_master_DEG_df <- readr::read_csv("/Users/ttm3567/Documents/February2024/NPSLE_DoubleScat_QuantCut_BinCut_20240213/double_scatter/CreCom_180508/CSV_CReCOM_180508_Control-CD11cC8_High-Low.csv")
annotated_master_DEG_df <- DEseq_BinaryScatter_output_CSV
  
# this is the FPKM table with the samples and data for visualization later
#FPKM_table <- readr::read_csv("/Users/ttm3567/Documents/February2024/NPSLE_DoubleScat_QuantCut_BinCut_20240213/double_scatter/Input_CSVs_DESeq_Info/CReCOM_HQ_FPKM_Data.csv")
FPKM_table <- FPKM
  
# samples to include (ONLY these samples will be in heatmap x-axis) and order by which they should appear - sample_list
#sample_list <- data.frame(readr::read_csv("/Users/ttm3567/Documents/February2024/NPSLE_DoubleScat_QuantCut_BinCut_20240213/double_scatter/Input_CSVs_DESeq_Info/crecom_180508_HQ_samples.csv"))
sample_list <- SAMPLES
print(FPKM_table[1:10,])
print(annotated_master_DEG_df[1:10,])
print(sample_list[1:10,])

#Experiment_name <- "CReCOM_180508_Control-CD11cC8"
#OUT <- "/Users/ttm3567/Documents/February2024/NPSLE_DoubleScat_QuantCut_BinCut_20240213/binary_HM/CReCOM_180508_CD11cC8_6_AnnotationGrps/"
model_name <- Experiment_name
output_directory <- OUT
# ============================================= read in code # =============================================

# ============================================= processing and plot code # =============================================

filtered_FPKM_counts <- FPKM_table[FPKM_table$Gene_Symbol %in% annotated_master_DEG_df$Genes,] 
print("filtered FPKM counts")

unique_FPKM_counts <- filtered_FPKM_counts %>% distinct(Gene_Symbol, .keep_all = TRUE)
print("removed duplicated FPKM counts")

# the genes between are the same between the two tables now

# Prepare the column names to include 'Symbol' and 'Gene_Symbol'
selected_columns <- c("Symbol", "Gene_Symbol", sample_list$hq_samples)

# Filter and reorder the columns based on selected_columns
# Note: Using `any_of()` to avoid errors in case some sample names in sample_list$hq_samples are not column names in unique_FPKM_counts_df
filtered_FPKM_counts_df <- unique_FPKM_counts %>%
  select(any_of(selected_columns))
print("HQ Sample Names and FPKM counts columns match now")

# create extra annotation groups before the row annotation table is made =================

annotated_master_DEG_df <- update_annotation(annotated_master_DEG_df)
print(annotated_master_DEG_df[1:5,])

# next, match annotations
#unique_FPKM_counts$Annotation <- annotated_master_DEG_df$Annotation

annotated_master_DEG_df <- annotated_master_DEG_df %>% 
  rename_with(~"Gene_Symbol", starts_with("Genes"))
print("annotated FPKM counts")

colnames(annotated_master_DEG_df)



Full_FPKM_annotation_df <- full_join(x = filtered_FPKM_counts_df,
                                     y = annotated_master_DEG_df,
                                     by = "Gene_Symbol")
print("full joined the annotated dataframe and FPKM counts dataframe")


Full_FPKM_order_by_annotation <- Full_FPKM_annotation_df %>%
  arrange(Annotation)
print("sorted joined table by gene group annotation")


Full_FPKM_order_by_annotation <- data.frame(Full_FPKM_order_by_annotation)
print("FPKM counts converted to data.frame")

rownames(Full_FPKM_order_by_annotation) <- Full_FPKM_order_by_annotation$Gene_Symbol
print("FPKM counts have Gene Symbols as row names")

print(Full_FPKM_order_by_annotation[1:5,])

row_annotation_df <- data.frame(Gene_Groups = Full_FPKM_order_by_annotation$Annotation,
                                row.names = rownames(Full_FPKM_order_by_annotation))
print("FPKM counts have a row annotation dataframe")

#fpkm_counts_from_Full_df <- Full_FPKM_order_by_annotation[,3:22] - using indices is bad practice and does not work when datatables change
fpkm_counts_from_Full_celltype <- Full_FPKM_order_by_annotation[, match(sample_list$hq_samples, colnames(Full_FPKM_order_by_annotation))]

print(fpkm_counts_from_Full_celltype[1:10,])
print("FPKM counts no longer have any annotation info in the dataframe")

normalized_fpkm_counts_from_Full <- t(apply(fpkm_counts_from_Full_celltype, 1, min_max_normalization, min_value = -1, max_value = 1))
print("FPKM counts are min-maxed normalized")


# this heatmap is grouped by low and high cell type samples 

HeatMap_binary_celltype <- pheatmap(normalized_fpkm_counts_from_Full, 
                           annotation_row = row_annotation_df, 
                           cluster_cols = FALSE,
                           cluster_rows = FALSE,
                           fontsize_row = .5, 
                           main = paste0("Min-Max Binary Gene Assignments - ",model_name), 
                           treeheight_row = 30, 
                           angle_col = 315, 
                           fontsize_col = 7, 
                           width = 10, height = 12)

print("binary heatmap created, grouped by celltype")


ggsave(filename = paste0(output_directory,"Binary_Celltype_Heatmap_",model_name,".pdf"),
          plot = HeatMap_binary_celltype, width = 8, height = 12)
print("saved")

write_csv(Full_FPKM_order_by_annotation, 
          paste0(output_directory,"Binary_FPKM_Annotation_Data_",model_name,".csv"))
print("wrote CSV")


# ============================================= processing and plot code # =============================================
