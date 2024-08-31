# USAGE: 1. module load R/4.2.0
# 2. Rscript Relabeled_Volcano_Plot_CmdLine.R [1] CSV file containing DEseq File Paths, Reference condition, and Output File Name [2] Desired output directory location 

libs <- .libPaths("/home/ttm3567/63_tylert/DEseq_libpath_dir/")

library(EnhancedVolcano)
library(DESeq2)
library(purrr)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(tidyr)

## Read command line arguments ------------------------
args <- commandArgs(TRUE)
DEseq_file_list_csv <- as.character(args[1])# include '/' at the beginning and end of the file path
out_directory <- as.character(args[2]) # include '/' at the beginning and end of the file path

outs_directory <- out_directory
# Define the function
create_volcano_plot <- function(input_file_cmd, reference_condition, output_file) {
  DEseq_input <- read.csv(input_file_cmd)
  comparisons2 <- sub("^log2FoldChange.*","log2FoldChange",colnames(DEseq_input))
  DEseq_input <- set_names(DEseq_input, nm = comparisons2)
  # Remove duplicate rows
  DEseq_input <- DEseq_input %>% distinct(Gene_Symbol, .keep_all = TRUE)
  pos_logchange_sig_genes <- DEseq_input %>% filter(log2FoldChange > 1) %>% filter(padj < 0.05)
  pos_sorted_genes <- arrange(pos_logchange_sig_genes, padj)
  pos_selected_genes <- head(pos_sorted_genes,10)
  neg_logchange_sig_genes <- DEseq_input %>% filter(log2FoldChange < -1) %>% filter(padj < 0.05)
  neg_sorted_genes <- arrange(neg_logchange_sig_genes, padj)
  neg_selected_genes <- head(neg_sorted_genes,10)
  selected_sig_genes <- rbind(pos_selected_genes, neg_selected_genes)
  # Set the row names as the Gene Symbol
  rownames(DEseq_input) <- DEseq_input$Gene_Symbol
  rownames(selected_sig_genes) <- selected_sig_genes$Gene_Symbol
  
  plot <- EnhancedVolcano(DEseq_input,
                          lab = rownames(DEseq_input),
                          selectLab = rownames(selected_sig_genes),
                          x = 'log2FoldChange',
                          y = 'pvalue',
                          title = output_file, # name of DEseq comparison
                          subtitle = paste0("Reference Condition: ", reference_condition),
                          titleLabSize = 10,
                          pointSize = 1.5,
                          labSize = 3,
                          drawConnectors = TRUE,
                          cutoffLineType = 'blank',
                          pCutoff = .05,
                          pCutoffCol = 'padj',
                          legendLabels = c("NS", expression(Log[2] ~ FC), "adjusted p-value", expression(adjusted ~ p - value ~ and ~ log[2] ~ FC)),
                          max.overlaps = Inf)
  
  ggsave(paste0(outs_directory, output_file, "_Volcano_Plot.pdf"), plot, width = 22, height = 15, units = "cm")
}

# read in a CSV that contains multiple DEseq files
input_volcano_CSV <- read.csv(DEseq_file_list_csv)

# make three lists so each once can be iterated through

fp_list <- as.list(input_volcano_CSV$file_path)
ref_list <- as.list(input_volcano_CSV$reference_condition)
output_list <- as.list(input_volcano_CSV$output_file)

# List of input files - turn into an LAPPLY
input_file_cmd <- lapply(seq_along(fp_list), function(x) {
  list(file_path = as.character(fp_list[[x]]), reference_condition = as.character(ref_list[[x]]), output_file = as.character(output_list[[x]]))
})


# Apply the function to the list of input files
lapply(input_file_cmd, function(x) {
  create_volcano_plot(x$file_path, x$reference_condition, x$output_file)
})
