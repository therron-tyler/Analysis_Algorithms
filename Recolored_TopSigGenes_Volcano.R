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
  pos_logchange_sig_genes <<- DEseq_input %>% filter(log2FoldChange > 1) %>% filter(padj < 0.05)
  pos_sorted_genes <- arrange(pos_logchange_sig_genes, padj)
  pos_selected_genes <- head(pos_sorted_genes,10)
  neg_logchange_sig_genes <<- DEseq_input %>% filter(log2FoldChange < -1) %>% filter(padj < 0.05)
  neg_sorted_genes <- arrange(neg_logchange_sig_genes, padj)
  neg_selected_genes <- head(neg_sorted_genes,10)
  selected_sig_genes <- rbind(pos_selected_genes, neg_selected_genes)
  
# The data is cleaned for the plot and labels are now created for each DEseq comparison ran through the function
  output_file_spaces <- gsub("_", " ", output_file)
  reference_condition_spaces <- gsub("_", " ", reference_condition)
  down_genes_custom_label <- as.character(paste0('Down-regulated \nGenes in reference to\n', reference_condition_spaces,'\n(n = ',nrow(neg_logchange_sig_genes),')'))
  up_genes_custom_label <- as.character(paste0('Up-regulated \nGenes in reference to\n', reference_condition_spaces,'\n(n = ',nrow(pos_logchange_sig_genes),')'))
# Define the keyvalue pairs for the data coloring since the original is no longer used
  keyval_data <- DEseq_input
  keyvals <- ifelse(
    keyval_data$padj < 0.05 & keyval_data$log2FoldChange < -1, 'royalblue',
    ifelse(keyval_data$padj < 0.05 & keyval_data$log2FoldChange > 1, 'red',
           'grey50'))
  keyvals[is.na(keyvals)] <- 'grey50'
  names(keyvals)[keyvals == 'royalblue'] <- down_genes_custom_label
  names(keyvals)[keyvals == 'grey50'] <- 'Not \nSignificant'
  names(keyvals)[keyvals == 'red'] <- up_genes_custom_label

# Set the row names as the Gene Symbol
  rownames(DEseq_input) <- DEseq_input$Gene_Symbol
  rownames(selected_sig_genes) <- selected_sig_genes$Gene_Symbol
  rownames(keyval_data) <- DEseq_input$Gene_Symbol # matches the 'selectLab' arg to the 'lab' arg 
  
  plot <- EnhancedVolcano(keyval_data,
                          lab = rownames(keyval_data),
                          selectLab = rownames(selected_sig_genes),
                          x = 'log2FoldChange',
                          y = 'pvalue',
                          title = output_file_spaces, # name of DEseq comparison
                          subtitle = paste0("Reference Condition: ", reference_condition_spaces),
                          titleLabSize = 10,
                          pointSize = 1.5,
                          labSize = 3,
                          drawConnectors = TRUE,
                          cutoffLineType = 'blank',
                          pCutoff = .05,
                          pCutoffCol = 'padj',
                          max.overlaps = Inf,
                          colAlpha = .9,
                          colCustom = keyvals,
                          legendPosition = 'bottom',
                          legendLabSize = 8,
                          subtitleLabSize = 8,
                          labCol = 'black')
  #return(plot)


  ggsave(paste0(outs_directory, output_file, "_SigGenes_Volcano_Plot.pdf"), plot, width = 22, height = 15, units = "cm")
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
  # Count the number of significantly up and down-regulated genes
  #num_upregulated <- nrow(pos_logchange_sig_genes)
  #num_downregulated <- nrow(neg_logchange_sig_genes)
  #num_upregulated
  #num_downregulated
  # Add the annotations to the plot and save it
  #final_plot <- volcano_plot +
  #  annotate("text", x = -Inf, y = Inf, label = paste("Down-regulated genes:", num_downregulated), hjust = 0, vjust = 1, size = 4, color = "royalblue") +
  #  annotate("text", x = Inf, y = Inf, label = paste("Up-regulated genes:", num_upregulated), hjust = 1, vjust = 1, size = 4, color = "red")
  # Save the final plot to a file
  #ggsave(paste0(outs_directory, x$output_file, "SigGenes_Volcano_Plot.pdf"), final_plot, width = 22, height = 15, units = "cm")
})
