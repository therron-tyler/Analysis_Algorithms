# USAGE: 1. module load R/4.2.0

# 2. Rscript VennDiagram_CmdLine.R [1]path to CSV containing the DEseq file names and the path to their DEseq output [2]output directory

# the CSV input has four columns - file1, file2, path1, and path2
# file1 is the 'control' group name
# file2 is the 'experimental' group name
# path1 goes to the control group DEseq output
# path2 goes to the experimental group DEseq output

libs <- .libPaths("/home/ttm3567/63_tylert/DEseq_libpath_dir/")

library(ggplot2)
library(ggrepel)
library(dplyr)
library(DESeq2)
library(topGO)
library(VennDiagram)
library(gridExtra)                # Load gridExtra package
library(purrr)
## Read command line arguments ------------------------
args <- commandArgs(TRUE)
CSV_input_path <- as.character(args[1])
out_directory <- as.character(args[2]) # include '/' at the beginning and end of the file path

# read in the DEseq files ================================= =================================

outs_directory <- out_directory

CSV_input <- read.csv(CSV_input_path)

# labels for the Venn Diagram
cat_names <- c(list(CSV_input$file1),list(CSV_input$file2))

cat_names_edit <- lapply(cat_names, function(x) {
  sub("*_vs_*","\n_vs_\n", x)
})


# create two category name lists 
control_cat_names_edit <- cat_names_edit[[1]]
exp_cat_names_edit <- cat_names_edit[[2]]
control_cat_names <- cat_names[[1]]
exp_cat_names <- cat_names[[2]]
control_cat_names_edit
exp_cat_names_edit
control_cat_names
exp_cat_names

# read in the control DEseq dataframes ############## Differential UPregulated Genes ###############
controls_list <- list()
for(i in 1:nrow(CSV_input)){
  row <- read.csv(CSV_input$path1[i])
  controls_list[[i]] <- row
}
# change log2FoldChange column so it is the same for everything
controls_list <- lapply(controls_list, function(x) {
  comparisons2 <- sub("^log2FoldChange.*","log2FoldChange",colnames(x))
  set_names(x, nm = comparisons2)
})


# read in experimental DEseq dataframes and change log2FoldChange column names
experimental_list <- list()
for(i in 1:nrow(CSV_input)){
  row <- read.csv(CSV_input$path2[i])
  experimental_list[[i]] <- row
}
# change log2FoldChange column so it is the same for everything
experimental_list <- lapply(experimental_list, function(x) {
  comparisons2 <- sub("^log2FoldChange.*","log2FoldChange",colnames(x))
  set_names(x, nm = comparisons2)
})


# the files are read in and need names 
controls_name_list <- list(CSV_input$file1)
for(i in 1:length(controls_name_list)){
  names(controls_list) <- controls_name_list[[i]]
}

experimental_name_list <- list(CSV_input$file2)
for(i in 1:length(experimental_name_list)){
  names(experimental_list) <- experimental_name_list[[i]]
}

# both conditions have been read in and are ready to be filtered

# Upregulated differential genes #####  control Condition ######
CD11c_ctrl <- lapply(controls_list, function(x) x %>% filter(log2FoldChange > 1))
CD11c_ctrl_pval_filtered <- lapply(CD11c_ctrl, function(x) x %>% filter(padj < 0.05))
CD11c_ctrl_genes_upreg <- lapply(CD11c_ctrl_pval_filtered, function(x) length(x$X)) 
CD11c_ctrl_genes_upreg

lapply(seq_along(CD11c_ctrl_pval_filtered), function(x) {if (CD11c_ctrl_genes_upreg[[x]] > 0) {
  write.csv(CD11c_ctrl_pval_filtered[[x]], 
            file = paste0(outs_directory, control_cat_names[[x]], "__", "significant_UPregulated_Genes",".csv"), 
            row.names = TRUE)
} else {
  print(CD11c_ctrl_genes_upreg[[x]])
}
})

# Upregulated differential genes #####  experimental Condition ######
CD11c_experimental <- lapply(experimental_list, function(x) x %>% filter(log2FoldChange > 1))
CD11c_experimental_pval_filtered <- lapply(CD11c_experimental, function(x) x %>% filter(padj < 0.05))
CD11c_experimental_genes_upreg <- lapply(CD11c_experimental_pval_filtered, function(x) length(x$X)) 
CD11c_experimental_genes_upreg

lapply(seq_along(CD11c_experimental_pval_filtered), function(x) {if (CD11c_experimental_genes_upreg[[x]] > 0) {
  write.csv(CD11c_experimental_pval_filtered[[x]], 
            file = paste0(outs_directory, exp_cat_names[[x]], "__", "significant_UPregulated_Genes",".csv"), 
            row.names = TRUE)
} else {
  print(CD11c_experimental_genes_upreg[[x]])
}
})

# Upregulated differential genes ##### Shared among CD11c ctrl and CD11c experimental ######

# USE seq_along and lapply ### 3/16/2023 ##
shared_upregulated_genes <- lapply(seq_along(CD11c_ctrl_pval_filtered), function(x){
  merge(CD11c_ctrl_pval_filtered[[x]],CD11c_experimental_pval_filtered[[x]], by="X")})

shared_up_diff_genes <- lapply(shared_upregulated_genes, function(x) length(x$X)) 
shared_up_diff_genes




# part II ###########################
venn_diagrams <- lapply(seq_along(CD11c_ctrl_genes_upreg), function(x){
  grid.newpage()
  draw.pairwise.venn(area1=CD11c_ctrl_genes_upreg[[x]], area2=CD11c_experimental_genes_upreg[[x]],cross.area=shared_up_diff_genes[[x]],
                     category=c(control_cat_names_edit[[x]],exp_cat_names_edit[[x]]),fill=c("Red","Yellow"), cat.pos = c(-6,6),fontfamily = rep("serif", 3), cat.cex = rep(.5, 1), cat.fontfamily = rep("serif", 1), cat.dist = rep(.05, 2),euler.d = FALSE, scaled = FALSE,ext.dist = rep(-2, 2))
})

venn_diagrams <- lapply(venn_diagrams, function(x){grid.arrange(gTree(children =x), # Add title & subtitle
                                                                top = "Venn Diagram of Significant Up-regulated Genes",
                                                                bottom = "CD11chi versus CD11clo")})


lapply(seq_along(venn_diagrams), function(x) {ggsave(filename = paste0(outs_directory,control_cat_names[[x]],"__",exp_cat_names[[x]],"_significant_UPregulated_Genes_VennDiagram",".pdf"),plot = venn_diagrams[[x]], width = 10, height = 8)})

# write the csv only if it has shared genes 
lapply(seq_along(venn_diagrams), function(x) {if (shared_up_diff_genes[[x]] > 0) {
  write.csv(shared_upregulated_genes[[x]], 
            file = paste0(outs_directory, control_cat_names[[x]], "__", exp_cat_names[[x]], "_significant_UPregulated_Genes",".csv"), 
            row.names = TRUE)
} else {
  print(shared_up_diff_genes[[x]])
}
})


#### Differential DOWNregulated Genes #################################### #################################### #################################### #################################### #################################### ####################################


# both conditions have been read in at beginning of the script and are ready to be filtered

# DOWNregulated differential genes #####  control Condition ######
CD11c_ctrl <- lapply(controls_list, function(x) x %>% filter(log2FoldChange < -1))
CD11c_ctrl_pval_filtered <- lapply(CD11c_ctrl, function(x) x %>% filter(padj < 0.05))
CD11c_ctrl_genes_Dreg <- lapply(CD11c_ctrl_pval_filtered, function(x) length(x$X)) 
CD11c_ctrl_genes_Dreg

lapply(seq_along(CD11c_ctrl_pval_filtered), function(x) {if (CD11c_ctrl_genes_Dreg[[x]] > 0) {
  write.csv(CD11c_ctrl_pval_filtered[[x]], 
            file = paste0(outs_directory, control_cat_names[[x]], "__", "significant_DOWNregulated_Genes",".csv"), # 3/27/2023 - write down-reg sig genes
            row.names = TRUE)
} else {
  print(CD11c_ctrl_genes_Dreg[[x]])
}
})


# Upregulated differential genes #####  experimental Condition ######
CD11c_experimental <- lapply(experimental_list, function(x) x %>% filter(log2FoldChange < -1))
CD11c_experimental_pval_filtered <- lapply(CD11c_experimental, function(x) x %>% filter(padj < 0.05))
CD11c_experimental_genes_Dreg <- lapply(CD11c_experimental_pval_filtered, function(x) length(x$X)) 
CD11c_experimental_genes_Dreg

lapply(seq_along(CD11c_experimental_pval_filtered), function(x) {if (CD11c_experimental_genes_Dreg[[x]] > 0) {
  write.csv(CD11c_experimental_pval_filtered[[x]], 
            file = paste0(outs_directory, exp_cat_names[[x]], "__", "significant_DOWNregulated_Genes",".csv"), # 3/27/2023 - write down-reg sig genes
            row.names = TRUE)
} else {
  print(CD11c_experimental_genes_Dreg[[x]])
}
})

# Upregulated differential genes ##### Shared among CD11c ctrl and CD11c experimental ######

shared_Dregulated_genes <- lapply(seq_along(CD11c_ctrl_pval_filtered), function(x){
  merge(CD11c_ctrl_pval_filtered[[x]],CD11c_experimental_pval_filtered[[x]], by="X")})

shared_D_diff_genes <- lapply(shared_Dregulated_genes, function(x) length(x$X)) 
shared_D_diff_genes
# labels for the Venn Diagram
cat_names <- c(list(CSV_input$file1),list(CSV_input$file2))

cat_names_edit <- lapply(cat_names, function(x) {
  sub("*_vs_*","\n_vs_\n", x)
})


# part II ###########################
venn_diagrams <- lapply(seq_along(CD11c_ctrl_genes_Dreg), function(x){
  grid.newpage()
  draw.pairwise.venn(area1=CD11c_ctrl_genes_Dreg[[x]], area2=CD11c_experimental_genes_Dreg[[x]],cross.area=shared_D_diff_genes[[x]],
                     category=c(control_cat_names_edit[[x]],exp_cat_names_edit[[x]]),fill=c("Red","Yellow"), cat.pos = c(-6,6),fontfamily = rep("serif", 3), cat.cex = rep(.5, 1), cat.fontfamily = rep("serif", 1), cat.dist = rep(.05, 2),euler.d = FALSE, scaled = FALSE,ext.dist = rep(-2, 2))
})

venn_diagrams <- lapply(venn_diagrams, function(x){grid.arrange(gTree(children =x), # Add title & subtitle
                                                                top = "Venn Diagram of Significant Down-regulated Genes",
                                                                bottom = "CD11chi versus CD11clo")})

# write pdf of shared genes

lapply(seq_along(venn_diagrams), function(x) {ggsave(filename = paste0(outs_directory,control_cat_names[[x]],"__",exp_cat_names[[x]],"_significant_DOWNregulated_Genes_VennDiagram",".pdf"),plot = venn_diagrams[[x]], width = 10, height = 8)})



# write the csv only if it has shared genes 
lapply(seq_along(venn_diagrams), function(x) {if (shared_D_diff_genes[[x]] > 0) {
  write.csv(shared_Dregulated_genes[[x]], 
            file = paste0(outs_directory, control_cat_names[[x]], "__", exp_cat_names[[x]], "_significant_DOWNregulated_Genes",".csv"), 
            row.names = TRUE)
} else {
  print(shared_D_diff_genes[[x]])
}
})


