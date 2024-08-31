# USAGE: 1. module load R/4.2.0

# 2. Rscript VennDiagram_CmdLine.R [1]path to CSV containing the DEseq file names and the path to their DEseq output [2]output directory [3]LFC cutoff value

# the CSV input has four columns - file1, file2, file3, path1, path2, path3
# file1 is the 'control' group name  ----> LFCabove1_<DEseqfile1>
# file2 is the 'experimental1' group name ----> LFCabove1_<DEseqfile2>
# file3 is the 'experimental2' group name ----> LFCbelow-1_<DEseqfile1>
# file4 is the 'experimental2' group name ----> LFCbelow-1_<DEseqfile2>
# path1 goes to the LFC>1 group1 DEseq output ----> path1 and path3 are the same; filtering is different
# path2 goes to the LFC>1 group2 DEseq output ----> path2 and path4 are the same; filtering is different
# path3 goes to the LFC<-1 group1 DEseq output
# path4 goes to the LFC<-1 group2 DEseq output

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
LFC_input <- as.numeric(args[3])
# read in the DEseq files ================================= =================================

outs_directory <- out_directory

CSV_input <- read.csv(CSV_input_path)


LFC <- LFC_input
# labels for the Venn Diagram
cat_names <- c(list(CSV_input$file1),list(CSV_input$file2),list(CSV_input$file3),list(CSV_input$file4))

cat_names_edit <- lapply(cat_names, function(x) {
  sub("*_vs_*","\n_vs_\n", x)
})

# create two category name lists 
control_cat_names_edit <- cat_names_edit[[1]]
exp1_cat_names_edit <- cat_names_edit[[2]]
exp2_cat_names_edit <- cat_names_edit[[3]]
exp3_cat_names_edit <- cat_names_edit[[4]]

control_cat_names <- cat_names[[1]]
exp_cat_names1 <- cat_names[[2]]
exp_cat_names2 <- cat_names[[3]]
exp_cat_names3 <- cat_names[[4]]

control_cat_names_edit
exp1_cat_names_edit
exp2_cat_names_edit
exp3_cat_names_edit

control_cat_names
exp_cat_names1
exp_cat_names2
exp_cat_names3

# change the names for the output files so they are not that lengthy 
file_label_deseq1 <- lapply(control_cat_names, function(x) {
  sub("LFC_>_1_*","", x)})
file_label_deseq2 <- lapply(exp_cat_names1, function(x) {
  sub("LFC_>_1_*","", x)})
file_label_deseq1
file_label_deseq2
file_label_deseq3 <- lapply(exp_cat_names2, function(x) {
  sub("LFC_<_-1_*","", x)})
file_label_deseq4 <- lapply(exp_cat_names3, function(x) {
  sub("LFC_<_-1_*","", x)})
file_label_deseq3
file_label_deseq4

# remove underscores from the venn diagram

# Replace underscores with spaces
control_cat_names_edit <- lapply(control_cat_names_edit, function(x) gsub("_", " ", x))
control_cat_names_edit
exp1_cat_names_edit <- lapply(exp1_cat_names_edit, function(x) gsub("_", " ", x))
exp1_cat_names_edit
exp2_cat_names_edit <- lapply(exp2_cat_names_edit, function(x) gsub("_", " ", x))
exp2_cat_names_edit
exp3_cat_names_edit <- lapply(exp3_cat_names_edit, function(x) gsub("_", " ", x))
exp3_cat_names_edit

# read in the control DEseq dataframes ############## Differential UPregulated Genes ###############
controls_list <- list() # - CONTROL
for(i in 1:nrow(CSV_input)){
  row <- read.csv(CSV_input$path1[i])
  controls_list[[i]] <- row
}
# change log2FoldChange column so it is the same for everything
controls_list <- lapply(controls_list, function(x) {
  comparisons2 <- sub("^log2FoldChange.*","log2FoldChange",colnames(x))
  set_names(x, nm = comparisons2)
})


# read in experimental DEseq dataframes and change log2FoldChange column names - EXP1
experimental_list1 <- list()
for(i in 1:nrow(CSV_input)){
  row <- read.csv(CSV_input$path2[i])
  experimental_list1[[i]] <- row
}
# change log2FoldChange column so it is the same for everything
experimental_list1 <- lapply(experimental_list1, function(x) {
  comparisons2 <- sub("^log2FoldChange.*","log2FoldChange",colnames(x))
  set_names(x, nm = comparisons2)
})

# EXP2
experimental_list2 <- list()
for(i in 1:nrow(CSV_input)){
  row <- read.csv(CSV_input$path3[i])
  experimental_list2[[i]] <- row
}
# change log2FoldChange column so it is the same for everything
experimental_list2 <- lapply(experimental_list2, function(x) {
  comparisons2 <- sub("^log2FoldChange.*","log2FoldChange",colnames(x))
  set_names(x, nm = comparisons2)
})


# EXP3
experimental_list3 <- list()
for(i in 1:nrow(CSV_input)){
  row <- read.csv(CSV_input$path4[i])
  experimental_list3[[i]] <- row
}
# change log2FoldChange column so it is the same for everything
experimental_list3 <- lapply(experimental_list3, function(x) {
  comparisons2 <- sub("^log2FoldChange.*","log2FoldChange",colnames(x))
  set_names(x, nm = comparisons2)
})


# the file groups are read in and need names ++++++++++++++++++++
#CONTROL
controls_name_list <- list(CSV_input$file1)
for(i in 1:length(controls_name_list)){
  names(controls_list) <- controls_name_list[[i]]
}

# EXP1
experimental1_name_list <- list(CSV_input$file2)
for(i in 1:length(experimental1_name_list)){
  names(experimental_list1) <- experimental1_name_list[[i]]
}

#EXP2
experimental2_name_list <- list(CSV_input$file3)
for(i in 1:length(experimental2_name_list)){
  names(experimental_list2) <- experimental2_name_list[[i]]
}

#EXP2
experimental3_name_list <- list(CSV_input$file4)
for(i in 1:length(experimental3_name_list)){
  names(experimental_list3) <- experimental3_name_list[[i]]
}
# both conditions have been read in and are ready to be filtered

# LFC is the variable that holds the cutoff value for filtering, adjusted p-value will stay at 0.05 typically. 

# Upregulated differential genes #####  control Condition ######
CD11c_ctrl <- lapply(controls_list, function(x) x %>% filter(log2FoldChange > LFC))
CD11c_ctrl_pval_filtered <- lapply(CD11c_ctrl, function(x) x %>% filter(padj < 0.05))
CD11c_ctrl_genes_upreg <- lapply(CD11c_ctrl_pval_filtered, function(x) length(x$X)) 
CD11c_ctrl_genes_upreg

lapply(seq_along(CD11c_ctrl_pval_filtered), function(x) {if (CD11c_ctrl_genes_upreg[[x]] > 0) {
  write.csv(CD11c_ctrl_pval_filtered[[x]], 
            file = paste0(outs_directory, "LFC_>_1_Genes_",file_label_deseq1[[x]], ".csv"), 
            row.names = TRUE)
} else {
  print(CD11c_ctrl_genes_upreg[[x]])
}
})

# Upregulated differential genes ##### 1st experimental Condition ######
CD11c_experimental1 <- lapply(experimental_list1, function(x) x %>% filter(log2FoldChange > LFC))
CD11c_experimental_pval_filtered1 <- lapply(CD11c_experimental1, function(x) x %>% filter(padj < 0.05))
CD11c_experimental_genes_upreg1 <- lapply(CD11c_experimental_pval_filtered1, function(x) length(x$X)) 
CD11c_experimental_genes_upreg1

lapply(seq_along(CD11c_experimental_pval_filtered1), function(x) {if (CD11c_experimental_genes_upreg1[[x]] > 0) {
  write.csv(CD11c_experimental_pval_filtered1[[x]], 
            file = paste0(outs_directory, "LFC_>_1_Genes_", file_label_deseq2[[x]],".csv"), 
            row.names = TRUE)
} else {
  print(CD11c_experimental_genes_upreg1[[x]])
}
})


# Upregulated differential genes ##### 2nd experimental Condition ###### EXP2
CD11c_experimental2 <- lapply(experimental_list2, function(x) x %>% filter(log2FoldChange < -LFC))
CD11c_experimental_pval_filtered2 <- lapply(CD11c_experimental2, function(x) x %>% filter(padj < 0.05))
CD11c_experimental_genes_upreg2 <- lapply(CD11c_experimental_pval_filtered2, function(x) length(x$X)) 
CD11c_experimental_genes_upreg2

lapply(seq_along(CD11c_experimental_pval_filtered2), function(x) {if (CD11c_experimental_genes_upreg2[[x]] > 0) {
  write.csv(CD11c_experimental_pval_filtered2[[x]], 
            file = paste0(outs_directory, "LFC_<_-1_Genes_",file_label_deseq3[[x]], ".csv"), 
            row.names = TRUE)
} else {
  print(CD11c_experimental_genes_upreg2[[x]])
}
})


# Upregulated differential genes ##### 3rd experimental Condition ###### EXP3
CD11c_experimental3 <- lapply(experimental_list3, function(x) x %>% filter(log2FoldChange < -LFC))
CD11c_experimental_pval_filtered3 <- lapply(CD11c_experimental3, function(x) x %>% filter(padj < 0.05))
CD11c_experimental_genes_upreg3 <- lapply(CD11c_experimental_pval_filtered3, function(x) length(x$X)) 
CD11c_experimental_genes_upreg3

lapply(seq_along(CD11c_experimental_pval_filtered3), function(x) {if (CD11c_experimental_genes_upreg3[[x]] > 0) {
  write.csv(CD11c_experimental_pval_filtered3[[x]], 
            file = paste0(outs_directory, "LFC_<_-1_Genes_",file_label_deseq4[[x]], ".csv"), 
            row.names = TRUE)
} else {
  print(CD11c_experimental_genes_upreg3[[x]])
}
})
# Upregulated differential genes ##### Shared among CD11c ctrl and CD11c experimental ######

# merge the three data frames by the "X" column (Ensemble ID)
shared_upregulated_genes <- lapply(seq_along(CD11c_ctrl_pval_filtered), function(x){merge(merge(merge(CD11c_ctrl_pval_filtered[[x]], CD11c_experimental_pval_filtered1[[x]], by = "X"), CD11c_experimental_pval_filtered2[[x]], by = "X"), CD11c_experimental_pval_filtered3[[x]], by = "X")
})
shared_up_diff_genes <- lapply(shared_upregulated_genes, function(x) length(x$X)) 
shared_up_diff_genes

# merge area1 and area2
n12_upregulated_genes <- lapply(seq_along(CD11c_ctrl_pval_filtered), function(x){merge(CD11c_ctrl_pval_filtered[[x]], CD11c_experimental_pval_filtered1[[x]], by = "X")
})
n12_up_diff_genes <- lapply(n12_upregulated_genes, function(x) length(x$X)) 
n12_up_diff_genes
lapply(seq_along(n12_upregulated_genes), function(x) {if (n12_up_diff_genes[[x]] > 0) {
  write.csv(n12_upregulated_genes[[x]], 
            file = paste0(outs_directory,"Shared_Genes1_", file_label_deseq1[[x]],"__",file_label_deseq2[[x]],".csv"), 
            row.names = TRUE)
} else {
  print(n12_up_diff_genes[[x]])
}
})

# merge area2 and area3
n23_upregulated_genes <- lapply(seq_along(CD11c_ctrl_pval_filtered), function(x){merge(CD11c_experimental_pval_filtered1[[x]], CD11c_experimental_pval_filtered2[[x]], by = "X")
})
n23_up_diff_genes <- lapply(n23_upregulated_genes, function(x) length(x$X)) 
n23_up_diff_genes
lapply(seq_along(n23_upregulated_genes), function(x) {if (n23_up_diff_genes[[x]] > 0) {
  write.csv(n23_upregulated_genes[[x]], 
            file = paste0(outs_directory, "Shared_Genes2_", file_label_deseq2[[x]],"__",file_label_deseq3[[x]],".csv"), 
            row.names = TRUE)
} else {
  print(n23_up_diff_genes[[x]])
}
})

# merge area3 and area1
n13_upregulated_genes <- lapply(seq_along(CD11c_ctrl_pval_filtered), function(x){merge(CD11c_ctrl_pval_filtered[[x]], CD11c_experimental_pval_filtered2[[x]], by = "X")
})
n13_up_diff_genes <- lapply(n13_upregulated_genes, function(x) length(x$X)) 
n13_up_diff_genes
lapply(seq_along(n13_upregulated_genes), function(x) {if (n13_up_diff_genes[[x]] > 0) {
  write.csv(n13_upregulated_genes[[x]], 
            file = paste0(outs_directory,"Shared_Genes3_",file_label_deseq1[[x]],"__",file_label_deseq3[[x]],".csv"), 
            row.names = TRUE)
} else {
  print(n13_up_diff_genes[[x]])
}
})

# merge area1 and area4
n14_upregulated_genes <- lapply(seq_along(CD11c_ctrl_pval_filtered), function(x){merge(CD11c_ctrl_pval_filtered[[x]], CD11c_experimental_pval_filtered3[[x]], by = "X")
})
n14_up_diff_genes <- lapply(n14_upregulated_genes, function(x) length(x$X)) 
n14_up_diff_genes
lapply(seq_along(n14_upregulated_genes), function(x) {if (n14_up_diff_genes[[x]] > 0) {
  write.csv(n14_upregulated_genes[[x]], 
            file = paste0(outs_directory,"Shared_Genes4_", file_label_deseq1[[x]],"__",file_label_deseq4[[x]],".csv"), 
            row.names = TRUE)
} else {
  print(n14_up_diff_genes[[x]])
}
})

# merge area2 and area4
n24_upregulated_genes <- lapply(seq_along(CD11c_ctrl_pval_filtered), function(x){merge(CD11c_experimental_pval_filtered1[[x]], CD11c_experimental_pval_filtered3[[x]], by = "X")
})
n24_up_diff_genes <- lapply(n24_upregulated_genes, function(x) length(x$X)) 
n24_up_diff_genes
lapply(seq_along(n24_upregulated_genes), function(x) {if (n24_up_diff_genes[[x]] > 0) {
  write.csv(n24_upregulated_genes[[x]], 
            file = paste0(outs_directory,"Shared_Genes5_", file_label_deseq2[[x]],"__",file_label_deseq4[[x]],".csv"), 
            row.names = TRUE)
} else {
  print(n24_up_diff_genes[[x]])
}
})

# merge area3 and area4
n34_upregulated_genes <- lapply(seq_along(CD11c_ctrl_pval_filtered), function(x){merge(CD11c_experimental_pval_filtered2[[x]], CD11c_experimental_pval_filtered3[[x]], by = "X")
})
n34_up_diff_genes <- lapply(n34_upregulated_genes, function(x) length(x$X)) 
n34_up_diff_genes
lapply(seq_along(n34_upregulated_genes), function(x) {if (n34_up_diff_genes[[x]] > 0) {
  write.csv(n34_upregulated_genes[[x]], 
            file = paste0(outs_directory,"Shared_Genes6_", file_label_deseq3[[x]],"__",file_label_deseq4[[x]],".csv"), 
            row.names = TRUE)
} else {
  print(n34_up_diff_genes[[x]])
}
})

# merge area1, area2, and area4
n124_upregulated_genes <- lapply(seq_along(CD11c_ctrl_pval_filtered), function(x){merge(merge(CD11c_ctrl_pval_filtered[[x]], CD11c_experimental_pval_filtered1[[x]], by = "X"), CD11c_experimental_pval_filtered3[[x]], by = "X")
})
n124_up_diff_genes <- lapply(n124_upregulated_genes, function(x) length(x$X)) 
n124_up_diff_genes
lapply(seq_along(n124_upregulated_genes), function(x) {if (n124_up_diff_genes[[x]] > 0) {
  write.csv(n124_upregulated_genes[[x]], 
            file = paste0(outs_directory,"Shared_Genes7_", file_label_deseq1[[x]],"__",file_label_deseq2[[x]],"__",file_label_deseq4[[x]],".csv"), 
            row.names = TRUE)
} else {
  print(n124_up_diff_genes[[x]])
}
})

# merge area1, area2, and area3
n123_upregulated_genes <- lapply(seq_along(CD11c_ctrl_pval_filtered), function(x){merge(merge(CD11c_ctrl_pval_filtered[[x]], CD11c_experimental_pval_filtered1[[x]], by = "X"), CD11c_experimental_pval_filtered2[[x]], by = "X")
})
n123_up_diff_genes <- lapply(n123_upregulated_genes, function(x) length(x$X)) 
n123_up_diff_genes
lapply(seq_along(n123_upregulated_genes), function(x) {if (n123_up_diff_genes[[x]] > 0) {
  write.csv(n123_upregulated_genes[[x]], 
            file = paste0(outs_directory,"Shared_Genes8_", file_label_deseq1[[x]],"__",file_label_deseq2[[x]],"__",file_label_deseq3[[x]],".csv"), 
            row.names = TRUE)
} else {
  print(n123_up_diff_genes[[x]])
}
})

# merge area2, area3, and area4
n234_upregulated_genes <- lapply(seq_along(CD11c_ctrl_pval_filtered), function(x){merge(merge(CD11c_experimental_pval_filtered1[[x]], CD11c_experimental_pval_filtered2[[x]], by = "X"), CD11c_experimental_pval_filtered3[[x]], by = "X")
})
n234_up_diff_genes <- lapply(n234_upregulated_genes, function(x) length(x$X)) 
n234_up_diff_genes
lapply(seq_along(n234_upregulated_genes), function(x) {if (n234_up_diff_genes[[x]] > 0) {
  write.csv(n234_upregulated_genes[[x]], 
            file = paste0(outs_directory,"Shared_Genes9_", file_label_deseq2[[x]],"__",file_label_deseq3[[x]],"__",file_label_deseq4[[x]],".csv"), 
            row.names = TRUE)
} else {
  print(n234_up_diff_genes[[x]])
}
})

# merge area1, area3, and area4
n134_upregulated_genes <- lapply(seq_along(CD11c_ctrl_pval_filtered), function(x){merge(merge(CD11c_ctrl_pval_filtered[[x]], CD11c_experimental_pval_filtered2[[x]], by = "X"), CD11c_experimental_pval_filtered3[[x]], by = "X")
})
n134_up_diff_genes <- lapply(n134_upregulated_genes, function(x) length(x$X)) 
n134_up_diff_genes
lapply(seq_along(n134_upregulated_genes), function(x) {if (n134_up_diff_genes[[x]] > 0) {
  write.csv(n134_upregulated_genes[[x]], 
            file = paste0(outs_directory,"Shared_Genes10_",file_label_deseq1[[x]],"__",file_label_deseq3[[x]],"__",file_label_deseq4[[x]],".csv"), 
            row.names = TRUE)
} else {
  print(n134_up_diff_genes[[x]])
}
})
# part II ###########################

### diagram not labeling ## potential solution below
venn_diagrams <- lapply(seq_along(CD11c_ctrl_genes_upreg), function(x){
  grid.newpage()
  draw.quad.venn(
    area1 = CD11c_ctrl_genes_upreg[[x]], 
    area2 = CD11c_experimental_genes_upreg1[[x]],
    area3 = CD11c_experimental_genes_upreg2[[x]],
    area4 = CD11c_experimental_genes_upreg3[[x]],
    n12 = n12_up_diff_genes[[x]],
    n23 = n23_up_diff_genes[[x]],
    n13 = n13_up_diff_genes[[x]],
    n123 = n123_up_diff_genes[[x]],
    n14 = n14_up_diff_genes[[x]],
    n24 = n24_up_diff_genes[[x]],
    n34 = n34_up_diff_genes[[x]],
    n124 = n124_up_diff_genes[[x]],
    n234 = n234_up_diff_genes[[x]],
    n134 = n134_up_diff_genes[[x]],
    n1234 = shared_up_diff_genes[[x]],
    category = c(control_cat_names_edit[[x]],exp1_cat_names_edit[[x]],exp2_cat_names_edit[[x]],exp3_cat_names_edit[[x]]),
    fill = c("Red","Yellow","Green","Blue"),
    label.col = rep("black",15), # specify label color
    cex = rep(2, 15), # size of area's label
    cat.cex = rep(0.5, 4), # specify category label size
    cat.fontfamily = rep("sans", 4), # specify category label font family
    lwd = rep(2, 4), 
    lty = rep("solid", 4), 
    alpha = rep(0.5, 4), 
    col = rep("black", 4),
    ind = TRUE,
    fontface = rep("plain", 15), 
    fontfamily = rep("serif", 15),
    cat.pos = c(-10, 5, -5, -5), 
    cat.dist = c(0.21, 0.21, 0.11, 0.11),
    print.mode = "raw",
    cat.just = rep(list(c(0.5, 0.5)), 4),
    rotation.degree = 0,
    rotation.centre = c(0.5, 0.5),
    cex.prop = NULL
  )
})
venn_diagrams <- lapply(venn_diagrams, function(x){grid.arrange(gTree(children =x), # Add title & subtitle
                                                                top = "Venn Diagram of Significant Differentially Expressed Genes",
                                                                bottom = paste0("CD11chi versus CD11clo\n","LFC cutoffs: ",LFC," and ",-LFC,"\n","Adjusted P-value cutoff: 0.05"))})


lapply(seq_along(venn_diagrams), function(x) {ggsave(filename = paste0(outs_directory,"Significant_Genes_VennDiagram_",file_label_deseq1[[x]],"__",file_label_deseq2[[x]],".pdf"),plot = venn_diagrams[[x]], width = 10, height = 8)})

# write the csv only if it has shared genes 
lapply(seq_along(venn_diagrams), function(x) {if (shared_up_diff_genes[[x]] > 0) {
  write.csv(shared_upregulated_genes[[x]], 
            file = paste0(outs_directory, "Shared_Genes_All_", file_label_deseq1[[x]], "__", file_label_deseq2[[x]], "__", file_label_deseq3[[x]],"__", file_label_deseq4[[x]],".csv"), 
            row.names = TRUE)
} else {
  print(shared_up_diff_genes[[x]])
}
})

