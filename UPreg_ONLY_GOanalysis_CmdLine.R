#USAGE: Rscript GOanalysis_CmdLine.R [1]Path to CSV containing locations of DEseq files [2]Path to out_directory


libs <- .libPaths("/home/ttm3567/63_tylert/DEseq_libpath_dir/")

#BiocManager::install("org.Mm.eg.db")
#BiocManager::install("topGO")
#install.packages("rlang")
library(readr)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(DESeq2)
library(topGO)
library(rlang)
#library(biomaRt) - holding off on Gene Names for now

## Read command line arguments ------------------------
args <- commandArgs(TRUE)
DEseq_file_list_path <- as.character(args[1])# include '/' at the beginning and end of the file path
out_directory <- as.character(args[2]) # include '/' at the beginning and end of the file path

# UPregulated genes section =========================================================================================================

# the first step is to read in a list of csv files that are the DEseq output comparisons

comparisons <- read.csv(DEseq_file_list_path)

# directly write files to the remote sever.
outs_directory <- out_directory

# the list of file paths and comparison names should be read in as tables now 
comparison_table_list <- list()
for(i in 1:nrow(comparisons)){
  row <- read.csv(comparisons$comp_file[i])
  comparison_table_list[[i]] <- row
}
# have a list of comparison dataframes and each object in the list needs a name now
comp_name_list <- list(comparisons$comp_name)
for(i in 1:length(comp_name_list)){
  names(comparison_table_list) <- comp_name_list[[i]]
}
##### where the log2 foldchange column gets renamed for GO analysis #######

comparison_table_list <- lapply(comparison_table_list, function(x) {
  comparisons2 <- sub("^log2FoldChange.*","log2FoldChange",colnames(x))
  set_names(x, nm = comparisons2)
  })

##################################################################
comparison_table_list_filtered <- lapply(comparison_table_list, function(x) x %>% filter(log2FoldChange != 0))

comparison_table_list_pval <- lapply(comparison_table_list_filtered, function(x) x$padj)

comparison_table_list_pval_GOdata <- comparison_table_list_pval
for(i in 1:length(comparison_table_list_pval_GOdata)){
  as.numeric(comparison_table_list_pval_GOdata[[i]])
  #colnames(comparison_table_list_pval[[i]]) <- c("allScore", "allLogFoldChange")
  names(comparison_table_list_pval_GOdata[[i]]) <- comparison_table_list_filtered[[i]]$X
}#ensembleID = 'X'
# 3/8/23 - changing inputs from gene symbol to ensemble ID
# 1. comparison table list filtered is a list of df's that can be used to serve as the background gene set

bkgrd_geneNames <- lapply(comparison_table_list_filtered, function(x) x$X)#ensembleID

# 2. make my upregulated genes of interest from log2foldchange > 1 and adjusted pvalue < 0.05
GoI_table_list_filtered <- lapply(comparison_table_list, function(x) x %>% filter(log2FoldChange > 1 & padj < 0.05))

GoI <- lapply(GoI_table_list_filtered, function(x) x$X) # switch to ensembleID 

# mapply attempt
gl_fxn <- function(x,y){
  factor(as.integer(x %in% y))}

genelist_mapply <- as.list(mapply(gl_fxn,x=bkgrd_geneNames,y=GoI))

# ==================
names_gl_fxn <- function(x,y){
  names(x) <- y
  str(x)
  return(x)}
genelist_names_mapply <- as.list(mapply(names_gl_fxn,x=genelist_mapply,y=bkgrd_geneNames))
# this is the solution to the naming issue encountered earlier # ===================
# read in mapping Db
allGO2genes <- annFUN.org(whichOnto="BP", feasibleGenes=NULL, mapping="org.Mm.eg.db", ID="ensembl")

# create the topGOdata object

GOdata_list <- list()
for(i in 1:length(genelist_names_mapply)){
  GOdata <- new("topGOdata",
                ontology="BP",
                allGenes=genelist_names_mapply[[i]],
                annot=annFUN.GO2genes,
                GO2genes=allGO2genes)
  GOdata_list[[i]] <- GOdata
}
# stat analysis #### ==============================================================================

resultFisher_list <- list()
for(i in 1:length(GOdata_list)){
  resultFisher <- runTest(GOdata_list[[i]], algorithm = "classic", statistic = "fisher")
  resultFisher_list[[i]] <- resultFisher
}
resultFisher_list
# take the raw stats results and output to the GenTable then make into a dataframe and append to a new list 
result_table_list <- list()
for(i in 1:length(resultFisher_list)){
  allRes <- as.data.frame(GenTable(GOdata_list[[i]], classicFisher = resultFisher_list[[i]], orderBy = "classicFisher", ranksOf = "classicFisher", topNodes =150))
  result_table_list[[i]] <- allRes
}
# ================================= Created topGOdata object, and result table ==============================
# ========= Next part is adding the extra columns because then it will more closely match GOrilla output =================

####### 3/7/2023 UPDATED background genes ##################
bkgrdgenes_fxn <- function(x,y){
  backgroundgenesII <- nrow(as.data.frame(x@feasible) %>% filter(x@feasible == TRUE))# gets the actual number of background genes 
  y$backgroundgenesII <-backgroundgenesII
  return(y)
}
# put the inputs through the function
bkgrdgenes_list_mapply <- mapply(bkgrdgenes_fxn,x=GOdata_list,y=result_table_list)
bkgrd_genes_matrix <- as.data.frame(bkgrdgenes_list_mapply)
backgroundgenes <- lapply(bkgrd_genes_matrix, function(x) {(x$backgroundgenesII)})
names(backgroundgenes) <- names(comparison_table_list_filtered)

########## significant genes 3/2/2023 ########
siggenes_fxn <- function(x,y) {
  siggenes <- as.data.frame(slot(x, "geneData"))[2,]
  y$siggenes <- siggenes
  return(y)
}

siggenes_list_mapply <- mapply(siggenes_fxn,x=resultFisher_list,y=result_table_list)
sig_genes_matrix <- as.data.frame(siggenes_list_mapply)
colnames(sig_genes_matrix) <- names(comparison_table_list_filtered)
significantgenes <- lapply(sig_genes_matrix, function(x) {(x$siggenes)})
## end of sig genes ##

#### Enrichment Scores #####
lapply_enrich <- lapply(result_table_list, function(x) {
  Enrichment_Scores <- vector("numeric", nrow(x))
  for(row in 1:nrow(x)){
    Enrichment_Scores[row] <- x$Significant[row]/x$Expected[row]}
  return(Enrichment_Scores)})
names(lapply_enrich) <- names(comparison_table_list_filtered)
#####


# Next step is to combine these extra columns into one dataframe list to be written as GO analysis CSVs ===========
# cbind(), mapply() 

# lists to join: backgroundgenes, significantgenes, lapply_enrich, result_table_list 

# defining function to join elements of GO analysis together
df_out_fxn <- function(Background_Genes,Significant_Genes,Enrichment_Score,r) {
  CSV_out_df <- cbind(r,Enrichment_Score,Significant_Genes,Background_Genes)
  return(CSV_out_df)}

CSV_out_mapply <- as.data.frame(mapply(df_out_fxn, Background_Genes=backgroundgenes, Significant_Genes=significantgenes, r=result_table_list, Enrichment_Score=lapply_enrich, SIMPLIFY = TRUE))
CSV_GOanalysis_list <- as.list(CSV_out_mapply)
# this is the point where the dataframes for each input file are constructed, but still need their names changed - and have to add the p-adjusted column and gene symbol associated with each GOterm column ==========================================================================================
# library(biomaRt)
# ensembl = useMart("ensembl",dataset="mmusculus_gene_ensembl")
# # this can error out due to time limit reached for curl command ^^^^^^^^^^^^^^^
# # match the GO term to the Gene Symbol
# GOterm_gene <- lapply(CSV_GOanalysis_list, function(x) {
#   goID_togene_df <- as.list(x$GO.ID)
#   lapply(goID_togene_df, function(x) {
#     go = c(x)
#     Gene_Symbol_Entrez2 <- getBM(c("external_gene_name"),
#                                  filters = "go",
#                                  values = c(x),
#                                  mart = ensembl)})})

# this lapply function will take the old names of the results and update them
GOrenamed <- lapply(CSV_GOanalysis_list, function(x) {
  names(x) <- c("GO_Term","Term_Description","B","Significant","Expected","Fisher_p-value","Enrichment_Score","n","N")
  return(x)})

# GOterms_genesymbols <- lapply(GOterm_gene, function(y) { lapply(y, function(x){
#   return(paste(x, collapse = ","))})})
# 
# GOterms_genesymbols_fxn <- function(x,y) {
#   names(x) <- y$GO_Term
#   go_genesymbol <- t(as.data.frame(x))
# }
# # use mapply to combine the input lists so that the GOterm to gene symbol column can be made
# GOterms_genesymbols_mapply <- mapply(GOterms_genesymbols_fxn, x=GOterms_genesymbols, y=GOrenamed)
# GOterms_genesymbols_list <- as.list(as.data.frame(GOterms_genesymbols_mapply))


# ADD in adjusted p-value column
adjusted_pval_list <- lapply(resultFisher_list, function(x) {
  p_value <- x@score
  p = as.data.frame(p_value)
  p_arranged <- p %>% arrange(p_value)
  sorted_p2 <- as.list(p_arranged)
  BH_Adjusted_pvalue <- p.adjust(sorted_p2$p_value, method = "BH", n = length(sorted_p2$p_value))
  BH_p2 <- as.data.frame(BH_Adjusted_pvalue)
  filtered_BH_pval2 <- head(BH_p2,150)})

# combining the genes associated with each GOterm as one list for each file being analyzed
dfconcat_goterm_genesymbol_fxn <- function(x,BH_Adjusted_pvalue) {
  cbind(as.data.frame(x),BH_Adjusted_pvalue)}

dfconcat_goterm_genesymbol_mapply <- mapply(dfconcat_goterm_genesymbol_fxn, x=GOrenamed, BH_Adjusted_pvalue=adjusted_pval_list)
dfconcat_goterm_genesymbol_list <- as.list(as.data.frame(dfconcat_goterm_genesymbol_mapply))
# at this point there should be 10 columns and have been renamed


# the extra columns are added to the GOanalysis outputs and can be written to a CSV file. 
lapply(1:length(dfconcat_goterm_genesymbol_list), function(i) write.csv(dfconcat_goterm_genesymbol_list[[i]], 
                                                                        file = paste0(outs_directory,names(comparison_table_list[i]),"_UPregulated_GOterms",".csv"),
                                                                        row.names = TRUE))

#REVIGO input
revigo_list <- lapply(dfconcat_goterm_genesymbol_list, function(x) {
  as.data.frame(x) %>% dplyr::select(-Term_Description) %>% dplyr::select(-B) %>% dplyr::select(-N) %>% dplyr::select(-n) %>% dplyr::select(-Significant) %>% dplyr::select(-Expected) %>% dplyr::select(-Fisher_p.value) %>% dplyr::select(-Enrichment_Score)
})


lapply(1:length(revigo_list), function(i) write_tsv(revigo_list[[i]], 
                                                    file = paste0(outs_directory,names(comparison_table_list[i]),"_REVIGO_UPregulated_GOterms",".tsv")))
# ==============================================================================================================================================
# DOWNREGULATED GENES

