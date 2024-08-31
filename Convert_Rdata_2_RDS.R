#LIBRARY PATH
libs <- .libPaths("/home/ttm3567/63_tylert/DEseq_libpath_dir/")
#install.packages("Seurat")
library(Seurat)
library(dplyr)
library(tidyr)
library(readr)

load("/home/ttm3567/63_tylert/Sarcoidosis/sarcoidosis.harmonyObj.RData")


saveRDS(harmonyObj, file = "/home/ttm3567/63_tylert/Sarcoidosis/sarcoidosis.harmonyObj.rds")

load("/home/ttm3567/63_tylert/Sarcoidosis/sarcoidosis.macroHarmony.RData")


saveRDS(macroHarmony, file = "/home/ttm3567/63_tylert/Sarcoidosis/sarcoidosis.macroHarmony.rds")
