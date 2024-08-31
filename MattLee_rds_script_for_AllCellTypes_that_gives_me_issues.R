.libPaths("/home/ttm3567/63_tylert/Analysis_Algorithms/Seurat_CellBrowser_Renv_ONLY/")
rds <- readRDS("/home/ttm3567/63_tylert/Sarcoidosis/redo/Matts_RDS/Sarcoidosis_All_CellTypes_MattLee.rds")
View(rds)
rds@assays$RNA@counts <- rds@assays[["RNA"]]@data
library(Seurat)

rds@meta.data[["Gender"]] <- as.factor(rds@meta.data[["Gender"]])
rds@active.ident <- rds@meta.data[["Gender"]]
Idents(rds) <- rds@meta.data[["Gender"]]
rds@misc$markers <- FindAllMarkers(rds)
# install.packages('devtools')
# devtools::install_github('immunogenomics/presto')
