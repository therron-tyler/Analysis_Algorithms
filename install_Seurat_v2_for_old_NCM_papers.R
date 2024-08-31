#LIBRARY PATH
libs <- .libPaths("/home/ttm3567/63_tylert/Analysis_Algorithms/Seurat_v2_CellBrowser_Renv_20240205/")

#install.packages("devtools")

#library(remotes)
#library(devtools)
# 
# devtools::install_version(
#   "hdf5r",
#   dependencies = FALSE)
# 
# devtools::install_version(
#   "metap",
#   dependencies = FALSE)
# 
# devtools::install_version(
#   "SDMTools",
#   dependencies = FALSE)
# 
# devtools::install_version(
#   "Seurat",
#   version = "2.3.4",
#   dependencies = FALSE)

#source("https://z.umn.edu/archived-seurat")

# Set CRAN mirror
options(repos = c(CRAN = "https://cran.r-project.org"))

install.packages('remotes')
# Replace '2.3.0' with your desired version
remotes::install_version(package = 'Seurat', version = package_version('2.3.0'))
library(Seurat)
