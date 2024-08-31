# Tyler 
#LIBRARY PATH
libs <- .libPaths("/home/ttm3567/63_tylert/Analysis_Algorithms/Seurat_v2_CellBrowser_Renv_20240205/")
# Set CRAN mirror
options(repos = c(CRAN = "https://cran.r-project.org"))

#install.packages("devtools")
library(devtools)

#devtools::install_version("Seurat", version = "2.3.4")


install.packages("/home/ttm3567/63_tylert/Analysis_Algorithms/adehabitat_1.8.20.tar", repos = NULL, type = "source")
install.packages("/home/ttm3567/63_tylert/Analysis_Algorithms/hdf5r_1.3.9.tar", repos = NULL, type = "source")
install.packages("/home/ttm3567/63_tylert/Analysis_Algorithms/SDMTools_1.1.tar", repos = NULL, type = "source")

#devtools::install_version("adehabitat", version = "1.8.20")
#devtools::install_version("SDMTools", version = "1.1")
#devtools::install_version("metap", version = "1.0")
#devtools::install_version("hdf5r", version = "1.3.9")

# # Function to install a package version and its dependencies
# install_package_version <- function(package, version) {
#   dependencies <- remotes::package_deps(package, repos = getOption("repos"), type = "source")
#   for (dep in dependencies$package) {
#     remotes::install_cran(dep, type = "source")
#   }
#   remotes::install_version(package, version = version, type = "source", upgrade = "never")
# }
# 
# # Specify required packages and their versions
# required_packages <- c('mixtools', 'fpc', 'VGAM', 'igraph', 'FNN', 'irlba', 'tclust', 'ranger', 'SDMTools', 'diffusionMap', 'Hmisc', 'metap', 'png')
# 
# # Attempt to install dependencies
# for (pkg in required_packages) {
#   install_package_version(pkg, version = "NA") # NA means install latest version compatible with R version
# }
# 
# # Finally, install Seurat version 2.3.0
# install_package_version('Seurat', '2.3.0')

# List of dependencies you might need to install manually
# required_packages <- c('mixtools', 'fpc', 'VGAM', 'igraph', 'FNN', 'irlba', 'tclust', 'ranger', 'SDMTools', 'diffusionMap', 'Hmisc', 'metap', 'png')
# 
# # Install dependencies
# for (pkg in required_packages) {
#   if (!requireNamespace(pkg, quietly = FALSE)) {
#     install.packages(pkg)
#   }
# }
# 
# # Install Seurat 2.3.0
# if (!requireNamespace("Seurat", quietly = FALSE)) {
#   remotes::install_version('Seurat', version = '2.3.0')
# }
#source("https://z.umn.edu/archived-seurat")

