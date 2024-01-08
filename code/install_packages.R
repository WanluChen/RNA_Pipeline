#########################
## install_packages.R  ##
#########################

# Input Yes if you encounter such question:
# Do you want to attempt to install these from sources? (Yes/no/cancel)

cran.package.list <- c("BiocManager", "glue", "data.table", "magrittr", "rjson",
                       "readxl", "ggplot2", "ggfortify", "ggrepel", "dplyr", "circlize")

for(package in cran.package.list) {
  tryCatch({
    if(!require(package, character.only = TRUE)) install.packages(package)
    library(package, character.only = T)
  }, error=function(e){cat("ERROR: package", package, "download failed.\n")})
}

bioc.package.list <- c("DESeq2", "ComplexHeatmap", "scater", "GenomicRanges", "limma",
                       "tximeta", "SummarizedExperiment", "biomaRt", "DBI", "GSVA")

for(package in bioc.package.list) {
  tryCatch({
    if(!require(package, character.only = TRUE)) BiocManager::install(package)
    library(package, character.only = T)
  }, error=function(e){cat("ERROR: package", package, "download failed.\n")})
}
