

bioconductor_prerequisites <- function() {

  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }

  if(!requireNamespace("Biobase", quietly = TRUE)) {
    BiocManager::install()
    BiocManager::install("Biobase", version = "3.8")
  }
}
