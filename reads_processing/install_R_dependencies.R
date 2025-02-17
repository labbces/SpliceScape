# Load or install BiocManager
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}

# Install Bioconductor and CRAN packages
BiocManager::install(c(
    "SGSeq",
    "txdbmaker",
    "GenomicFeatures",
    "DBI",
    "RSQLite"
))
