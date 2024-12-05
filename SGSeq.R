# Install and load packages
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("SGSeq", "txdbmaker", "GenomicFeatures"))

if (!require("DBI")) install.packages("DBI")
if (!require("RSQLite")) install.packages("RSQLite")

library(GenomicFeatures)
library(SGSeq)
library(txdbmaker)
library(DBI)
library(RSQLite)

# Define the arguments directly in the script
sra_id <- "SRR25663954"                                           # Sample ID
bam_file_path <- "/home/bia/LandscapeSplicingGrasses/Results/star/Athaliana_447Aligned.sortedByCoord.out.bam" # BAM file
gff_path <- "/home/bia/LandscapeSplicingGrasses/data/Phytozome/PhytozomeV12/early_release/Athaliana_447_Araport11/annotation/Athaliana_447_Araport11.gene.gff3.gz" # GFF file
output_bed_path <- "/home/bia/LandscapeSplicingGrasses/Results/testing_sgseq"   # Output directory
num_cores <- 1                                                   # Number of cores for processing
output_gff_path <- "/home/bia/LandscapeSplicingGrasses/Results/testing_sgseq"   # Output directory

# Connect to DB to retrieve metadata
conn <- dbConnect(RSQLite::SQLite(), db_path)

get_sample_metadata <- function(sra_id) {
  query <- sprintf("SELECT * FROM sra_metadata WHERE sra_id = '%s'", sra_id)
  sample_data <- dbGetQuery(conn, query)
  return(sample_data)
}

# Load annotation file
txdb <- txdbmaker::makeTxDbFromGFF(gff_path, format = "gff3")
txFeatures <- convertToTxFeatures(txdb)

# Get information from BAM - Output from STAR
sample_info <- data.frame(
  sample_name = c(sra_id),
  file_bam = bam_file_path
)

BAM <- getBamInfo(sample_info, yieldSize = NULL, cores = num_cores)

# Analyze features
analysis_results <- analyzeFeatures(BAM, which = NULL, features = txFeatures,
                                    alpha = 2, psi = 0, beta = 0.2,
                                    gamma = 0.2, min_junction_count = NULL,
                                    min_anchor = 1, min_n_sample = 1,
                                    min_overhang = NA, annotation = txFeatures,
                                    max_complexity = NA, verbose = TRUE, cores = num_cores)

# Analyze variants
analysis_variants <- analyzeVariants(analysis_results, maxnvariant = 20, include = "default",
                                     min_denominator = NA, min_anchor = 1, cores = num_cores)


# Annotate
annotated_variants <- annotate(analysis_variants, txFeatures)



# Save as bed
exportSGVariantsToBED <- function(sg_variant_counts, file_name, sra_id) {
  variant_data <- as.data.frame(rowData(sg_variant_counts))
  bed_data <- data.frame(
    chrom = sub("^[^:]+:([^:]+):.*$", "\\1", variant_data$from),
    start = as.integer(sub("^[^:]+:[^:]+:(\\d+):.*$", "\\1", variant_data$from)) - 1,
    end = as.integer(sub("^[^:]+:[^:]+:(\\d+):.*$", "\\1", variant_data$to)),
    name = variant_data$variantName,
    score = 0,
    strand = sub(".*:([^:]+)$", "\\1", variant_data$from),
    sra_id = sra_id  
  )
  
  write.table(bed_data, file = file_name, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
}

output <- paste(output_bed_path,"/","SGSeq","_",sra_id,".bam")
exportSGVariantsToBED(annotated_variants, output, sra_id)



