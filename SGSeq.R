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

# Getting arguments from terminal
args <- commandArgs(trailingOnly = TRUE)
db_path <- args[1]            # Path to DB (grasses.db)
sra_id <- args[2]             # sample ID (SRA ID)
bam_file_path <- args[3]      # Path to BAM - Output from STAR 
gff_path <- args[4]           # Path to GFF - from Phytozome
output_bed_path <- args[5]    # Path to output - BED (path, to directory)
num_cores <- as.integer(args[6]) # Max cores for processing

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


# Save as bed
exportSGVariantsToBED <- function(sg_variant_counts, file_name, sra_id) {
  read_counts <- rowSums(assay(sg_variant_counts))
  variant_data <- as.data.frame(rowData(sg_variant_counts))
  
  gene_names <- ifelse(!is.null(variant_data$geneName), variant_data$geneName, "NA")
  transcript_names <- ifelse(!is.null(variant_data$txName), variant_data$txName, "NA")
  
    bed_data <- data.frame(
    chrom = tolower(sub("^[^:]+:([^:]+):.*$", "\\1", variant_data$from)),  
    start = as.integer(sub("^[^:]+:[^:]+:(\\d+):.*$", "\\1", variant_data$from)),
    end = as.integer(sub("^[^:]+:[^:]+:(\\d+):.*$", "\\1", variant_data$to)),
    name = paste0(
      "ID=", variant_data$variantName,
      ";read_count=", read_counts,
      ";SRA_ID=", sra_id,
      ";gene_name=", gene_names,
      ";transcript_name=", transcript_names
    ),    
    score = 0,
    strand = sub(".*:([^:]+)$", "\\1", variant_data$from)
  )
  
  bed_data <- within(bed_data, {
    temp_start <- pmin(start, end)
    temp_end <- pmax(start, end)
    start <- temp_start
    end <- temp_end
    
    thickStart <- start 
    thickEnd <- start
    
    
    rm(temp_start, temp_end) 
  })
  
  track_line <- 'track name="Splicing Variants - SGSeq" description="Splicing variants identified with SGSeq" visibility=3 color=85,107,47 type=bed'
  write(track_line, file = file_name) 
  
  write.table(bed_data, file = file_name, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE, append = TRUE)
}

output <- paste0(output_bed_path, "/SGSeq_", sra_id, ".bed")

exportSGVariantsToBED(annotated_variants, output, sra_id)



# Running example
# Rscript SGSeq.R /home/bia/LandscapeSplicingGrasses/data/grasses.db SRR25663954 /home/bia/LandscapeSplicingGrasses/Results/star/Athaliana_447Aligned.sortedByCoord.out.bam /home/bia/LandscapeSplicingGrasses/data/Phytozome/PhytozomeV12/early_release/Athaliana_447_Araport11/annotation/Athaliana_447_Araport11.gene.gff3.gz /home/bia/LandscapeSplicingGrasses/Results/test_SGSeq.bed 1



