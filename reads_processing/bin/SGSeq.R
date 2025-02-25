#!/usr/bin/env Rscript
# Loading packages 
suppressMessages(library("optparse"))
suppressMessages(library("GenomicFeatures"))
suppressMessages(library("SGSeq"))
suppressMessages(library("txdbmaker"))


# Getting arguments from terminal
option_list = list(
  make_option(c("--gff"), type='character', default = NULL, help="Annotation GFF file from Phytozome", metavar='character'),
  make_option(c("--cores"), type='integer', default = 1, help="Number of cores"),
  make_option(c("--path_to_bam"), type='character', default = NULL, help="Path to BAM file(s)", metavar='character'),
  make_option(c("--sra_id"), type='character', default = NULL, help="Sample ID - as it is in database", metavar='character'),
  make_option(c("--out"), type='character', default = "", help="Output folder", metavar='character')
)
opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

#sra_id <- "SRR25663954"             # sample ID (SRA ID)
#bam_file_path <- "/home/bia/LandscapeSplicingGrasses/Results/star/Athaliana_447Aligned.sortedByCoord.out.bam"      # Path to BAM - Output from STAR 
#gff_path <-      "/home/bia/LandscapeSplicingGrasses/data/Phytozome/PhytozomeV12/early_release/Athaliana_447_Araport11/annotation/Athaliana_447_Araport11.gene_exons.gff3"      # Path to GFF - from Phytozome
#output_bed_path <- "/home/bia/LandscapeSplicingGrasses/Results/teste"    # Path to output - BED (path, to directory)
#num_cores <- 1

# Checking inputs
if (is.null(opt$gff) || is.null(opt$path_to_bam) || is.null(opt$sra_id))  {
  print_help(opt_parser)
  stop("Please provide valid -gff, --path_to_bam, and --sra_id arguments.", call.=FALSE)
}

cores = opt$cores
dir.create(opt$out, showWarnings = FALSE, recursive = TRUE)

# Load annotation file
txdb <- txdbmaker::makeTxDbFromGFF(opt$gff, format = "gff3")
txFeatures <- convertToTxFeatures(txdb)


# Get information from BAM - Output from STAR
sample_info <- data.frame(
  sample_name = c(opt$sra_id),
  file_bam = opt$path_to_bam
)

bam <- getBamInfo(sample_info, yieldSize = NULL, cores = cores)

# Analyze features e variants
analysis_features <- analyzeFeatures(bam, features = txFeatures, cores = cores)
analysis_variants <- analyzeVariants(analysis_features, cores = cores)

# Export results to CSV
sgvc_fpkm <- variantFreq(analysis_variants)

# Flatten the metadata columns
flat_mcols <- as.data.frame(mcols(analysis_variants))
flat_mcols[] <- lapply(flat_mcols, function(x) if (is.list(x)) sapply(x, toString) else x)

# Flatten the FPKM data
flat_fpkm <- as.data.frame(sgvc_fpkm)
flat_fpkm[] <- lapply(flat_fpkm, function(x) if (is.list(x)) sapply(x, toString) else x)

# Combine data
result <- cbind(flat_mcols, flat_fpkm)

# Write to CSV
output <- paste0(opt$out, "/SGSeq_", opt$sra_id, ".csv")
write.csv(result, file.path(output), row.names = FALSE)

# Running example
# Rscript SGSeq.R --gff /home/bia/LandscapeSplicingGrasses/data/Phytozome/PhytozomeV12/early_release/Athaliana_447_Araport11/annotation/Athaliana_447_Araport11.gene_exons.gff3 --cores 1 --path_to_bam /home/bia/LandscapeSplicingGrasses/Results/star/Athaliana_447Aligned.sortedByCoord.out.bam --sra_id SRR25663954 --out teste