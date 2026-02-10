#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(GenomicFeatures)
  library(SGSeq)
  library(txdbmaker)
  library(GenomicRanges)
})

# -----------------------
# User inputs
# -----------------------
sra_id <- "SRR12642286"
bam_file_path <- "/Storage/data2/riano/testManualSpliceScape/SRR12642286/03_star/ath_SRR12642286.Aligned.sortedByCoord.out.bam"
gtf_path <- "/Storage/data2/riano/data/Athaliana_447_Araport11.gene_exons.gtf"
cores <- 30
output_dir <- "/Storage/data2/riano/testManualSpliceScape/SRR12642286/05_sgseq/"

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

#Validate user inputs
if (!file.exists(bam_file_path)) {
  stop("BAM file not found: ", bam_file_path)
}
if (!file.exists(gtf_path)) {
  stop("GTF file not found: ", gtf_path)
}

# -----------------------
# Helpers
# -----------------------

# Flatten list-like columns from mcols export
flatten_mcols <- function(x) {
  df <- as.data.frame(x)
  df[] <- lapply(df, function(v) if (is.list(v)) vapply(v, toString, character(1)) else v)
  df
}

# Parse "D:Chr5:26882759:-" or "A:Chr5:26882588:-" (vectorized)
parse_DA_vec <- function(x) {
  p <- strsplit(as.character(x), ":", fixed = TRUE)
  # Validar número de campos
  n_fields <- vapply(p, length, integer(1))
  if (any(n_fields != 4)) {
    warning("Some entries don't have 4 fields. Found: ", 
            paste(unique(n_fields), collapse = ", "))
  }
  data.frame(
    tag    = vapply(p, `[[`, character(1), 1),   # "D" or "A"
    chr    = vapply(p, `[[`, character(1), 2),
    pos    = as.integer(vapply(p, `[[`, character(1), 3)),
    strand = vapply(p, `[[`, character(1), 4),
    stringsAsFactors = FALSE
  )
}

# Build junction catalog from analysis_features
make_Jdf <- function(analysis_features) {
  gr <- granges(analysis_features)
  J  <- gr[SGSeq::type(gr) == "J"]

  if (length(J) == 0) {
    stop("No junction features found in analysis_features")
  }

  Jdf <- data.frame(
    seqnames  = as.character(seqnames(J)),
    start     = start(J),
    end       = end(J),
    strand    = as.character(strand(J)),
    featureID = SGSeq::featureID(J),
    stringsAsFactors = FALSE
  )
  stopifnot("featureID" %in% colnames(Jdf), nrow(Jdf) > 0)

  Jdf$key <- paste(Jdf$seqnames, Jdf$strand, Jdf$start, Jdf$end, sep = "|")
  Jdf
}

# Generic: add junction coords + featureID by matching from/to against Jdf
add_junction_match <- function(df, Jdf, prefix = "junc") {
  pf <- parse_DA_vec(df$from) # donor
  pt <- parse_DA_vec(df$to)   # acceptor

  js <- pmin(pf$pos, pt$pos)
  je <- pmax(pf$pos, pt$pos)

  k <- paste(pf$chr, pf$strand, js, je, sep = "|")
  idx <- match(k, Jdf$key)

  df[[paste0(prefix, "_start")]] <- js
  df[[paste0(prefix, "_end")]]   <- je
  df[[paste0(prefix, "_featureID")]] <- ifelse(is.na(idx), NA_integer_, Jdf$featureID[idx])

  df
}

# A5SS-specific: choose distal/prox junction within each event
add_A5SS_junctions <- function(df_A5SS, Jdf) {
  if (nrow(df_A5SS) == 0) {
    warning("No A5SS events to process")
    return(df_A5SS)
  }
  df_A5SS$junc_start <- NA_integer_
  df_A5SS$junc_end   <- NA_integer_

  pf <- parse_DA_vec(df_A5SS$from)
  pt <- parse_DA_vec(df_A5SS$to)

  key_vec <- paste(df_A5SS$geneID, df_A5SS$eventID, sep = "_")
  keys <- unique(key_vec)

  for (k in keys) {
    rows <- which(key_vec == k)
    if (!length(rows)) next

    chr <- pf$chr[rows[1]]
    strand_ev <- pf$strand[rows[1]]
    acc <- pt$pos[rows[1]]

    if (strand_ev == "+") {
      # acceptor fixed at end; donor varies at start
      cand <- Jdf[Jdf$seqnames == chr & Jdf$strand == strand_ev & Jdf$end == acc, , drop = FALSE]
      if (nrow(cand) < 1) next
      j_dist <- cand[which.min(cand$start), , drop = FALSE]  # distal = smaller donor coord
      j_prox <- cand[which.max(cand$start), , drop = FALSE]  # proximal = larger donor coord
    } else {
      # strand '-' : acceptor fixed at start; donor varies at end
      cand <- Jdf[Jdf$seqnames == chr & Jdf$strand == strand_ev & Jdf$start == acc, , drop = FALSE]
      if (nrow(cand) < 1) next
      j_dist <- cand[which.max(cand$end), , drop = FALSE]    # distal = larger donor coord
      j_prox <- cand[which.min(cand$end), , drop = FALSE]    # proximal = smaller donor coord
    }

    for (r in rows) {
      sel <- if (df_A5SS$variantType[r] == "A5SS:D") j_dist else j_prox
      df_A5SS$junc_start[r] <- sel$start
      df_A5SS$junc_end[r]   <- sel$end
    }
  }

  df_A5SS
}

# A3SS-specific: choose distal/prox junction within each event
add_A3SS_junctions <- function(df_A3SS, Jdf) {
  if (nrow(df_A3SS) == 0) {
    warning("No A3SS events to process")
    return(df_A3SS)
  }
  df_A3SS$junc_start <- NA_integer_
  df_A3SS$junc_end   <- NA_integer_

  pf <- parse_DA_vec(df_A3SS$from)  # donor
  pt <- parse_DA_vec(df_A3SS$to)    # acceptor (not always needed, but parsed for completeness)

  key_vec <- paste(df_A3SS$geneID, df_A3SS$eventID, sep = "_")
  keys <- unique(key_vec)

  for (k in keys) {
    rows <- which(key_vec == k)
    if (!length(rows)) next

    chr <- pf$chr[rows[1]]
    strand_ev <- pf$strand[rows[1]]
    don <- pf$pos[rows[1]]  # donor fixed for A3SS

    if (strand_ev == "+") {
      # '+' : donor fixed at start; acceptor varies at end
      cand <- Jdf[Jdf$seqnames == chr & Jdf$strand == strand_ev & Jdf$start == don, , drop = FALSE]
      if (nrow(cand) < 1) next

      # distal acceptor = larger genomic coord (more downstream)
      j_dist <- cand[which.max(cand$end), , drop = FALSE]
      j_prox <- cand[which.min(cand$end), , drop = FALSE]

    } else {
      # '-' : donor fixed at end; acceptor varies at start
      cand <- Jdf[Jdf$seqnames == chr & Jdf$strand == strand_ev & Jdf$end == don, , drop = FALSE]
      if (nrow(cand) < 1) next

      # distal acceptor = smaller genomic coord (more upstream in transcript)
      j_dist <- cand[which.min(cand$start), , drop = FALSE]
      j_prox <- cand[which.max(cand$start), , drop = FALSE]
    }

    for (r in rows) {
      sel <- if (df_A3SS$variantType[r] == "A3SS:D") j_dist else j_prox
      df_A3SS$junc_start[r] <- sel$start
      df_A3SS$junc_end[r]   <- sel$end
    }
  }

  df_A3SS
}

# Safely pull a sample column from an assay matrix-like object
get_sample_vec <- function(assay_obj, sample_name) {
  if (is.null(assay_obj)) {
    warning("Assay object is NULL for sample: ", sample_name)
    return(NULL)
  }
  
  m <- assay_obj
  # m can be matrix, dgCMatrix, DelayedArray, etc.
  if (!is.null(dim(m))) {
    cn <- colnames(m)
    if (!is.null(cn) && sample_name %in% cn) {
      return(as.numeric(m[, sample_name]))
    }
    # single-sample case with no colnames
    if (ncol(m) == 1) {
      return(as.numeric(m[, 1]))
    }
  }
  NULL
}

# -----------------------
# Build TxDb / TxFeatures
# -----------------------
txdb <- txdbmaker::makeTxDbFromGFF(gtf_path, format = "gtf")
txFeatures <- convertToTxFeatures(txdb)

# -----------------------
# Sample info + BAM
# -----------------------
si <- data.frame(
  sample_name = sra_id,
  file_bam    = bam_file_path,
  stringsAsFactors = FALSE
)

bam <- getBamInfo(si, yieldSize = NULL, cores = cores)

# -----------------------
# SGSeq analysis
# -----------------------
analysis_features <- analyzeFeatures(
  bam,
  annotation = txFeatures,
  predict  = TRUE,
  cores    = cores
)

analysis_variants <- analyzeVariants(
  analysis_features,
  cores = cores
)

# Save features object
saveRDS(analysis_features, file = file.path(output_dir, "sgseq_features.rds"))
saveRDS(analysis_variants, file = file.path(output_dir, "sgseq_variants.rds"))


# -----------------------
# Build junction catalog + export mcols
# -----------------------
Jdf <- make_Jdf(analysis_features)

flat_mcols <- flatten_mcols(mcols(analysis_variants))

# pull assays
quant <- assays(analysis_variants)

flat_mcols$countsVariant5p <- get_sample_vec(quant$countsVariant5p, sra_id) # supporting reads for 5' choice (donor side)
flat_mcols$countsVariant3p <- get_sample_vec(quant$countsVariant3p, sra_id) # supporting reads for 3' choice (acceptor side)
flat_mcols$countsEvent5p   <- get_sample_vec(quant$countsEvent5p,   sra_id) # total reads for event (5' side)
flat_mcols$countsEvent3p   <- get_sample_vec(quant$countsEvent3p,   sra_id) # total reads for event (3' side)
flat_mcols$PSI             <- get_sample_vec(quant$variantFreq,     sra_id) # PSI / variant frequency

# sanity checks
stopifnot(nrow(flat_mcols) == length(flat_mcols$PSI))
head(flat_mcols[, c("variantType","variantName","countsVariant5p","countsEvent5p","countsVariant3p","countsEvent3p","PSI")])

# -----------------------
# A5SS: junction coords
# -----------------------
flat_mcols_A5SS <- flat_mcols[flat_mcols$variantType %in% c("A5SS:D", "A5SS:P"), , drop = FALSE]
flat_mcols_A5SS <- add_A5SS_junctions(flat_mcols_A5SS, Jdf)

write.table(
  flat_mcols_A5SS,
  file = file.path(output_dir, "sgseq_variants_A5SS_with_junctions.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

# Quick inspect
cat("A5SS rows:", nrow(flat_mcols_A5SS), "\n")
cat("A5SS with coords:", sum(!is.na(flat_mcols_A5SS$junc_start) & !is.na(flat_mcols_A5SS$junc_end)), "\n")

# -----------------------
# A3SS: junction coords
# -----------------------
flat_mcols_A3SS <- flat_mcols[flat_mcols$variantType %in% c("A3SS:D", "A3SS:P"), , drop = FALSE]
flat_mcols_A3SS <- add_A3SS_junctions(flat_mcols_A3SS, Jdf)

write.table(
  flat_mcols_A3SS,
  file = file.path(output_dir, "sgseq_variants_A3SS_with_junctions.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

# Quick inspect
cat("A3SS rows:", nrow(flat_mcols_A3SS), "\n")
cat("A3SS with coords:", sum(!is.na(flat_mcols_A3SS$junc_start) & !is.na(flat_mcols_A3SS$junc_end)), "\n")

# -----------------------
# RI: junction coords + featureID match
# -----------------------
flat_mcols_RI <- flat_mcols[flat_mcols$variantType %in% c("RI:E", "RI:R"), , drop = FALSE]
flat_mcols_RI <- add_junction_match(flat_mcols_RI, Jdf, prefix = "ri_junc")

# Optional: define "new exon" coords for RI:R as the intron interval (exonized intron)
is_RIR <- flat_mcols_RI$variantType == "RI:R"
flat_mcols_RI$new_exon_seqnames <- NA_character_
flat_mcols_RI$new_exon_start    <- NA_integer_
flat_mcols_RI$new_exon_end      <- NA_integer_
flat_mcols_RI$new_exon_strand   <- NA_character_
flat_mcols_RI$new_exon_type     <- NA_character_

if (any(is_RIR)) {
  pf_RI <- parse_DA_vec(flat_mcols_RI$from)
  flat_mcols_RI$new_exon_seqnames[is_RIR] <- pf_RI$chr[is_RIR]
  flat_mcols_RI$new_exon_strand[is_RIR]   <- pf_RI$strand[is_RIR]
  flat_mcols_RI$new_exon_start[is_RIR]    <- pmin(flat_mcols_RI$ri_junc_start[is_RIR], flat_mcols_RI$ri_junc_end[is_RIR])
  flat_mcols_RI$new_exon_end[is_RIR]      <- pmax(flat_mcols_RI$ri_junc_start[is_RIR], flat_mcols_RI$ri_junc_end[is_RIR])
  flat_mcols_RI$new_exon_type[is_RIR]     <- "E_new"
}

write.table(
  flat_mcols_RI,
  file = file.path(output_dir, "sgseq_variants_RI_with_junctions.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

cat("RI rows:", nrow(flat_mcols_RI), "\n")
cat("RI matched junction featureID:", sum(!is.na(flat_mcols_RI$ri_junc_featureID)), "\n")

# -----------------------
# Done
# -----------------------
cat("Outputs written to:", output_dir, "\n")