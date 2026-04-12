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
output_dir <- "/home/bia.estevam/landscapeSplicingGrasses/testing_splicescape_improvements/SGSeq_parser/output"

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
# "D:Chr5:26882759:-"  →  c("D", "Chr5", "26882759", "-")
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
# junction and coordinates
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
# head(flat_mcols[, c("variantType","variantName","countsVariant5p","countsEvent5p","countsVariant3p","countsEvent3p","PSI")])

# -----------------------
# Initial Parser
# -----------------------
# -----------------------
# A5SS: junction coords
# -----------------------
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

flat_mcols_A5SS <- flat_mcols[flat_mcols$variantType %in% c("A5SS:D", "A5SS:P"), , drop = FALSE]
flat_mcols_A5SS <- add_A5SS_junctions(flat_mcols_A5SS, Jdf)

write.table(
  flat_mcols_A5SS,
  file = file.path(output_dir, "sgseq_A5SS_with_junctions.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

# Quick inspect
cat("A5SS rows:", nrow(flat_mcols_A5SS), "\n")
cat("A5SS with coords:", sum(!is.na(flat_mcols_A5SS$junc_start) & !is.na(flat_mcols_A5SS$junc_end)), "\n")

# -----------------------
# A3SS: junction coords
# -----------------------
# choose distal/prox junction within each event
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

flat_mcols_A3SS <- flat_mcols[flat_mcols$variantType %in% c("A3SS:D", "A3SS:P"), , drop = FALSE]
flat_mcols_A3SS <- add_A3SS_junctions(flat_mcols_A3SS, Jdf)

write.table(
  flat_mcols_A3SS,
  file = file.path(output_dir, "sgseq_A3SS_with_junctions.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

# Quick inspect
# cat("A3SS rows:", nrow(flat_mcols_A3SS), "\n")
# cat("A3SS with coords:", sum(!is.na(flat_mcols_A3SS$junc_start) & !is.na(flat_mcols_A3SS$junc_end)), "\n")

# -----------------------
# RI: junction coords + featureID match
# -----------------------
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
  file = file.path(output_dir, "sgseq_RI_with_junctions.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

# cat("RI rows:", nrow(flat_mcols_RI), "\n")
# cat("RI matched junction featureID:", sum(!is.na(flat_mcols_RI$ri_junc_featureID)), "\n")

# -----------------------
# SE + S2E junction extraction
# -----------------------

flat_mcols_SE_raw <- flat_mcols[
  flat_mcols$variantType %in% c("SE:S", "SE:I", "S2E:S", "S2E:I"),
  , drop = FALSE
]

# Skipping junctions (SE:S / S2E:S)
# featureID5p stores the skipping junction featureID directly
se_skip <- flat_mcols_SE_raw[
  flat_mcols_SE_raw$variantType %in% c("SE:S", "S2E:S"), ]


skip_fID <- as.integer(se_skip$featureID5p)
idx_skip <- match(skip_fID, Jdf$featureID)

Jdf_SE_skipping <- data.frame(
  seqnames   = Jdf$seqnames[idx_skip],
  start      = Jdf$start[idx_skip],
  end        = Jdf$end[idx_skip],
  strand     = Jdf$strand[idx_skip],
  featureID  = skip_fID,
  key        = Jdf$key[idx_skip],
  PSI        = se_skip$PSI,
  countsVariant5p = se_skip$countsVariant5p,
  countsVariant3p = se_skip$countsVariant3p,
  event_type = "SE:S",
  geneID     = se_skip$geneID,
  eventID    = se_skip$eventID,
  geneName   = se_skip$geneName,
  stringsAsFactors = FALSE
)

# Inclusion junctions (SE:I / S2E:I)
# Use featureID5p + featureID3p of the :I row, minus the skipping
# junction featureID taken from the matching :S row

se_incl <- flat_mcols_SE_raw[
  flat_mcols_SE_raw$variantType %in% c("SE:I", "S2E:I"), ]

# Lookup table: skipping featureID per event, from the SE:S rows
se_skip_lookup <- se_skip[, c("geneID", "eventID", "featureID5p")]

rows_list <- vector("list", nrow(se_incl))

for (i in seq_len(nrow(se_incl))) {
  # Get the skipping junction featureID from the matching SE:S row
  match_skip <- se_skip_lookup[
    se_skip_lookup$geneID  == se_incl$geneID[i] &
    se_skip_lookup$eventID == se_incl$eventID[i],
  ]
  skip_fID_i <- as.integer(trimws(match_skip$featureID5p))

  # All featureIDs from 5p and 3p of this SE:I row
  fIDs_5p  <- as.integer(trimws(strsplit(se_incl$featureID5p[i], ",")[[1]]))
  fIDs_3p  <- as.integer(trimws(strsplit(se_incl$featureID3p[i], ",")[[1]]))
  all_fIDs <- unique(c(fIDs_5p, fIDs_3p))

  # Remove only the skipping junction featureID
  incl_fIDs <- setdiff(all_fIDs, skip_fID_i)

  if (length(incl_fIDs) == 0) next

  idx   <- match(incl_fIDs, Jdf$featureID)
  valid <- !is.na(idx)
  if (!any(valid)) next

  rows_list[[i]] <- data.frame(
    seqnames   = Jdf$seqnames[idx[valid]],
    start      = Jdf$start[idx[valid]],
    end        = Jdf$end[idx[valid]],
    strand     = Jdf$strand[idx[valid]],
    featureID  = incl_fIDs[valid],
    key        = Jdf$key[idx[valid]],
    PSI        = se_incl$PSI[i],
    countsVariant5p = se_incl$countsVariant5p[i],
    countsVariant3p = se_incl$countsVariant3p[i],
    event_type = "SE:I",
    geneID     = se_incl$geneID[i],
    eventID    = se_incl$eventID[i],
    geneName   = se_incl$geneName[i],
    stringsAsFactors = FALSE
  )
}

Jdf_SE_inclusion <- do.call(rbind, rows_list[!sapply(rows_list, is.null)])

# Add event_id using skipping junction coordinates
event_id_lookup <- data.frame(
  geneID   = Jdf_SE_skipping$geneID,
  eventID  = Jdf_SE_skipping$eventID,
  event_id = paste(
    Jdf_SE_skipping$seqnames,
    Jdf_SE_skipping$strand,
    Jdf_SE_skipping$start,
    Jdf_SE_skipping$end,
    sep = "_"
  ),
  stringsAsFactors = FALSE
)

# Combine, add event_id, and save
Jdf_SE <- rbind(Jdf_SE_skipping, Jdf_SE_inclusion)
Jdf_SE <- Jdf_SE[order(Jdf_SE$geneID, Jdf_SE$eventID, Jdf_SE$event_type), ]

Jdf_SE$event_id <- event_id_lookup$event_id[
  match(
    paste(Jdf_SE$geneID, Jdf_SE$eventID, sep = "_"),
    paste(event_id_lookup$geneID, event_id_lookup$eventID, sep = "_")
  )
]

write.table(
  Jdf_SE,
  file = file.path(output_dir, "sgseq_SE_junctions.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

# -----------------------
# MXE
# -----------------------
flat_mcols_MXE <- flat_mcols[flat_mcols$variantType == "MXE", , drop = FALSE]

rows_list <- vector("list", nrow(flat_mcols_MXE))

for (i in seq_len(nrow(flat_mcols_MXE))) {
  fID_5p <- as.integer(trimws(flat_mcols_MXE$featureID5p[i]))
  fID_3p <- as.integer(trimws(flat_mcols_MXE$featureID3p[i]))

  junc_5p <- Jdf[Jdf$featureID == fID_5p, , drop = FALSE]
  junc_3p <- Jdf[Jdf$featureID == fID_3p, , drop = FALSE]

  if (nrow(junc_5p) == 0 || nrow(junc_3p) == 0) {
    warning("MXE row ", i, ": featureID not found in Jdf")
    next
  }

  # One row per junction (long format)
  # Both junctions share the same PSI, counts, geneID, eventID
  rows_list[[i]] <- rbind(
    data.frame(
      seqnames        = junc_5p$seqnames,
      start           = junc_5p$start,
      end             = junc_5p$end,
      strand          = junc_5p$strand,
      featureID       = junc_5p$featureID,
      key             = junc_5p$key,
      PSI             = flat_mcols_MXE$PSI[i],
      countsVariant5p = flat_mcols_MXE$countsVariant5p[i],
      countsVariant3p = flat_mcols_MXE$countsVariant3p[i],
      event_type      = "MXE",
      geneID          = flat_mcols_MXE$geneID[i],
      eventID         = flat_mcols_MXE$eventID[i],
      geneName        = flat_mcols_MXE$geneName[i],
      stringsAsFactors = FALSE
    ),
    data.frame(
      seqnames        = junc_3p$seqnames,
      start           = junc_3p$start,
      end             = junc_3p$end,
      strand          = junc_3p$strand,
      featureID       = junc_3p$featureID,
      key             = junc_3p$key,
      PSI             = flat_mcols_MXE$PSI[i],
      countsVariant5p = flat_mcols_MXE$countsVariant5p[i],
      countsVariant3p = flat_mcols_MXE$countsVariant3p[i],
      event_type      = "MXE",
      geneID          = flat_mcols_MXE$geneID[i],
      eventID         = flat_mcols_MXE$eventID[i],
      geneName        = flat_mcols_MXE$geneName[i],
      stringsAsFactors = FALSE
    )
  )
}

Jdf_MXE <- do.call(rbind, rows_list[!sapply(rows_list, is.null)])

# event_id: shared anchor coordinates (min start and max end of the event)
# Built per geneID+eventID grouping
Jdf_MXE$event_id <- NA_character_
event_keys <- unique(paste(Jdf_MXE$geneID, Jdf_MXE$eventID, sep = "_"))

for (k in event_keys) {
  rows <- which(paste(Jdf_MXE$geneID, Jdf_MXE$eventID, sep = "_") == k)
  Jdf_MXE$event_id[rows] <- paste(
    Jdf_MXE$seqnames[rows[1]],
    Jdf_MXE$strand[rows[1]],
    min(Jdf_MXE$start[rows]),
    max(Jdf_MXE$end[rows]),
    sep = "_"
  )
}

Jdf_MXE <- Jdf_MXE[order(Jdf_MXE$geneID, Jdf_MXE$eventID, Jdf_MXE$start), ]

write.table(
  Jdf_MXE,
  file = file.path(output_dir, "sgseq_MXE_junctions.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

# cat("=== MXE junction table summary ===\n")
# cat("Total junction rows:", nrow(Jdf_MXE), "\n")
# cat("Total events:       ", length(unique(paste(Jdf_MXE$geneID, Jdf_MXE$eventID))), "\n")
# cat("Output written to:  ", file.path(output_dir, "sgseq_MXE_junctions.tsv"), "\n")

# -----------------------
# AFE
# -----------------------
flat_mcols_AFE <- flat_mcols[flat_mcols$variantType == "AFE", , drop = FALSE]

rows_list <- vector("list", nrow(flat_mcols_AFE))

for (i in seq_len(nrow(flat_mcols_AFE))) {
  # Split comma-separated featureID3p values
  fIDs <- as.integer(trimws(strsplit(flat_mcols_AFE$featureID3p[i], ",")[[1]]))
  fIDs <- fIDs[!is.na(fIDs)]

  if (length(fIDs) == 0) next

  idx   <- match(fIDs, Jdf$featureID)
  valid <- !is.na(idx)

  if (!any(valid)) next

  rows_list[[i]] <- data.frame(
    seqnames   = Jdf$seqnames[idx[valid]],
    start      = Jdf$start[idx[valid]],
    end        = Jdf$end[idx[valid]],
    strand     = Jdf$strand[idx[valid]],
    featureID  = fIDs[valid],
    key        = Jdf$key[idx[valid]],
    PSI        = flat_mcols_AFE$PSI[i],
    countsVariant5p = flat_mcols_AFE$countsVariant5p[i],
    countsVariant3p = flat_mcols_AFE$countsVariant3p[i],
    event_type = "AFE",
    geneID     = flat_mcols_AFE$geneID[i],
    eventID    = flat_mcols_AFE$eventID[i],
    geneName   = flat_mcols_AFE$geneName[i],
    stringsAsFactors = FALSE
  )
}

Jdf_AFE <- do.call(rbind, rows_list[!sapply(rows_list, is.null)])

Jdf_AFE <- Jdf_AFE[order(Jdf_AFE$geneID, Jdf_AFE$eventID, Jdf_AFE$start), ]

write.table(
  Jdf_AFE,
  file = file.path(output_dir, "sgseq_AFE_junctions.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

# cat("=== AFE junction table summary ===\n")
# cat("Total junction rows:", nrow(Jdf_AFE), "\n")
# cat("Total events:       ", length(unique(paste(Jdf_AFE$geneID, Jdf_AFE$eventID))), "\n")
# cat("Output written to:  ", file.path(output_dir, "sgseq_AFE_junctions.tsv"), "\n")

# -----------------------
# ALE
# -----------------------
flat_mcols_ALE <- flat_mcols[flat_mcols$variantType == "ALE", , drop = FALSE]

rows_list <- vector("list", nrow(flat_mcols_ALE))

for (i in seq_len(nrow(flat_mcols_ALE))) {
  # Split comma-separated featureID5p values
  fIDs <- as.integer(trimws(strsplit(flat_mcols_ALE$featureID5p[i], ",")[[1]]))
  fIDs <- fIDs[!is.na(fIDs)]

  if (length(fIDs) == 0) next

  idx   <- match(fIDs, Jdf$featureID)
  valid <- !is.na(idx)

  if (!any(valid)) next

  rows_list[[i]] <- data.frame(
    seqnames   = Jdf$seqnames[idx[valid]],
    start      = Jdf$start[idx[valid]],
    end        = Jdf$end[idx[valid]],
    strand     = Jdf$strand[idx[valid]],
    featureID  = fIDs[valid],
    key        = Jdf$key[idx[valid]],
    PSI        = flat_mcols_ALE$PSI[i],
    countsVariant5p = flat_mcols_ALE$countsVariant5p[i],
    countsVariant3p = flat_mcols_ALE$countsVariant3p[i],
    event_type = "ALE",
    geneID     = flat_mcols_ALE$geneID[i],
    eventID    = flat_mcols_ALE$eventID[i],
    geneName   = flat_mcols_ALE$geneName[i],
    stringsAsFactors = FALSE
  )
}

Jdf_ALE <- do.call(rbind, rows_list[!sapply(rows_list, is.null)])
Jdf_ALE <- Jdf_ALE[order(Jdf_ALE$geneID, Jdf_ALE$eventID, Jdf_ALE$start), ]

write.table(
  Jdf_ALE,
  file = file.path(output_dir, "sgseq_ALE_junctions.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

# -----------------------
# Constitutive junctions
# -----------------------
# Helper: parse all comma-separated featureIDs from a column
parse_all_fIDs <- function(x) {
  unique(na.omit(as.integer(trimws(unlist(strsplit(x, ","))))))
}

# Collect all featureIDs from every featureID column in flat_mcols
# This captures ALL alternative events including combined types like "SE:S, A3SS:D"
alt_fIDs <- unique(c(
  parse_all_fIDs(flat_mcols$featureID5p),
  parse_all_fIDs(flat_mcols$featureID3p),
  parse_all_fIDs(flat_mcols$featureID5pEvent),
  parse_all_fIDs(flat_mcols$featureID3pEvent)
))

# Junctions in Jdf NOT in any alternative event
Jdf_constitutive <- Jdf[!Jdf$featureID %in% alt_fIDs, ]

# Add PSI = 1 and event_type
Jdf_constitutive$PSI        <- 1
Jdf_constitutive$event_type <- "constitutive"

# Pull raw junction counts from analysis_features
# For constitutive junctions there are no variant counts (no alternative event),
# so we use the raw counts assay from analysis_features instead
feat_counts <- assays(analysis_features)$counts
feat_fIDs   <- SGSeq::featureID(analysis_features)
idx_const   <- match(Jdf_constitutive$featureID, feat_fIDs)

Jdf_constitutive$counts         <- as.numeric(feat_counts[idx_const, sra_id])

write.table(
  Jdf_constitutive,
  file = file.path(output_dir, "sgseq_constitutive_junctions.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

# -----------------------
# Done
# -----------------------
cat("Outputs written to:", output_dir, "\n")