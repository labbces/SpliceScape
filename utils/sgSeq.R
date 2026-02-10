library("GenomicFeatures")
library("SGSeq")
library("txdbmaker")

sra_id <- "SRR12642286"
bam_file_path <- "/Storage/data2/riano/testManualSpliceScape/SRR12642286/03_star/ath_SRR12642286.Aligned.sortedByCoord.out.bam"
gtf_path<-"/Storage/data2/riano/data/Athaliana_447_Araport11.gene_exons.gtf"
cores <- 30
output_dir <- "/Storage/data2/riano/testManualSpliceScape/SRR12642286/05_sgseq/"
txdb <- txdbmaker::makeTxDbFromGFF(gtf_path, format = "gtf")
txFeatures <- convertToTxFeatures(txdb)

# --- Sample info
si <- data.frame(
  sample_name = sra_id,
  file_bam    = bam_file_path,
  stringsAsFactors = FALSE
)
# Column type for an SGFeatures object takes values
# J (splice junction)
# E (disjoint exon bin)
# D (splice donor site)
# A (splice acceptor site)

bam <- getBamInfo(si, yieldSize = NULL, cores = cores)
# --- Analyzing features and variants

# 
analysis_features <- analyzeFeatures(
  bam,
  features = txFeatures,   # usa anotação
  annotation = txFeatures, # usa anotação
  predict  = TRUE,         # adiciona de novo
  cores    = cores
)
analysis_variants <- analyzeVariants(
    analysis_features,
    cores = cores
)


# --- junction catalog from SGFeatures (analysis_features)
gr <- granges(analysis_features)
J  <- gr[SGSeq::type(gr) == "J"]

Jdf <- data.frame(
  seqnames  = as.character(seqnames(J)),
  start     = start(J),   # genomic start (always < end)
  end       = end(J),     # genomic end
  strand    = as.character(strand(J)),
  featureID = SGSeq::featureID(J),
  stringsAsFactors = FALSE
)

#Check that we have data in Jdf
stopifnot("featureID" %in% colnames(Jdf))
stopifnot(nrow(Jdf) > 0)

# --- parse helpers for "D:Chr5:26882759:-" and "A:Chr5:26882588:-"
parse_from <- function(x) {
  p <- strsplit(x, ":", fixed = TRUE)[[1]]
  list(chr = p[2], pos = as.integer(p[3]), strand = p[4])
}
parse_to <- function(x) {
  p <- strsplit(x, ":", fixed = TRUE)[[1]]
  list(pos = as.integer(p[3]))
}

flat_mcols <- as.data.frame(mcols(analysis_variants))
flat_mcols[] <- lapply(flat_mcols, function(x) if (is.list(x)) sapply(x, toString) else x)


###Selecting only A5SS events and assigning junction coordinates

flat_mcols_A5SS <- flat_mcols[flat_mcols$variantType %in% c("A5SS:D","A5SS:P"), ]
# --- output columns
flat_mcols_A5SS$junc_start <- NA_integer_
flat_mcols_A5SS$junc_end   <- NA_integer_

# --- per-event processing
key_vec <- paste(flat_mcols_A5SS$geneID, flat_mcols_A5SS$eventID, sep = "_")
keys <- unique(key_vec)

for (k in keys) {
  rows <- which(key_vec == k)
  if (length(rows) == 0) next

  # event metadata
  f <- parse_from(flat_mcols_A5SS$from[rows[1]])
  t <- parse_to(flat_mcols_A5SS$to[rows[1]])

  chr        <- f$chr
  acc        <- t$pos
  strand_ev  <- f$strand

  # pick candidate junctions + choose distal/prox correctly for strand
  if (strand_ev == "+") {
    # acceptor is at end; donor varies at start
    cand <- Jdf[Jdf$seqnames == chr & Jdf$strand == strand_ev & Jdf$end == acc, , drop = FALSE]
    if (nrow(cand) < 1) next

    j_dist <- cand[which.min(cand$start), , drop = FALSE]  # distal = smaller donor coord
    j_prox <- cand[which.max(cand$start), , drop = FALSE]  # proximal = larger donor coord

  } else {
    # strand '-' : acceptor is at start; donor varies at end
    cand <- Jdf[Jdf$seqnames == chr & Jdf$strand == strand_ev & Jdf$start == acc, , drop = FALSE]
    if (nrow(cand) < 1) next

    # distal/prox in transcript sense: distal corresponds to larger genomic donor coord on '-'
    j_dist <- cand[which.max(cand$end), , drop = FALSE]
    j_prox <- cand[which.min(cand$end), , drop = FALSE]
  }

  # assign one junction per variant
  for (r in rows) {
    vt <- flat_mcols_A5SS$variantType[r]
    sel <- if (vt == "A5SS:D") j_dist else j_prox
    flat_mcols_A5SS$junc_start[r] <- sel$start
    flat_mcols_A5SS$junc_end[r]   <- sel$end
  }
}

# --- result (inspect)
flat_mcols_A5SS[, c("geneName","eventID","variantType","from","to","junc_start","junc_end")]


###Selecting only RI events and assigning junction coordinates
#TODO: Check 
# RI:E (spliced) are events with removed introns, splicing happened.
#  the coordinates should be the last/first nucleotide of the exon
# RI:R (retained intron) are the retained introns and
#  the coordinates should be first/last nucleotide of the intron 
flat_mcols_RI <- flat_mcols[flat_mcols$variantType %in% c("RI:E","RI:R"), ]

# parse from/to com strsplit (sem regex)
parse_from_vec <- function(x) {
  p <- strsplit(as.character(x), ":", fixed=TRUE)
  chr <- vapply(p, `[[`, character(1), 2)
  pos <- as.integer(vapply(p, `[[`, character(1), 3))
  str <- vapply(p, `[[`, character(1), 4)
  data.frame(chr=chr, D=pos, strand=str, stringsAsFactors=FALSE)
}
parse_to_vec <- function(x) {
  p <- strsplit(as.character(x), ":", fixed=TRUE)
  pos <- as.integer(vapply(p, `[[`, character(1), 3))
  data.frame(A=pos, stringsAsFactors=FALSE)
}

pf <- parse_from_vec(flat_mcols_RI$from)
pt <- parse_to_vec(flat_mcols_RI$to)

# 4) calcular junction coords (start<end sempre)
flat_mcols_RI$ri_junc_start <- pmin(pf$D, pt$A)
flat_mcols_RI$ri_junc_end   <- pmax(pf$D, pt$A)

# 5) match por chave
Jdf$key <- paste(Jdf$seqnames, Jdf$strand, Jdf$start, Jdf$end, sep="|")
k <- paste(pf$chr, pf$strand, flat_mcols_RI$ri_junc_start, flat_mcols_RI$ri_junc_end, sep="|")

idx <- match(k, Jdf$key)
flat_mcols_RI$ri_junc_featureID <- ifelse(is.na(idx), NA_integer_, Jdf$featureID[idx])

# 6) diagnóstico rápido
cat("RI rows:", nrow(flat_mcols_RI), "\n")
cat("Matched featureID:", sum(!is.na(flat_mcols_RI$ri_junc_featureID)), "\n")

# veja alguns exemplos que bateram / não bateram
head(flat_mcols_RI[, c("from","to","ri_junc_start","ri_junc_end","ri_junc_featureID")], 10)
head(flat_mcols_RI[is.na(flat_mcols_RI$ri_junc_featureID),
                   c("from","to","ri_junc_start","ri_junc_end")], 10)


saveRDS(analysis_features, file = file.path(output_dir, "sgseq_features.rds"))