#!/usr/bin/env Rscript
# RDS-backed, read-only SGSeq A3SS/A5SS export audit.
# Usage: Rscript audit_sgseq_a3a5_rds.R MANIFEST NEW_OUTPUT_DIR
# Requires SGSeq and GenomicRanges (same environment as SGSeq exporter).
args <- commandArgs(trailingOnly=TRUE)
if (length(args)!=2L) stop("Usage: Rscript audit_sgseq_a3a5_rds.R manifest.tsv NEW_OUTPUT_DIR")
manifest_path <- args[[1]]; output <- args[[2]]
if (file.exists(output)) stop("Refusing existing output directory: ",output)
suppressPackageStartupMessages({library(SGSeq);library(GenomicRanges)})
manifest <- read.delim(manifest_path,header=FALSE,sep="\t",stringsAsFactors=FALSE,quote="")
if (ncol(manifest)!=3L || anyDuplicated(manifest[[1]])) stop("Expected unique sample IDs and 3 manifest columns")
parse_node <- function(x) {
  v <- strsplit(as.character(x),":",fixed=TRUE)[[1]]
  if (length(v)!=4L || !(v[[1]] %in% c("D","A")) || !(v[[4]] %in% c("+","-")))
    return(NULL)
  pos <- suppressWarnings(as.integer(v[[3]]))
  if (is.na(pos)) return(NULL)
  list(tag=v[[1]],chr=v[[2]],pos=pos,strand=v[[4]])
}
records <- list(); n <- 0L; errors <- character()
for (i in seq_len(nrow(manifest))) {
  sample <- manifest[i,1]; dir <- manifest[i,3]
  fp <- file.path(dir,"sgseq_features.rds")
  vp <- file.path(dir,"sgseq_variants.rds")
  if (!file.exists(fp) || !file.exists(vp)) {
    errors <- c(errors,paste(sample,"missing RDS",fp,vp)); next
  }
  f <- readRDS(fp); v <- readRDS(vp)
  gr <- granges(f)
  j <- gr[SGSeq::type(gr)=="J"]
  catalog <- data.frame(chr=as.character(seqnames(j)),start=start(j),end=end(j),
                        strand=as.character(strand(j)),featureID=SGSeq::featureID(j),
                        stringsAsFactors=FALSE)
  variants <- as.data.frame(mcols(v),stringsAsFactors=FALSE)
  required <- c("variantType","from","to","eventID","geneID")
  if (!all(required %in% names(variants))) {
    errors <- c(errors,paste(sample,"missing variant columns",paste(setdiff(required,names(variants)),collapse=",")));next
  }
  for (typ in c("A3SS","A5SS")) {
    src <- variants[as.character(variants$variantType) %in% paste0(typ,c(":P",":D")),,drop=FALSE]
    p <- file.path(dir,paste0("sgseq_",typ,".tsv"))
    if (!file.exists(p)) {errors<-c(errors,paste(sample,"missing TSV",p));next}
    exported <- read.delim(p,sep="\t",header=TRUE,check.names=FALSE,stringsAsFactors=FALSE)
    if (!all(c("geneID","eventID","variantType","from","to","junc_start","junc_end") %in% names(exported))) {
      errors<-c(errors,paste(sample,typ,"missing TSV columns"));next
    }
    # Each variant's own graph endpoints specify the corresponding junction.
    # Compare both against the sample-wide extrema selected by the exporter.
    for (k in seq_len(nrow(exported))) {
      r <- exported[k,,drop=FALSE]
      a<-parse_node(r$from);b<-parse_node(r$to)
      if (is.null(a)||is.null(b)||a$chr!=b$chr||a$strand!=b$strand) next
      expected <- sort(c(a$pos,b$pos))
      chosen <- suppressWarnings(c(as.integer(r$junc_start),as.integer(r$junc_end)))
      cand <- if (typ=="A3SS") {
        if (a$strand=="+") catalog[catalog$chr==a$chr & catalog$strand==a$strand & catalog$start==a$pos,,drop=FALSE]
        else catalog[catalog$chr==a$chr & catalog$strand==a$strand & catalog$end==a$pos,,drop=FALSE]
      } else {
        if (a$strand=="+") catalog[catalog$chr==a$chr & catalog$strand==a$strand & catalog$end==b$pos,,drop=FALSE]
        else catalog[catalog$chr==a$chr & catalog$strand==a$strand & catalog$start==b$pos,,drop=FALSE]
      }
      # Does the original SGSeq RDS contain this exact variant graph definition?
      idx <- which(as.character(src$geneID)==as.character(r$geneID) &
                   as.character(src$eventID)==as.character(r$eventID) &
                   as.character(src$variantType)==as.character(r$variantType) &
                   as.character(src$from)==as.character(r$from) &
                   as.character(src$to)==as.character(r$to))
      n<-n+1L
      records[[n]] <- data.frame(sample=sample,type=typ,geneID=as.character(r$geneID),
        eventID=as.character(r$eventID),variantType=as.character(r$variantType),
        from=as.character(r$from),to=as.character(r$to),
        expected_start=expected[[1]],expected_end=expected[[2]],
        exported_start=chosen[[1]],exported_end=chosen[[2]],
        matches_graph=all(!is.na(chosen))&&identical(as.integer(chosen),as.integer(expected)),
        n_matching_rds_variants=length(idx),n_catalog_candidates=nrow(cand),
        exact_junction_in_catalog=any(cand$start==expected[[1]] & cand$end==expected[[2]]),
        exported_junction_in_catalog=all(!is.na(chosen))&&any(cand$start==chosen[[1]] & cand$end==chosen[[2]]),
        stringsAsFactors=FALSE)
    }
  }
}
dir.create(output,recursive=TRUE,showWarnings=FALSE)
result <- if(length(records)) do.call(rbind,records) else data.frame()
write.table(result,file=file.path(output,"a3_a5_rds_comparison.tsv"),
            sep="\t",row.names=FALSE,quote=FALSE)
writeLines(errors,file.path(output,"missing_or_invalid_inputs.txt"))
if(nrow(result)) {
 print(aggregate(cbind(matches_graph,n_matching_rds_variants,exact_junction_in_catalog,
                       exported_junction_in_catalog)~type,data=result,FUN=sum))
 print(table(result$type,result$n_catalog_candidates>2))
}
cat("Records:",nrow(result)," Input warnings:",length(errors),"\n")
