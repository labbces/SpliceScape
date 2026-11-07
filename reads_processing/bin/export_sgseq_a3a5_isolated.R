#!/usr/bin/env Rscript
# Read-only experimental A3SS/A5SS exporter from native SGSeq RDS.
# Usage: Rscript export_sgseq_a3a5_isolated.R --self-test
#        Rscript export_sgseq_a3a5_isolated.R manifest.tsv NEW_OUTPUT_DIR
# No genome FASTA/GFF, BAM, production TSV rewrite, or SQLite required.

parse_node <- function(x) {
  a <- strsplit(as.character(x), ":", fixed=TRUE)[[1]]
  if (length(a)!=4L || !(a[1] %in% c("D","A")) ||
      !(a[4] %in% c("+","-"))) return(NULL)
  pos <- suppressWarnings(as.integer(a[3]))
  if(is.na(pos)) return(NULL)
  list(tag=a[1],chr=a[2],pos=pos,strand=a[4])
}
parse_feature_ids <- function(x) {
  if(is.list(x)) x <- unlist(x,recursive=TRUE,use.names=FALSE)
  ids <- unlist(strsplit(paste(as.character(x),collapse=","),",",fixed=TRUE),
                use.names=FALSE)
  ids <- trimws(ids)
  ids <- ids[nzchar(ids) & !is.na(ids) & ids!="NA"]
  if(!length(ids) || any(!grepl("^[0-9]+$",ids))) return(integer())
  unique(as.integer(ids))
}
# Pure mapper: source variants and source junction catalog are passed explicitly.
# Only one-to-one, graph-consistent junction variants are exported as 'mapped'.
map_variants <- function(v,catalog) {
  rv <- c("geneID","eventID","variantID","variantType","type","featureID","from","to")
  rc <- c("chr","start","end","strand","featureID")
  if(!all(rv %in% names(v))) stop("Missing variant columns: ",
                                    paste(setdiff(rv,names(v)),collapse=","))
  if(!all(rc %in% names(catalog))) stop("Missing catalog columns")
  ind <- which(grepl("^A(3|5)SS:(P|D)$",as.character(v$variantType)))
  rows <- vector("list",length(ind))
  for(k in seq_along(ind)) {
    i<-ind[k];typ<-as.character(v$variantType[i])
    kind<-sub(":.*","",typ);role<-sub(".*:","",typ)
    a<-parse_node(v$from[i]);b<-parse_node(v$to[i])
    graph_ok<-!is.null(a) && !is.null(b) &&
       a$tag=="D" && b$tag=="A" &&
       a$chr==b$chr && a$strand==b$strand
    ft<-as.character(v$type[i])
    ids<-parse_feature_ids(v$featureID[[i]])
    status<-"ok";hit<-catalog[FALSE,,drop=FALSE]
    if(!graph_ok) status<-"invalid_graph"
    else if(ft!="J") status<-"non_junction_variant"
    else if(length(ids)!=1L) status<-"non_unique_feature_reference"
    else {
      hit<-catalog[!is.na(catalog$featureID)&catalog$featureID==ids[1],,drop=FALSE]
      if(nrow(hit)==0L) status<-"feature_not_found"
      else if(nrow(hit)>1L) status<-"duplicate_catalog_feature_id"
      else if(hit$chr[1]!=a$chr || hit$strand[1]!=a$strand)
        status<-"chromosome_or_strand_mismatch"
      else {
        observed<-as.integer(c(hit$start[1],hit$end[1]))
        expected<-as.integer(sort(c(a$pos,b$pos)))
        if(role=="D" && !identical(observed,expected))
          status<-"distal_endpoint_mismatch"
        if(role=="P") {
          fixed<-if(kind=="A3SS") a$pos else b$pos
          site<-if(kind=="A3SS") {
            if(a$strand=="+") observed[1] else observed[2]
          } else if(a$strand=="+") observed[2] else observed[1]
          if(site!=fixed) status<-"proximal_fixed_site_mismatch"
        }
      }
    }
    rows[[k]]<-data.frame(
      geneID=as.character(v$geneID[i]),eventID=as.character(v$eventID[i]),
      variantID=as.character(v$variantID[i]),variantType=typ,
      from=as.character(v$from[i]),to=as.character(v$to[i]),
      variantFeatureType=ft,variantFeatureIDs=paste(ids,collapse=","),
      junctionFeatureID=if(nrow(hit)==1L) as.integer(hit$featureID[1]) else NA_integer_,
      seqnames=if(graph_ok) a$chr else NA_character_,
      strand=if(graph_ok) a$strand else NA_character_,
      junc_start=if(nrow(hit)==1L) as.integer(hit$start[1]) else NA_integer_,
      junc_end=if(nrow(hit)==1L) as.integer(hit$end[1]) else NA_integer_,
      mapping_status=status,mapped=status=="ok",rds_row=i,
      stringsAsFactors=FALSE)
  }
  if(!length(rows)) return(data.frame())
  do.call(rbind,rows)
}
run_tests <- function() {
  j<-data.frame(chr=c("Chr1","Chr1","Chr1","Chr2"),
    start=c(100L,100L,100L,300L),end=c(200L,220L,260L,500L),
    strand=c("+","+","+","-"),featureID=c(1L,2L,3L,4L))
  v<-data.frame(geneID=c("G","G","G","G2"),eventID=c("E","E","E","E2"),
    variantID=c("1","2","3","4"),variantType=c("A3SS:P","A3SS:D","A3SS:D","A5SS:D"),
    type=rep("J",4),featureID=c("1","2","3","4"),
    from=c(rep("D:Chr1:100:+",3),"D:Chr2:500:-"),
    to=c("A:Chr1:260:+","A:Chr1:220:+","A:Chr1:260:+","A:Chr2:300:-"))
  out<-map_variants(v,j)
  stopifnot(all(out$mapped),identical(out$junc_end,c(200L,220L,260L,500L)))
  # Extra junction in whole catalog never changes variant-specific assignment.
  j2<-rbind(j,data.frame(chr="Chr1",start=100L,end=999L,strand="+",featureID=5L))
  stopifnot(identical(map_variants(v,j)$junc_end,map_variants(v,j2)$junc_end))
  v$featureID[2]<-"5"
  stopifnot(map_variants(v,j2)$mapping_status[2]=="distal_endpoint_mismatch")
  v$type[2]<-"JEJ"
  stopifnot(map_variants(v,j2)$mapping_status[2]=="non_junction_variant")
  message("PASS: RDS feature identity, proximal/distal mapping, strand, extra catalog features")
}
get_assay <- function(x,nm,sample,n) {
  m<-x[[nm]]
  if(is.null(m)) return(rep(NA_real_,n))
  if(is.null(dim(m)) || nrow(m)!=n) stop("Assay dimensions invalid: ",nm)
  col<-if(!is.null(colnames(m)) && sample %in% colnames(m)) {
    match(sample,colnames(m))
  } else if(ncol(m)==1L) 1L else stop("Ambiguous sample in assay: ",nm)
  as.numeric(m[,col])
}
save_tsv <- function(x,path) {
  write.table(x,path,sep="\t",quote=FALSE,row.names=FALSE,na="NA")
}
run_export <- function(manifest_path,outdir) {
  if(file.exists(outdir)) stop("Refusing to overwrite existing output: ",outdir)
  suppressPackageStartupMessages({
    library(SGSeq);library(GenomicRanges);library(S4Vectors);library(SummarizedExperiment)
  })
  mf<-read.delim(manifest_path,header=FALSE,sep="\t",quote="",stringsAsFactors=FALSE)
  if(ncol(mf)!=3L || nrow(mf)==0L || anyDuplicated(mf[[1L]]))
    stop("Expected headerless sample/MAJIQ_dir/SGSeq_dir manifest with unique samples")
  all<-list();warnings<-character()
  for(i in seq_len(nrow(mf))) {
    sample<-mf[i,1];dir<-mf[i,3]
    fp<-file.path(dir,"sgseq_features.rds");vp<-file.path(dir,"sgseq_variants.rds")
    if(!file.exists(fp)||!file.exists(vp)) {
      warnings<-c(warnings,paste(sample,"missing RDS",fp,vp));next
    }
    f<-readRDS(fp);v<-readRDS(vp)
    j<-GenomicRanges::granges(f)
    j<-j[SGSeq::type(j)=="J"]
    catalog<-data.frame(chr=as.character(GenomicRanges::seqnames(j)),
      start=as.integer(GenomicRanges::start(j)),end=as.integer(GenomicRanges::end(j)),
      strand=as.character(GenomicRanges::strand(j)),
      featureID=as.integer(SGSeq::featureID(j)))
    meta<-as.data.frame(S4Vectors::mcols(v),stringsAsFactors=FALSE)
    rec<-map_variants(meta,catalog)
    if(!nrow(rec)) next
    q<-SummarizedExperiment::assays(v)
    for(nm in c("variantFreq","countsVariant5p","countsVariant3p",
                "countsEvent5p","countsEvent3p")) {
      val<-get_assay(q,nm,sample,nrow(meta))
      rec[[if(nm=="variantFreq") "PSI" else nm]]<-val[rec$rds_row]
    }
    rec$sample<-sample
    all[[length(all)+1L]]<-rec
  }
  res<-if(length(all)) do.call(rbind,all) else data.frame()
  dir.create(outdir,recursive=TRUE)
  save_tsv(res,file.path(outdir,"a3a5_variants_junctions.tsv"))
  save_tsv(if(nrow(res)) res[!res$mapped,,drop=FALSE] else data.frame(),
           file.path(outdir,"mapping_issues.tsv"))
  writeLines(warnings,file.path(outdir,"missing_inputs.txt"))
  if(nrow(res)) {
    print(table(res$variantType,res$mapping_status))
    cat("Records:",nrow(res)," mapped:",sum(res$mapped),"\n")
  }
}
if(sys.nframe()==0L) {
  a<-commandArgs(trailingOnly=TRUE)
  if(length(a)==1L && a[1L]=="--self-test") run_tests()
  else if(length(a)==2L) run_export(a[1L],a[2L])
  else stop("Usage: Rscript export_sgseq_a3a5_isolated.R --self-test\n",
            "or: Rscript export_sgseq_a3a5_isolated.R MANIFEST NEW_OUTPUT_DIR")
}
