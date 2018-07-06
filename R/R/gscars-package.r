#' gscars
#'
#' @name gscars
#' @docType package
#' @useDynLib gscars
#' @importFrom Rcpp sourceCpp
#' @importFrom parallel mclapply detectCores
#' @importFrom rtracklayer import export
#' @importFrom Rsamtools BamFile
#' @importFrom BSgenome.Hsapiens.UCSC.hg38 BSgenome.Hsapiens.UCSC.hg38
#' @importFrom GenomeInfoDb keepStandardChromosomes Seqinfo dropSeqlevels seqlevels keepSeqlevels seqlevels<- seqinfo<- seqnames seqnames<-
#' @importFrom GenomicRanges tileGenome reduce granges findOverlaps width pintersect mcols mcols<- start end start<- end<- split gaps strand
#' @importFrom S4Vectors queryHits subjectHits DataFrame %in% endoapply
#' @importFrom IRanges %over% 
#' @importFrom data.table data.table setkey as.data.table fread fwrite setDT rbindlist dcast.data.table copy melt
#' @importFrom stringr str_sub str_match str_replace
#' @importFrom Biostrings getSeq letterFrequency
#' @importFrom VariantAnnotation readVcf ScanVcfParam info geno header info<- header<- info<- geno<- vcfWhich<- vcfWhich qual qual<-
#' @importFrom limma loessFit
#' @importFrom DNAcopy CNA smooth.CNA
#' @importFrom DelayedArray rowRanges
#' @importFrom jointseg jointSeg estimateSd
#' @importFrom raster raster focal Which
NULL
