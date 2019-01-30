#' cnvex
#'
#' @name cnvex
#' @docType package
#' @importFrom parallel mclapply detectCores
#' @importFrom rtracklayer import export
#' @importFrom matrixStats rowMedians rowMads rowMins rowMaxs colMedians
#' @importFrom Rsamtools BamFile
#' @importFrom BSgenome.Hsapiens.UCSC.hg38 BSgenome.Hsapiens.UCSC.hg38
#' @importFrom GenomeInfoDb keepStandardChromosomes Seqinfo dropSeqlevels seqlevels keepSeqlevels seqlevels<- seqinfo<- seqnames seqnames<-
#' @importFrom GenomicRanges tileGenome reduce granges findOverlaps width pintersect mcols mcols<- start end start<- end<- split gaps strand promoters GRangesList
#' @importFrom S4Vectors queryHits subjectHits DataFrame %in% endoapply elementNROWS
#' @importFrom IRanges %over% 
#' @importFrom data.table data.table setkey as.data.table fread fwrite setDT rbindlist dcast.data.table copy melt :=
#' @importFrom stringr str_sub str_match str_replace str_split
#' @importFrom Biostrings getSeq letterFrequency
#' @importFrom VariantAnnotation readVcf ScanVcfParam info geno header info<- header<- info<- geno<- vcfWhich<- vcfWhich qual qual<- fixed
#' @importFrom limma loessFit
#' @importFrom DNAcopy CNA smooth.CNA
#' @importFrom DelayedArray rowRanges
#' @importFrom jointseg jointSeg estimateSd
#' @importFrom raster raster focal Which
#' @importFrom robfilter hybrid.filter
#' @importFrom fastICA fastICA
#' @importFrom MASS ginv
NULL
