addBaf <- function(cnv, opts) {
    cnv$tile <- .getBaf(cnv$tile, cnv$var, opts)
    return(cnv)
}

addCopy <- function(cnv, opts) {
    if (is.null(opts$sex)) {
        sex <- .detect.sex(cnv$var, cnv$tile)
    } else {
        sex <- opts$sex
    }
    copy <- rep(2L, length(cnv$tile))
    if (sex=="male") {
        x.copy <- 1L
        y.copy <- 1L
    } else if (sex=="female"){
        x.copy <- 2L
        y.copy <- 0L
    } else {
        stop("wrong sex.")
    }
    copy[as.logical(seqnames(gt) %in% c("chrX"))] <- x.copy
    copy[as.logical(seqnames(gt) %in% c("chrY"))] <- y.copy
    cnv$tile$nC <- copy
    return(cnv)
}

addGene <- function(cnv, opts) {
    gtf <- import(opts$gene.fn)
    gtf.gene <- gtf[gtf$type=="gene",c("source", "type", "gene_id", "gene_name")]
    cnv$gene <- gtf.gene
    return(cnv)
}
