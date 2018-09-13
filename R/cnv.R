.getCopy <- function(gt, sex) {
    copy <- rep(2L, length(gt))
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
    return(copy)
}


.getGene <- function(gtf.fn) {
    gtf <- import(gtf.fn)
    gtf.gene <- gtf[gtf$type=="gene",c("source", "type", "gene_id", "gene_name")]
    return(gtf.gene)
}
