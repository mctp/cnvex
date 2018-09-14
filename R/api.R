addCorrectLogRatio <- function(cnv, opts) {
    cnv$tile$lr.gc <- .correctGC(cnv$tile, opts)
    return(cnv)
}

addSmoothLogRatio <- function(cnv, opts) {
    cnv$tile$lr.smooth <- cnv$tile$lr.gc
    if (opts$lr.smooth=="hybrid") {
        ## hybrid smooth only on targeted
        cnv$tile$lr.smooth[cnv$tile$target] <- .smoothLogRatio(cnv$tile[cnv$tile$target], opts)
    } else if (opts$lr.smooth=="outlier") {
        cnv$tile$lr.smooth <- .smoothOutliers(cnv$tile, opts)
    }
    return(cnv)
}

addJointSegment <- function(cnv, opts) {
    ## estimated standard-deviation
    sd.lr <- estimateSd(cnv$tile$lr.smooth) * opts$seg.sd.lr.penalty
    sd.baf <- estimateSd(cnv$tile$baf) * opts$seg.sd.baf.penalty
    ## create segmentation
    cnv$tile <- .addJointSeg(cnv$tile, sd.lr, sd.baf, opts$seg.method,
                             opts$tile.width, opts$seg.len.min, opts$seg.cbs.lr,
                             opts$seg.cbs.baf, opts$seg.rbs.selection, opts$seg.sd.prune, opts$seg.len.prune)
    cnv$seg <- unname(unlist(range(split(cnv$tile, cnv$tile$seg))))
    ## seg -> var
    cnv$var$seg <- findOverlaps(cnv$var, cnv$seg, select="first", maxgap = opts$tile.shoulder-1)
    ## seg -> gene
    hit.prom <- as.data.table(findOverlaps(promoters(cnv$gene, 0, 1), cnv$seg))
    hit.gene <- as.data.table(findOverlaps(cnv$gene, cnv$seg))
    hit.gene <- hit.gene[!(queryHits %in% hit.prom$queryHits)]
    hit.gene <- hit.gene[,.SD[1],by=queryHits]
    hit <- rbind(hit.prom, hit.gene)
    cnv$gene$seg <- NA_real_
    cnv$gene$seg[hit$queryHits] <- hit$subjectHits
    return(cnv)
}

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
