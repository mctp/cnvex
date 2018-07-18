.normCoverage <- function(t.cov, n.cov, gt, opts) {
    mx.cov <- cbind(t.cov, n.cov)*width(gt)
    mx.cov.tar <- mx.cov[ gt$target,]
    mx.cov.tar <- t(t(mx.cov.tar)/(colSums(mx.cov.tar)/1e6))
    mx.cov.tar <- 1000*mx.cov.tar/width(gt[ gt$target])
    gc.cov.tar <- cbind(gc=gt[gt$target]$gc, mx.cov.tar)
    mx.cov.off <- mx.cov[!gt$target,]
    mx.cov.off <- t(t(mx.cov.off)/(colSums(mx.cov.off)/1e6))
    mx.cov.off <- 1000*mx.cov.off/width(gt[!gt$target])
    mx.cov[ gt$target,] <- mx.cov.tar
    mx.cov[!gt$target,] <- mx.cov.off
    return(mx.cov)
}

rawLogRatio <- function(t.bam, n.bam, tile, opts) {
    ## normalize by sequencing depth
    t.cov <- .runMosdepthTile(t.bam, tile, opts$cores)
    n.cov <- .runMosdepthTile(n.bam, tile, opts$cores)
    mx.cov <- .normCoverage(t.cov, n.cov, tile, opts)
    lr <- log2(mx.cov[,1] / mx.cov[,2])
    lr.raw <- ifelse(is.finite(lr), lr, NA_real_)
    mcols(tile) <- cbind(mcols(tile), cbind(mx.cov, lr.raw))
    return(tile)
}
