.fixOutliers <- function(gt, opts) {
    ## remove gross outliers
    if (any(gt$target)) {
        gt$lr[ gt$target] <- .smooth.outliers.gr(gt[ gt$target], "lr")
    }
    if (any(!gt$target)) {
        gt$lr[!gt$target] <- .smooth.outliers.gr(gt[!gt$target], "lr")
    }
    return(gt)
}

.smoothLogRatio <- function(tile, opts) {
    lr.smooth <- unname(unlist(lapply(split(tile$lr.gc, tile$arm), function(y) {
        if (length(y) >= opts$lr.smooth.window) {
            y.smooth <- suppressWarnings(
                hybrid.filter(y, opts$lr.smooth.window, method=c("MEAN", "MED"))$level[["MEAN"]]
            )
        } else {
            y.smooth <- y
        }
        return(y.smooth)
    })))
    return(lr.smooth)
}

smoothLogRatio <- function(cnv, opts) {
    cnv$tile$lr.smooth <- cnv$tile$lr.gc
    cnv$tile$lr.smooth[tile$target] <- .smoothLogRatio(cnv$tile[cnv$tile$target], opts)
    return(cnv)
}
