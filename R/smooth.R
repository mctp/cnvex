.smooth.outliers.gr <- function(gr, data.col) {
    obj <- CNA(mcols(gr)[[data.col]], as.integer(seqnames(gr)), floor((start(gr)+end(gr))/2),
              data.type = "logratio",
              sampleid = "sample")
    adj <- smooth.CNA(obj)$sample
    return(adj)
}

.smoothOutliers <- function(tile, opts) {
    ## remove gross outliers
    lr.smooth <- tile$lr.gc
    if (any(tile$target)) {
        lr.smooth[ tile$target] <- .smooth.outliers.gr(tile[ tile$target], "lr.gc")
    }
    if (any(!tile$target)) {
        lr.smooth[!tile$target] <- .smooth.outliers.gr(tile[!tile$target], "lr.gc")
    }
    return(lr.smooth)
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
    if (opts$lr.smooth=="hybrid") {
        ## hybrid smooth only on targeted
        cnv$tile$lr.smooth[cnv$tile$target] <- .smoothLogRatio(cnv$tile[cnv$tile$target], opts)
    } else if (opts$lr.smooth=="outlier") {
        cnv$tile$lr.smooth <- .smoothOutliers(cnv$tile, opts)
    }
    return(cnv)
}
