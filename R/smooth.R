.smoothOutliers <- function(y, tile, opts) {
    ## remove gross outliers
    y.tar.n <- sum(!is.na(y[ tile$target]))
    if (y.tar.n>1) {
        y[ tile$target] <- .runSmooth(y[ tile$target], tile[ tile$target])
    }
    y.off.n <- sum(!is.na(y[!tile$target]))
    if (y.off.n>1) {
        y[!tile$target] <- .runSmooth(y[!tile$target], tile[!tile$target])
    }
    return(y)
}

.rawLogRatio <- function(t.cov, n.cov, opts) {
    lr.raw <- log2(t.cov / n.cov)
    lr.raw <- ifelse(is.finite(lr.raw), lr.raw, NA_real_)
    return(lr.raw)
}

.scaleLogRatio <- function(lr, gt, opts) {
    weight <- width(gt)
    target <- gt$target
    lr[ target] <- lr[ target] - weighted.mean(lr[ target], weight[ target], na.rm=TRUE)
    lr[!target] <- lr[!target] - weighted.mean(lr[!target], weight[!target], na.rm=TRUE)
    return(lr)
}

.gcLogRatio <- function(lr, gt, opts) {
    ## normalize, smooth, and gc-correct
    if (opts$gc.adjust.trend) {
        for (sel in unique(gt$target)) {
            lr.sel <- lr[gt$target==sel]
            gt.sel <- gt[gt$target==sel]
            weight.sel <- ifelse(gt.sel$blacklist==0 & gt.sel$unmasked>0.9, 1, 0)
            span.sel <- ifelse(sel, opts$gc.adjust.span.on, opts$gc.adjust.span.off)
            gc.residuals.sel <- limma::loessFit(y=lr.sel, x=gt.sel$gc,
                                                weight=weight.sel, min.weight=1e-9,
                                                span=span.sel)$residuals
            if (opts$gc.adjust.offset) {
                lr.offset.sel <- lm(gc.residuals.sel~lr.sel, weights=weight.sel)$coefficients[1]
                gc.residuals.sel <- gc.residuals.sel - lr.offset.sel
            }
            if (sel) {
                gc.range <- gt$gc > opts$gc.adjust.on[1] & gt$gc < opts$gc.adjust.on[2]
                lr[gt$target & gc.range] <- gc.residuals.sel[gc.range[gt$target]]
            } else {
                gc.range <- gt$gc > opts$gc.adjust.off[1] & gt$gc < opts$gc.adjust.off[2]
                lr[!gt$target & gc.range] <- gc.residuals.sel[gc.range[!gt$target]]
            }
        }
    }
    lr[!is.finite(lr)] <- NA_real_
    return(lr)
}

.smoothLogRatio <- function(lr, tile, opts) {
    lr.smooth <- unname(unlist(lapply(split(lr, tile$arm), function(y) {
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
