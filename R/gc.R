.correctGC <- function(gt, opts) {
    lr <- gt$lr.raw
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
    return(lr)
}

gcLogRatio <- function(cnv, opts) {
    cnv$tile$lr.gc <- .correctGC(cnv$tile, opts)
    return(cnv)
}
