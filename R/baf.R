.getTargetBaf <- function(gt, snp, opts) {
    hits <- findOverlaps(snp, gt, maxgap = opts$shoulder-1)
    hits <- hits[gt[subjectHits(hits)]$target] # prefer assignment to target
    hits <- hits[!duplicated(queryHits(hits))] # if snp close to two targets pick one
    bad=ifelse(snp[queryHits(hits)]$t.AF<0.5,
               round(   snp[queryHits(hits)]$t.AF  * snp[queryHits(hits)]$t.DP),
               round((1-snp[queryHits(hits)]$t.AF) * snp[queryHits(hits)]$t.DP))
    tmp <- data.table(
        idx=subjectHits(hits),
        bad=bad,
        depth=snp[queryHits(hits)]$t.DP
    )
    setkey(tmp, idx)
    tmp <- tmp[J(1:length(gt))]
    tmp <- tmp[,.(baf=sum(bad)/sum(depth)),by=idx]
    gt$baf <- tmp$baf
    gt$baf <- .smooth.outliers.gr(gt, "baf")
    return(gt)
}

getBaf <- function(cnv, opts) {
    if (opts$target=="genome") {
        target.fun <- .getGenomeBaf
    } else {
        target.fun <- .getTargetBaf
    }
    snp <- cnv$var[cnv$var$germline]
    cnv$tile <- target.fun(cnv$tile, snp, opts)
    return(cnv)
}
