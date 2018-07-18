.detect.sex <- function(var, tile) {
    ## placeholder implementation of sex detection
    chrx.sel <- tile$target & (seqnames(tile) %in% c("chrX"))
    chrx.cov <- median(tile[chrx.sel]$n.cov)
    chry.sel <- tile$target & (seqnames(tile) %in% c("chrY"))
    chry.cov <- median(tile[chry.sel]$n.cov)
    sex <- if (chrx.cov/chry.cov > 1.5) "female" else "male"
    return(sex)
}

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

.addCopy <- function(gr, sex) {
    gr$nC <- .getCopy(gr, sex)
    return(gr)
}

.addBaf <- function(gt, snp, opts) {
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

## ## import all data
## tile <- target.fun(target.fn, opts)
## tile <- .addCoverage(tile, t.bam, n.bam, opts)
## ## detect sex
## if (is.null(opts$sex)) {
##     sex <- .detect.sex(var, tile)
## } else {
##     sex <- opts$sex
## }
## tile <- .addCopy(tile, sex)
## var <- .addCopy(var, sex)
## ## add BAF
## snp <- snp.filter.fun(var, tile, opts)
## tile <- .addBaf(tile, snp, opts)
