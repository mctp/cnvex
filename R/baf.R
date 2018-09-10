.getTargetHits <- function(gt, snp, opts) {
    hits <- findOverlaps(snp, gt, maxgap = opts$tile.shoulder-1)
    hits <- hits[gt[subjectHits(hits)]$target] # prefer assignment to target
    hits <- hits[!duplicated(queryHits(hits))] # if snp close to two targets pick one
    return(hits)
}

.getGenomeHits <- function(gt, snp, opts) {
    hits <- findOverlaps(snp, gt, maxgap = opts$tile.shoulder-1)
    return(hits)
}

.filterGenomeGermlineHets <- function(gt, var, opts) {
    min.qual <- quantile(var$QUAL, 0.25, na.rm=TRUE)
    t.hi.dp <- quantile(var$t.DP, 0.95)
    t.lo.dp <- quantile(var$t.DP, 0.05)
    n.hi.dp <- quantile(var$n.DP, 0.95)
    n.lo.dp <- quantile(var$n.DP, 0.05)
    var$germline <- 
        ## germline
        !var$SOMATIC &
        ## simple variant
        var$TYPE=="SNV" &
        ## heterozygous
        var$n.GT %in% c("0/1", "1/0") &
        ## right range
        var$n.AF > 0.25 & var$n.AF < 0.75 &
        ## high-quality
        var$t.DP > t.lo.dp & var$t.DP < t.hi.dp &
        var$n.DP > n.lo.dp & var$n.DP < n.hi.dp &
        (!var$mask.strict | var$QUAL > min.qual)
    return(var)
}

.filterTargetGermlineHets <- function(gt, var, opts) {
    splash <- gt[gte$target] + opts$shoulder
    var$germline <- 
        var %over% splash &
        ## germline
        !var$SOMATIC &
        ## simple variant
        var$TYPE=="SNV" &
        ## heterozygous
        var$n.GT %in% c("0/1", "1/0") &
        ## right range
        var$n.AF > 0.25 & var$n.AF < 0.75 &
        ## high 
        var$n.DP > 16 &
        var$t.DP > 16 &
        ## mappability
        (!var$mask.strict)
    return(var)
}

.filterSnp <- function(tile, var, opts) {
    if (opts$target=="genome") {
        snp.filter.fun <- .filterGenomeGermlineHets
    } else {
        snp.filter.fun <- .filterTargetGermlineHets
    }
    var <- snp.filter.fun(gt, var, opts)
    snp <- var[var$germline]
    return(snp)
}

.getBaf <- function(gt, var, opts) {
    ## filter variants to SNPs
    snp <- .filterSnp(gt, var, opts)
    if (opts$target=="genome") {
        hits <- .getGenomeHits(gt, snp, opts)
    } else {
        hits <- .getTargetHits(gt, snp, opts)
    }
    bad=ifelse(snp[queryHits(hits)]$t.AF < 0.5,
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
    cnv$tile <- .getBaf(cnv$tile, cnv$var, opts)
    return(cnv)
}
