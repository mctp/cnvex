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
    min.qual <- quantile(var$QUAL, opts$baf.qual.min, na.rm=TRUE)
    t.hi.dp <- quantile(var$t.DP, 0.95)
    t.lo.dp <- quantile(var$t.DP, 0.05)
    n.hi.dp <- quantile(var$n.DP, 0.95)
    n.lo.dp <- quantile(var$n.DP, 0.05)
    germline <- 
        ## germline
        !var$SOMATIC &
        ## simple variant
        var$TYPE=="SNV" &
        ## heterozygous
        var$n.GT %in% c("0/1", "1/0") &
        var$n.AF > opts$baf.het.range[1] & var$n.AF < opts$baf.het.range[2] &
        ## sufficient coverage
        var$t.DP > opts$baf.min.t.dp &
        ## high-quality
        var$t.DP > t.lo.dp & var$t.DP < t.hi.dp &
        var$n.DP > n.lo.dp & var$n.DP < n.hi.dp &
        (!var$mask.loose & (!var$mask.strict | var$QUAL > min.qual))
    return(germline)
}

.filterTargetGermlineHets <- function(gt, var, opts) {
    min.qual <- quantile(var$QUAL, opts$baf.qual.min, na.rm=TRUE)
    t.hi.dp <- quantile(var$t.DP, 0.95)
    t.lo.dp <- quantile(var$t.DP, 0.05)
    n.hi.dp <- quantile(var$n.DP, 0.95)
    n.lo.dp <- quantile(var$n.DP, 0.05)
    splash <- gt[gt$target] + opts$tile.shoulder
    germline <-
        var %over% splash &
        ## germline
        !var$SOMATIC &
        ## simple variant
        var$TYPE=="SNV" &
        ## heterozygous
        var$n.GT %in% c("0/1", "1/0") &
        ## right range
        var$n.AF > opts$baf.het.range[1] & var$n.AF < opts$baf.het.range[2] &
        ## high 
        var$n.DP > opts$baf.min.n.dp &
        var$t.DP > opts$baf.min.t.dp &
        ## mappability
        var$t.DP > t.lo.dp & var$t.DP < t.hi.dp &
        var$n.DP > n.lo.dp & var$n.DP < n.hi.dp &
        (!var$mask.loose & (!var$mask.strict | var$QUAL > min.qual))
    return(germline)
}

.filterGermlineHets <- function(tile, var, opts) {
    var.genome <- var[var$SOURCE=="genome"]
    germline.genome <- .filterGenomeGermlineHets(tile, var.genome, opts)
    var.target <- var[var$SOURCE=="target"]
    germline.target <- .filterTargetGermlineHets(tile, var.target, opts)
    germline <- ifelse(var$SOURCE=="genome", germline.genome, germline.target)
    return(germline)
}

.getBaf <- function(gt, var, opts) {
    ## filter variants to SNPs
    germline <- .filterGermlineHets(gt, var, opts)
    snp <- var[germline]
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
    tmp <- tmp[,.(
        baf=sum(bad)/sum(depth),
        depth=sum(depth)
    ), by=idx]
    tmp <- tmp[,":="(
        ## weight proportional to inverse varince ~ 1/sqrt(n)
        baf.weight = 1/sqrt(0.25/pmin(pmax(1, depth), opts$baf.max.eff.dp))
    )]
    gt$baf <- tmp$baf
    gt$baf.weight <- tmp$baf.weight
    gt$baf <- .smooth.outliers.gr(gt, "baf")
    return(gt)
}
