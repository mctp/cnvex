.addBaf <- function(gt, snp) {
    hits <- findOverlaps(gt, snp)
    bad=ifelse(snp[subjectHits(hits)]$t.AF<0.5,
               round(   snp[subjectHits(hits)]$t.AF  * snp[subjectHits(hits)]$t.DP),
               round((1-snp[subjectHits(hits)]$t.AF) * snp[subjectHits(hits)]$t.DP))
    tmp <- data.table(
        idx=queryHits(hits),
        bad=bad,
        depth=snp[subjectHits(hits)]$t.DP
    )
    setkey(tmp, idx)
    tmp <- tmp[J(1:length(gt))]
    tmp <- tmp[,.(baf=sum(bad)/sum(depth)),by=idx]
    gt$baf <- tmp$baf
    return(gt)
}

.addJointSeg <- function(gt, ...) {
    gt <- unname(unlist(endoapply(split(gt, gt$arm), .jointSegArm, ...)))
    unname(sort(unlist(reduce(split(gt, paste(gt$arm, gt$seg))))))
    return(gt)
}

