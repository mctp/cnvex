.normCoverage <- function(cov, gt, opts) {
    ## normalized coverage
    cov.n <- cov*width(gt)
    cov.n.tar <- cov.n[ gt$target]
    cov.n.tar <- cov.n.tar/(sum(cov.n.tar)/1e6)
    cov.n.tar <- 1000*cov.n.tar/width(gt[ gt$target])
    cov.n.off <- cov.n[!gt$target]
    cov.n.off <- cov.n.off/(sum(cov.n.off)/1e6)
    cov.n.off <- 1000*cov.n.off/width(gt[!gt$target])
    cov.n[ gt$target] <- cov.n.tar
    cov.n[!gt$target] <- cov.n.off
    return(cov.n)
}
