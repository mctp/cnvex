## .normCoverage <- function(t.cov, n.cov, gt, opts) {
##     mx.cov <- cbind(t.cov, n.cov)*width(gt)
##     mx.cov.tar <- mx.cov[ gt$target,]
##     mx.cov.tar <- t(t(mx.cov.tar)/(colSums(mx.cov.tar)/1e6))
##     mx.cov.tar <- 1000*mx.cov.tar/width(gt[ gt$target])
##     mx.cov.off <- mx.cov[!gt$target,]
##     mx.cov.off <- t(t(mx.cov.off)/(colSums(mx.cov.off)/1e6))
##     mx.cov.off <- 1000*mx.cov.off/width(gt[!gt$target])
##     mx.cov[ gt$target,] <- mx.cov.tar
##     mx.cov[!gt$target,] <- mx.cov.off
##     return(mx.cov)
## }

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
    return(mx.cov)
}

.poolCoverage <- function(pool, cnv, opts) {
    cov <- cnv$tile$t.cov
    sex <- .detect.sex(cnv$var, cnv$tile)
    ## rank normals by variance in log2(tumor/normal)
    lr.pool.mx <- log2(cov/pool)
    lr.pool.mx <- ifelse(is.finite(lr.pool.mx), lr.pool.mx, NA_real_)
    lr.pool.mx.sd <- apply(lr.pool.mx, 2, estimateSd)
    lr.pool.mx.sd.rnk <- order(lr.pool.mx.sd)
    ## greedy search for k-best normals to pool
    k <- which.min(sapply(seq_len(min(25, ncol(pool))), function(i) {
        sel<- lr.pool.mx.sd.rnk[1:i]
        sel.pool <- pool[,sel,drop=FALSE]
        lr.pool <- log2(t.cov/rowMedians(sel.pool))
        lr.pool <- ifelse(is.finite(lr.pool), lr.pool, NA_real_)
        sd.pool <- estimateSd(lr.pool)
    }))
    p.cov <- rowMedians(pool[,lr.pool.mx.sd.rnk[1:k],drop=FALSE])
    p.cov <- .normCoverage(p.cov, cnv$tile, opts)
    return(p.cov)
}
