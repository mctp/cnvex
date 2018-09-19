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

.poolCoverage <- function(pool, cnv, opts) {
    sex <- .detect.sex(cnv$var, cnv$tile)
    pool.sex <- pool$cov[,pool$sex==sex]
    ## rank normals by variance in log2(tumor/normal)
    lr.pool.mx <- log2(cnv$tile$t.cov/pool.sex)
    lr.pool.mx <- ifelse(is.finite(lr.pool.mx), lr.pool.mx, NA_real_)
    lr.pool.sd <- apply(lr.pool.mx, 2, estimateSd)
    lr.pool.rk <- order(lr.pool.mx.sd)
    ## greedy search for k-best normals to pool
    k <- which.min(sapply(seq_len(min(25, ncol(pool.sex))), function(i) {
        n <- lr.pool.rk[1:i]
        n.pool <- pool.sex[,n,drop=FALSE]
        n.cov <- .normCoverage(rowMedians(n.pool), cnv$tile, opts)
        lr.pool <- log2(cnv$tile$t.cov/n.cov)
        lr.pool <- ifelse(is.finite(lr.pool), lr.pool, NA_real_)
        sd.pool <- estimateSd(lr.pool)
    }))
    ## compute final log-ratio
    n <- lr.pool.rk[1:k]
    n.pool <- pool.sex[,n,drop=FALSE]
    n.cov <- .normCoverage(rowMedians(n.pool), cnv$tile, opts)
    return(n.cov)
}
