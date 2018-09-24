.coverageFilterSel <- function(xs, lo, hi) {
    uxs <- rowMeans(xs)
    muxs <- median(uxs)
    f <- lo * muxs < uxs & uxs < hi * muxs
    return(f)
}

.zeroFilterSel <- function(xs, hi) {
    f <- rowMeans(xs==0) < hi & rowMads(xs)!=0
    return(f)
}

.coverageFilter <- function(x, s, lo, hi) {
    f <- logical(length(s))
    f[ s] <- .coverageFilterSel(x[ s,], lo, hi)
    f[!s] <- .coverageFilterSel(x[!s,], lo, hi)
    return(f)
}

.zeroFilter <- function(x, s, hi) {
    f <- logical(length(s))
    f[ s] <- .zeroFilterSel(x[ s,], hi)
    f[!s] <- .zeroFilterSel(x[!s,], hi)
    return(f)
}

.poolFilter <- function(cov, sex, target, lo.cov, hi.cov, hi.zero) {
    cfm <- .coverageFilter(cov[,sex==  "male"], target, lo.cov, hi.cov)
    cff <- .coverageFilter(cov[,sex=="female"], target, lo.cov, hi.cov)
    zfm <- .zeroFilter(cov[,sex==  "male"], target, hi.zero)
    zff <- .zeroFilter(cov[,sex=="female"], target, hi.zero)
    f <- (cfm | cff) & (zfm | zff)
    return(f)
}

.outlierMaskSel <- function(xs, sd) {
    ## identify outlier tiles
    xs.med <- rowMedians(xs)
    xs.mad <- rowMads(xs)
    xs.out <- (xs - xs.med) / xs.mad
    ## lo-threshold
    xs[xs.out < -sd] <- NA_real_
    xs.min <- rowMins(xs, na.rm = TRUE)
    xs <- apply(xs, 2, function(col) {
        col[is.na(col)] <- xs.min[is.na(col)]
        col
    })
    ## hi-threshold
    xs[xs.out > sd] <- NA_real_
    xs.max <- rowMaxs(xs, na.rm = TRUE)
    xs <- apply(xs, 2, function(col) {
        col[is.na(col)] <- xs.max[is.na(col)]
        col
    })
    return(xs)
}

.outlierMask <- function(cov, sex, filter, sd) {
    cov[ filter,sex==  "male"] <- .outlierMaskSel(cov[filter,sex==  "male"], sd)
    cov[ filter,sex=="female"] <- .outlierMaskSel(cov[filter,sex=="female"], sd)
    cov[!filter,] <- NA_real_
    return(cov)
}

.medianLogRatioSel <- function(xs, sex, sex.chr) {
    xsr <- xs
    ## autosomes
    xsa <- xs[!sex.chr,]
    xsam <- rowMedians(xsa)
    xsr[!sex.chr,] <- log2(xsa/xsam)
    ## sex chromosomes
    xsm <- xs[sex.chr, sex=="male"]
    xsmm <- rowMedians(xsm)
    xsr[sex.chr,sex=="male"] <- log2(xsm/xsmm)
    xsf <- xs[sex.chr, sex=="female"]
    xsfm <- rowMedians(xsf)
    xsr[sex.chr,sex=="female"] <- log2(xsf/xsfm)
    ## fix bounds
    xsrm <- max(min(xsr[is.finite(xsr)]), -4)
    xsrx <- min(max(xsr[is.finite(xsr)]),  4)
    xsr[is.nan(xsr)] <- 0    ## 0 / 0
    xsr[xsr < xsrm] <-  xsrm ## 0 / x
    xsr[xsr > xsrx] <-  xsrx ## x / 0
    return(xsr)
}

.medianLogRatio <- function(x, s, f, sex, sex.chr) {
    x[( s & f),] <- .medianLogRatioSel(x[( s & f),], sex, sex.chr[( s & f)])
    x[(!s & f),] <- .medianLogRatioSel(x[(!s & f),], sex, sex.chr[(!s & f)])
    x[f,] <- t(t(x[f,]) - colMedians(x[f,]))
    return(x)
}

.medianCoverage <- function(cov, sex, sex.chr) {
    m <- rowMedians(cov)
    mm <- rowMedians(cov[,sex==  "male"])
    mf <- rowMedians(cov[,sex=="female"])
    tile.median <- cbind(
          male=ifelse(sex.chr, mm, m),
        female=ifelse(sex.chr, mf, m)
    )
    return(tile.median)
}

.importPoolData <- function(cnv.fns, opts) {
    cnvs <- mclapply(cnv.fns, function(fn) {
        cnv <- readRDS(fn)
    }, mc.cores=opts$cores)
    cov <- do.call(cbind, lapply(cnvs, function(cnv) cnv$tile$n.cov))
    sex <- sapply(cnvs, function(cnv) .detect.sex(cnv$var, cnv$tile))
    sex.chr <- as.logical(as.character(seqnames(cnvs[[1]]$tile)) %in% c("chrX", "chrY"))
    target <- cnvs[[1]]$tile$target
    pd <- list(cov=cov, sex=sex, sex.chr=sex.chr, target=target)
    return(pd)
}

.createPool <- function(pd, opts) {
    filter <- .poolFilter(pd$cov, pd$sex, pd$target, opts$pool.lo.cov, opts$pool.hi.cov, opts$pool.hi.zero)
    cov <- .outlierMask(pd$cov, pd$sex, filter, opts$pool.sd.out)
    lr <- .medianLogRatio(cov, pd$target, filter, pd$sex, pd$sex.chr)
    mc <- .medianCoverage(cov, pd$sex, pd$sex.chr)
    if (opts$pool.method=="pca") {
        p <- svd(lr[filter,])$u
    }
    if (opts$pool.method=="ica") {
        p <- fastICA(lr[filter,], 40, method="C")$S
    }
    pool <- list(method=opts$pool.method, sex=pd$sex, target=pd$target, sex.chr=pd$sex.chr,
                 cov=cov, cov.med=mc, projection=p, filter=filter)
    return(pool)
}

.fixLogRatio <- function(lr, filter, opts) {
    lr[!filter] <- NA_real_
    rf <- lr[filter]
    rf[is.nan(rf)] <- 0 ## 0 / 0
    rf[rf < -4] <- -4 ## 0 / x
    rf[rf >  4] <-  4 ## x / 0
    return(rf)
}

.poolIcaDenoise <- function(cnv, pool, opts) {
    sex <- .detect.sex(cnv$var, cnv$tile)
    lr <- log2(cnv$tile$t.cov / pool$cov.med[,sex])
    rf <- .fixLogRatio(lr, pool$filter, opts)
    S <- pool$projection[,seq(1,opts$pool.n.comp)]
    lr[pool$filter] <- rf - as.vector(tcrossprod(rf %*% t(ginv(S)), S))
    return(lr)
}

.poolPcaDenoise <- function(cnv, pool, opts) {
    sex <- .detect.sex(cnv$var, cnv$tile)
    lr <- log2(cnv$tile$t.cov / pool$cov.med[,sex])
    rf <- .fixLogRatio(lr, pool$filter, opts)
    P <- pool$projection[,seq(1,opts$pool.n.comp)]
    lr[pool$filter] <- rf - as.vector(tcrossprod(rf %*% P, P))
    return(lr)
}

.poolFakeCoverage <- function(cnv, pool, opts) {
    sex <- .detect.sex(cnv$var, cnv$tile)
    pool.sex <- pool$cov[,pool$sex==sex]
    ## rank normals by variance in log2(tumor/normal)
    lr.pool.mx <- log2(cnv$tile$t.cov/pool.sex)
    lr.pool.mx <- ifelse(is.finite(lr.pool.mx), lr.pool.mx, NA_real_)
    lr.pool.sd <- apply(lr.pool.mx, 2, estimateSd)
    lr.pool.rk <- order(lr.pool.sd)
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
