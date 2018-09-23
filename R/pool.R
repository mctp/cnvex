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

.medianLogRatioSel <- function(xs) {
    xsm <- colMedians(xs)
    xsr <- log2(t(t(xs) / xsm))
    xsrx <- max(xsr[is.finite(xsr)])
    xsr[is.nan(xsr)] <- 0     ## 0 / 0
    xsr[xsr < -xsrx] <- -xsrx ## 0 / x
    xsr[xsr >  xsrx] <-  xsrx ## x / 0
    return(xsr)
}

.poolFilter <- function(x, s, lo.cov=0.20, hi.cov=5, hi.zero=0.25) {
    cf <- .coverageFilter(x, s, lo.cov, hi.cov)
    zf <- .zeroFilter(x, s, hi.zero)
    f <- cf & zf
    return(f)
}

.outlierMask <- function(x, s, sd) {
    x[ s,] <- .outlierMaskSel(x[ s,], sd)
    x[!s,] <- NA_real_
    return(x)
}

.medianLogRatio <- function(x, s, f) {
    x[( s & f),] <- .medianLogRatioSel(x[( s & f),])
    x[(!s & f),] <- .medianLogRatioSel(x[(!s & f),])
    return(x)
}


.createPool <- function(cov, target, opts) {
    f <- .poolFilter(cov, target, opts$pool.lo.cov, opts$pool.hi.cov, opts$pool.hi.zero)
    m <- .outlierMask(cov, f, opts$pool.sd.out)
    r <- .medianLogRatio(m, s, f)
    if (opts$pool.method=="pca") {
        p <- svd(r[f,])$u
    }
    if (opts$pool.method=="ica") {
    }
    out <- list(projection=p, filter=f)
    return(out)
}
