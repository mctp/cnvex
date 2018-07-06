## lr likelihood function
.llik.rC.p.D.inner <- function(rC, p, D) {
    rC[,
       ":="(
           ## normal
           norm = dnorm(lr, mean=log2((p*C+(1-p)*nC)/D), sd=sd, log=TRUE),
           unif = dunif(lr, min=-6, max=6, log=TRUE)*2
       )
       ]
    x <- rC[,.( ## sum by C and seg
        lr = mean(lr, na.rm=TRUE),
        norm = sum(norm, na.rm=TRUE),
        unif = sum(unif, na.rm=TRUE),
        n = n[1],
        n.lr = n.lr[1],
        nC = nC[1],
        len = sum(len),
        p=p, D=D
    ), by=.(C, seg)]
    return(x)
}

.llik.rC.p.D.full <- function(rC, p, D, max.sC, max.len.per.probe) {
    ## compute likelihood for each segment and each C
    x <- .llik.rC.p.D.inner(rC, p, D)
    ## define sub-clonal segments
    x[,subc:=(unif>norm)]
    ## ML sub-clonal copy-number
    x[,sC  :=pmin((2^(lr) * D)/p - ((nC * (1 - p))/p), max.sC)]
    ## define high-quality segments based on coverage
    x[,hq:=is.finite(lr) & (len / n.lr < max.len.per.probe)]
    return(x)
}

.llik.rC.p.D.best <- function(rC, p, D, max.sC, max.len.per.probe) {
    x <- .llik.rC.p.D.full(rC, p, D, max.sC, max.len.per.probe)
    ## for each segment pick C with highest llik
    x <- x[order(-norm),.SD[1],by=seg]
    return(x)
}

.llik.rC.p.D.wrap <- function(rC, p, D, max.sC, max.len.per.probe) {
    x <- .llik.rC.p.D.best(rC, p, D, max.sC, max.len.per.probe)
    ## result
    y <- list(
        L1 = x[(hq), sum(ifelse(subc, unif, norm))], ## segment likelihood
        S1 = x[(hq), weighted.mean(subc, len)], ## proportion subclonal
        D1 = x[(hq),      p[1]  * weighted.mean(ifelse(subc, sC, C), len) +
                     (1 - p[1]) * weighted.mean(nC, len)], ## total ploidy
        p0 = p, D0 = D
    )
    return(y)
}


## .llik.full <- function(data, pi, Di, maxC) {
##     rC <- .grid.rC(data$seg, data$lr, data$sd, data$len, maxC)
##     CL <- rbindlist(mclapply(seq_len(length(pi)), function(i) {
##         .llik.full.inner(rC, pi[i], Di[i])
##     }, mc.cores=detectCores()))
##     seg.stats <- data[,.(lr=mean(lr,na.rm=TRUE), n.lr=n.lr[1]), by=seg]
##     setkey(seg.stats, seg)
##     setkey(CL, seg)
##     CL <- CL[seg.stats]
##     CL[,sC := (2^(lr) * Di)/pi - ((2 * (1 - pi))/pi)]
##     return(CL)
## }


## .llik.best <- function(data, pi, Di, maxC) {
##     rC <- .grid.rC(data$seg, data$lr, data$sd, data$len, maxC)
##     CL <- rbindlist(mclapply(seq_len(length(pi)), function(i) {
##         .llik.best.inner(rC, pi[i], Di[i])
##     }, mc.cores=detectCores()))
##     seg.stats <- data[,.(lr=mean(lr, na.rm=TRUE), n.lr=n.lr[1]), by=seg]
##     setkey(seg.stats, seg)
##     setkey(CL, seg)
##     CL <- CL[seg.stats]
##     CL[,sC := (2^(lr) * Di)/pi - ((2 * (1 - pi))/pi)]
##     return(CL)
## }

## ##
## .cand.lr <- function(cand, data) {
##     seg <- data[,.(lr=mean(lr, na.rm=TRUE),
##                    n.lr=seg.n[1]
##                    ), by=seg]
##     setkey(cand, seg)
##     setkey(seg, seg)
##     cand <- cand[seg]
##     return(cand)
## }
