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

.llik.rC.p.D <- function(rC, p, D, max.sC, max.len.per.probe) {
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
    x <- .llik.rC.p.D(rC, p, D, max.sC, max.len.per.probe)
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

.llik.baf.MC <- function(snpt, MC) {
    tmp1 <- rbindlist(rep(list(snpt), nrow(MC)))[order(idx)]
    tmp2 <- rbindlist(rep(list(MC), nrow(snpt)))
    x <- cbind(tmp1, tmp2)
    x[,
      beta:=dbeta(
          Ef,
          shape1=   t.AF  * t.DP + 1,
          shape2=(1-t.AF) * t.DP + 1,
          log=TRUE
      )]
    return(x)
}
