## pD grid generation
.lr.grid.pD <- function(grid.n, p.lo, p.hi, D.lo, D.hi) {
    p <- seq(p.lo, p.hi, length.out=grid.n)
    D <- seq(D.lo, D.hi, length.out=grid.n)
    pD <- as.matrix(expand.grid(p=p, D=D))
    return(pD)
}

## rC grid generation
.lr.grid.rC <- function(seg, lr, sd, len, nC, max.C) {
    rC <- expand.grid(lr=lr, C=0:max.C)
    rC$seg <- seg
    rC$sd <- sd
    rC$len <- len
    rC$nC <- nC
    setDT(rC)
    setkey(rC, C, seg)
    rC[,":="(
        n=.N,
        n.lr=sum(!is.na(lr))
    ),by=seg]
    rC <- rC[nC>0] # skip Y if not present
    return(rC)
}

## inner likelihood function
lrGrid <- function(data, opts) {
    rC <- .lr.grid.rC(data$seg, data$lr, data$sd, data$len, data$nC, opts$max.C)
    pD <- .lr.grid.pD(opts$grid.n, opts$p.lo, opts$p.hi, opts$D.lo, opts$D.hi)
    grid <- rbindlist(mclapply(seq_len(nrow(pD)), function(i) {
        .llik.rC.p.D.wrap(rC, pD[i,1], pD[i,2], opts$max.sC, opts$max.len.per.probe)
    }, mc.cores=detectCores()))
    return(grid)
}
