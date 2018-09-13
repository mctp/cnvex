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

## prior for C
.priorC <- function(C, D, Cmax, e=1) {
    a <- 1/(abs((0:7)-D)+e)
    b <- log(a/sum(a))
    b[C+1]
}

## prior for K
.priorK <- function(K, C, e=1) {
    log(1/(floor((C)/2) + e))
}

## MC grid generation
.baf.grid.MC <- function(p, D, Cmax) {
    beta.grid <- as.data.table(expand.grid(M=0:Cmax,C=0:Cmax))[M<=C]
    beta.grid[,":="(
        Ef=(p * M + 1 * (1-p)) / (p * C + 2 * (1-p)),
        K=pmin(M,C-M)
    )]
    ## priors
    beta.grid[,":="(
        PC=.priorC(C, D, Cmax),
        PK=.priorK(K, C)
    )]
    ## M - number of chromosome with variant
    ## C - total number of chromosomes
    ## Ef - expected frequency of variant
    ## PK - prior probability of K given C
    ## PC - prior probability of C given D
    return(beta.grid[])
}
