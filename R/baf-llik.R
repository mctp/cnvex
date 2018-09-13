.priorC <- function(C, D, Cmax, e=1) {
    a <- 1/(abs((0:7)-D)+e)
    b <- log(a/sum(a))
    b[C+1]
}

.priorK <- function(K, C, e=1) {
    log(1/(floor((C)/2) + e))
}

.beta.grid <- function(p, D, Cmax) {
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
    ## M	C	Ef	K	PC	PK
    ## minor	total	?	?	?	?
    ## M - number of chromosome with variant
    ## C - total number of chromosomes
    ## Ef - expected frequency of variant
    ## PK - prior probability of K given C
    ## PC - prior probability of C given D
    return(beta.grid[])
}

.baf.grid.inner <- function(snpt, beta.grid) {
    tmp1 <- rbindlist(rep(list(snpt), nrow(beta.grid)))[order(idx)]
    tmp2 <- rbindlist(rep(list(beta.grid), nrow(snpt)))
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
