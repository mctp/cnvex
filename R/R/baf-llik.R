.priorC <- function(C, D, e=1) {
    a <- 1/(abs((0:7)-D)+e)
    b <- log(a/sum(a))
    b[C+1]
}

.priorK <- function(C, e=1) {
    log(1/(floor((C)/2) + e))
}

.beta.grid <- function(p, D, Cmax) {
    beta.grid <- as.data.table(expand.grid(M=0:Cmax,C=0:Cmax))[M<=C]
    beta.grid[,":="(
        Ef=(p * M + 1 * (1-p)) / (p * C + 2 * (1-p)),
        K=pmin(M,C-M),
        PC=.priorC(C, D),
        PK=.priorK(C)
    )]
    return(beta.grid)
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
