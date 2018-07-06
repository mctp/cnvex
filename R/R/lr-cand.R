getCand <- function(grid, opts) {
    x <- as.matrix(dcast.data.table(grid, D0~p0, value.var="L1")[,-1])
    r <- raster(x)
    ## nearest odd integer >= to 3
    rres <- max(2*floor((nrow(r)*opts$res)/2)+1, 3)
    cres <- max(2*floor((ncol(r)*opts$res)/2)+1, 3)
    wind <- matrix(1, nrow=rres, ncol=cres)
    localmax <- focal(r, fun = max.na.rm, w = wind, pad=TRUE, padValue=NA)
    cand <- grid[Which(localmax==r, cells=TRUE)]
    return(cand)
}

fixCand <- function(data, cand, grid.n=11, p.offset=0.05, lr.offset=0.2, max.C=9, max.sC=20, max.len.per.probe=1e6) {
    rC <- .lr.grid.rC(data$seg, data$lr, data$sd, data$len, data$nC, max.C)
    cand.opt <- rbindlist(lapply(seq_len(nrow(cand)), function(i) {
        #### pD candidate
        p0 <- cand[i,p0]
        D1 <- cand[i,D1]
        #### grid
        p.scan <- seq(-p.offset, p.offset, length.out=grid.n)
        lr.scan <- seq(-lr.offset, lr.offset, length.out=grid.n)
        off.grid <- expand.grid(lr=lr.scan, p=p.scan)
        tmp <- rbindlist(mclapply(seq_len(nrow(off.grid)), function(ioff) {
            p.fix <- off.grid[ioff,"p"]
            lr.fix <- off.grid[ioff,"lr"]
            rCO <- copy(rC)
            rCO[,lr:=lr+lr.fix]
            ret <- .llik.rC.p.D.wrap(rCO, p0+p.fix, D1, max.sC, max.len.per.probe)
            ret$p.fix <- p.fix
            ret$lr.fix <- lr.fix
            return(ret)
        }, mc.cores=detectCores()))
        tmp[order(-L1),.SD[1]]
    }))
    return(cand.opt)
}

optCand <- function(data, cand, grid.n=32, p.offset=0.025, max.iter=50, max.C=9, 
                    max.sC=20, max.len.per.probe=1e6
                     ) {
    rC <- .lr.grid.rC(data$seg, data$lr, data$sd, data$len, data$nC, max.C)
    ## optimize p
    rets <- rbindlist(mclapply(seq_len(nrow(cand)), function(i) {
        p.fix <- cand[i,p.fix]
        lr.fix <- cand[i,lr.fix]
        p0 <- cand[i,p0]
        D0 <- cand[i,D0]
        D1 <- cand[i,D1]
        S1 <- cand[i,S1]
        ## setup variables
        iter <- 0
        rCO <- copy(rC)
        rCO[,lr:=lr+lr.fix]
        pi <- p0 + p.fix
        Di <- D1
        Si <- S1
        ps <- c()
        Ds <- c()
        ## iterate till convergence
        status <- "running"
        while(TRUE) {
            iter <- iter + 1
            ## optimize p
            p.lo <- max(pi-p.offset, 0.01)
            p.hi <- min(pi+p.offset, 0.99)
            pD <- cbind(p=seq(p.lo, p.hi, length.out=grid.n), D=Di)
            tmp <- rbindlist(lapply(seq_len(nrow(pD)), function(i) {
                .llik.rC.p.D.wrap(rCO, pD[i,1], pD[i,2], max.sC, max.len.per.probe)
            }))
            optpD <- tmp[order(-L1), .SD[1]]
            pj <- optpD[,p0]
            Dj <- optpD[,D1]
            Sj <- optpD[,S1]
            if (((abs(pj-pi)<5e-3) && (abs(Dj-Di)<5e-2)) ||
                ((pj %in% ps) && (Dj %in% Ds))) {
                status <- "converged"
            } else if (pj<=0.01+1e-9 || pj>=0.99-1e-9 || Dj<1+1e-9 || Dj>6-1e-9) {
                status <- "diverged"
            } else if (iter==max.iter) {
                status <- "maxiter"
            } else if (Sj-S1>0.02) {
                status <- "aborted"
            }
            ps <- c(ps, pj)
            Ds <- c(Ds, Dj)
            pi <- pj
            Di <- Dj
            cat(sprintf("candidate:%s iter:%s p:%.4f D:%.4f llik:%.0f status:%s\n",
                        i, iter, pi, Di, optpD[,L1], status))
            if (status != "running") {
                break
            }
        }
        ret <- list(p0=p0, D0=D0, D1=D1, L1=cand[i,L1], S1=cand[i,S1], pi=pi, Di=Di,
                    Li=optpD[,L1], Si=optpD[,S1], lr.fix=lr.fix, iter=iter, status=status)
        return(ret)
    }, mc.cores=detectCores()))
    return(rets)
}
