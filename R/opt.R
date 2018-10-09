## data generation
covData <- function(cnv, opts) {
    tile <- cnv$tile
    segs <- cnv$seg
    tmp <- data.table(
        seg=tile$seg,
        lr=tile$lr,
        nC=tile$nC
    )
    ##
    tmp[,":="(
        len=width(segs)[seg]/.N
    ),by=seg]
    global.sd <- estimateSd(tmp$lr)
    if (opts$opt.local.sd) {
        tmp[,sd := sd(lr, na.rm=TRUE), by=seg]
        tmp[,sd := ifelse(n.lr>30, sd, global.sd)]
    } else {
        tmp[,sd := global.sd]
    }
    return(tmp[])
}

bafData <- function(cnv, opts) {
    tmp <- as.data.table(mcols(cnv$var)[,c("t.AF", "t.DP", "seg")])
    tmp[,idx:=1:nrow(tmp)]
    germline <- filterGermlineHets(cnv$tile, cnv$var, opts)
    tmp <- tmp[germline]
    return(tmp)
}

covGrid <- function(data, opts) {
    rC <- .lr.grid.rC(data$seg, data$lr, data$sd, data$len, data$nC, opts$max.C)
    pD <- .lr.grid.pD(opts$grid.n, opts$p.lo, opts$p.hi, opts$D.lo, opts$D.hi)
    grid <- rbindlist(mclapply(seq_len(nrow(pD)), function(i) {
        .llik.rC.p.D.wrap(rC, pD[i,1], pD[i,2], opts$max.sC, opts$max.len.per.probe)
    }, mc.cores=detectCores()))
    return(grid)
}
