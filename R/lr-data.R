## data generation
lrData <- function(cnv, opts) {
    tile <- cnv$tile
    segs <- cnv$seg
    tmp <- data.table(
        seg=tile$seg,
        lr=tile$lr.smooth,
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
