## data generation
lrData <- function(cnv, opts) {
    tile <- cnv$tile
    segs <- cnv$seg
    local.sd <- opts$local.sd
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
    if (local.sd) {
        tmp[,sd := sd(lr, na.rm=TRUE), by=seg]
        tmp[,sd := ifelse(n.lr>30, sd, global.sd)]
    } else {
        tmp[,sd := global.sd]
    }
    return(tmp)
}
