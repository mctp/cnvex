.detect.sex <- function(var, tile) {
    sel <- var$n.GT %in% c("0/1", "1/0") & var$n.DP>16
    chrx.snp <- length(var[sel & seqnames(var)=="chrX"])
    chr2.snp <- length(var[sel & seqnames(var)=="chr2"])
    sex.snp <-  if (chrx.snp / chr2.snp > 0.25) "female" else "male"
    return(sex.snp)
}

.detect.offset <- function(cnv, opts) {
    tmp <- as.data.table(mcols(cnv$tile)[,c("seg", "lr.smooth", "baf")])
    tmp <- tmp[,.(.N, lr=mean(lr.smooth, na.rm=TRUE), baf=mean(baf, na.rm=TRUE)), by=seg]
    tmp <- tmp[is.finite(lr) & is.finite(baf) & baf>0.38]
    tmp <- tmp[order(lr)][1:(nrow(tmp)%/%2)]
    tmp <- tmp[order(-N)][1:(nrow(tmp)%/%1.5)]
    offset <- median(tmp$lr)
    return(offset)
}
