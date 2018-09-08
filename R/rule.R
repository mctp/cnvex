.detect.sex <- function(var, tile) {
    ## placeholder implementation of sex detection
    chrx.sel <- tile$target & (seqnames(tile) %in% c("chrX"))
    chrx.cov <- median(tile[chrx.sel]$n.cov)
    chry.sel <- tile$target & (seqnames(tile) %in% c("chrY"))
    chry.cov <- median(tile[chry.sel]$n.cov)
    sex.copy <- if (chrx.cov/chry.cov > 1.5) "female" else "male"
    chrx.snp <- length(var[var$germline & seqnames(var)=="chrX"])
    chr2.snp <- length(var[var$germline & seqnames(var)=="chr2"])
    sex.snp <-  if (chrx.snp / chr2.snp > 0.25) "female" else "male"
    if (sex.snp == sex.copy) {
        sex = sex.snp
    } else {
        sex = "unknown"
    }
    return(sex)
}

.detect.offset <- function(cnv, opts) {
    tmp <- as.data.table(mcols(cnv$tile)[,c("seg", "lr.smooth", "baf")])
    tmp <- tmp[,.(.N, lr=mean(lr.smooth, na.rm=TRUE), baf=mean(baf, na.rm=TRUE)), by=seg]
    tmp <- tmp[is.finite(lr) & is.finite(baf)]
    tmp <- tmp[order(abs(lr))][1:(nrow(tmp)%/%2)]
    tmp <- tmp[order(   -baf)][1:(nrow(tmp)%/%2)]
    offset <- median(tmp$lr)
    return(offset)
}
