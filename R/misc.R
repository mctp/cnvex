.robust.import <- function(fn, seqi, skip.chr=NULL) {
    tmp <- import(fn)
    valid.seql <- setdiff(intersect(seqlevels(tmp), seqlevels(seqi)), skip.chr)
    tmp <- keepSeqlevels(tmp, valid.seql, pruning.mode="coarse")
    seqlevels(tmp) <- seqlevels(seqi)
    seqinfo(tmp) <- seqi
    return(tmp)
}

.log.sum.exp <- function(x) {
    offset <- max(x)
    log(sum(exp(x - offset))) + offset
}

.dtRound <- function(dt) {
    for (i in names(dt)) {
        if (class(dt[,get(i)])=="numeric") {
            dt[,":="((i), round(get(i), 6))]
        }
    }
    return(dt)
}

max.na.rm <- function(x) max(x, na.rm=TRUE)

merge.list <- function (x, y) {
    if (length(x) == 0) 
        return(y)
    if (length(y) == 0) 
        return(x)
    i = match(names(y), names(x))
    i = is.na(i)
    if (any(i)) 
        x[names(y)[which(i)]] = y[which(i)]
    return(x)
}

.absMedDiff <- function(x, y) {
    abd <- abs(median(x, na.rm=TRUE)-median(y, na.rm=TRUE))
    return(abd)
}

