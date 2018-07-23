bafData <- function(cnv, opts) {
    tmp <- as.data.table(mcols(cnv$var)[,c("t.AF", "t.DP", "seg")])
    tmp[,idx:=1:nrow(tmp)]
    return(tmp[])
}
