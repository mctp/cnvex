bafData <- function(cnv) {
    tmp <- as.data.table(mcols(cnv$snp)[,c("t.AF", "t.DP", "seg")])
    tmp[,idx:=1:nrow(tmp)]
    return(tmp)
}
