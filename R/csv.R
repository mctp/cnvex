.geneCSV <- function(cnv, opts) {
    ## data
    tile <- cnv$tile
    gene <- cnv$gene[cnv$gene %over% tile]
    gene.data <- as.data.table(mcols(gene))
    tile.data <- as.data.table(mcols(tile))

    ## segment
    tmp.seg <- tile.data[,.(seg, lr.smooth, baf)]
    tmp.seg <- tmp.seg[,.(.N, lr.seg=mean(lr.smooth, na.rm=TRUE), mzd.seg=mean(baf, na.rm=TRUE)),
                       by=seg]
    tmp.seg <- tmp.seg[gene$seg, .(n.seg=N, lr.seg, mzd.seg)]

    ## local
    tmp.loc <- as.data.table(findOverlaps(gene, tile, maxgap = opts$tileflank * opts$tile.width))
    tmp.loc <- cbind(tmp.loc,
                     tile.data[tmp.loc$subjectHits, .(lr.smooth, baf)]
                     )
    tmp.loc <- tmp.loc[,.(n.loc=.N, lr.loc=mean(lr.smooth, na.rm=TRUE), mzd.loc=mean(baf, na.rm=TRUE)),
                       by=.(queryHits)]

    ## combine
    tmp <- cbind(gene.data[,.(gene_name, gene_id, seg)],
                 tmp.seg[,.(n.seg, lr.seg, mzd.seg)],
                 tmp.loc[,.(n.loc, lr.loc, mzd.loc)])
    return(tmp)
}

.segmentCSV <- function(cnv, opts) {
    tmp1 <- as.data.table(mcols(cnv$tile))
    tmp1 <- tmp1[,.(lr=mean(lr.smooth, na.rm=TRUE), mzd=mean(baf, na.rm=TRUE), .N),by=seg]
    tmp2 <- as.data.table(cnv$seg)[,.(chr=seqnames, start, end, seg=.I)]
    setkey(tmp1, seg)
    setkey(tmp2, seg)
    seg.stat <- tmp1[tmp2,.(chr, start, end, N, lr, mzd)]
    return(seg.stat)
}
