.geneLogRatio <- function(cnv, opts) {
    ## data
    off <- .detect.offset(cnv, opts)
    tile <- cnv$tile
    gene <- cnv$gene[cnv$gene %over% tile]
    gene.data <- as.data.table(mcols(gene))
    tile.data <- as.data.table(mcols(tile))

    ## segment
    tmp.seg <- tile.data[,.(seg, lr.smooth, baf)]
    tmp.seg <- tmp.seg[,.(.N, lr.seg=mean(lr.smooth-off, na.rm=TRUE), mzd.seg=mean(baf, na.rm=TRUE)),
                       by=seg]
    tmp.seg <- tmp.seg[gene$seg, .(n.seg=N, lr.seg, mzd.seg)]

    ## local
    tmp.loc <- as.data.table(findOverlaps(gene, tile, maxgap = opts$tileflank * opts$tile.width))
    tmp.loc <- cbind(tmp.loc,
                     tile.data[tmp.loc$subjectHits, .(lr.smooth, baf)]
                     )
    tmp.loc <- tmp.loc[,.(n.loc=.N, lr.loc=mean(lr.smooth-off, na.rm=TRUE), mzd.loc=mean(baf, na.rm=TRUE)),
                       by=.(queryHits)]

    ## combine
    tmp <- cbind(gene.data[,.(gene_name, gene_id, seg)],
                 tmp.seg[,.(n.seg, lr.seg, mzd.seg)],
                 tmp.loc[,.(n.loc, lr.loc, mzd.loc)])
    return(tmp)
}
