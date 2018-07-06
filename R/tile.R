.normCoverage <- function(t.cov, n.cov, gt, opts) {
    mx.cov <- cbind(t.cov, n.cov)*width(gt)
    mx.cov.tar <- mx.cov[ gt$target,]
    mx.cov.tar <- t(t(mx.cov.tar)/(colSums(mx.cov.tar)/1e6))
    mx.cov.tar <- 1000*mx.cov.tar/width(gt[ gt$target])
    gc.cov.tar <- cbind(gc=gt[gt$target]$gc, mx.cov.tar)
    mx.cov.off <- mx.cov[!gt$target,]
    mx.cov.off <- t(t(mx.cov.off)/(colSums(mx.cov.off)/1e6))
    mx.cov.off <- 1000*mx.cov.off/width(gt[!gt$target])
    mx.cov[ gt$target,] <- mx.cov.tar
    mx.cov[!gt$target,] <- mx.cov.off
    return(mx.cov)
}

.covarageRatio <- function(mx.cov, gt, opts) {
    lr.raw <- log2(mx.cov[,1]/mx.cov[,2])
    lr <- ifelse(is.finite(lr.raw), lr.raw, NA_real_)
    if (opts$gc.adjust.trend) {
        for (sel in c(TRUE, FALSE)) {
            lr.sel <- lr[gt$target==sel]
            gt.sel <- gt[gt$target==sel]
            weight.sel <- ifelse(gt.sel$blacklist==0 & gt.sel$unmasked>0.9, 1, 0)
            gc.residuals.sel <- limma::loessFit(y=lr.sel, x=gt.sel$gc,
                                                weight=weight.sel, min.weight=1e-9, span=0.35)$residuals
            if (opts$gc.adjust.offset) {
                lr.offset.sel <- lm(gc.residuals.sel~lr.sel, weights=weight.sel)$coefficients[1]
                lr[gt$target==sel] <- gc.residuals.sel - lr.offset.sel
            } else {
                lr[gt$target==sel] <- gc.residuals.sel
            }
        }
    } 
    mx.lr <- cbind(lr.raw, lr)
    return(mx.lr)
}
    
.fixOutliers <- function(gt, opts) {
    ## remove gross outliers
    if (any(gt$target)) {
        gt$lr[ gt$target] <- .smooth.outliers.gr(gt[ gt$target], "lr")
    }
    if (any(!gt$target)) {
        gt$lr[!gt$target] <- .smooth.outliers.gr(gt[!gt$target], "lr")
    }
    return(gt)
}
    
.addCoverage <- function(gt, t.bam, n.bam, opts) {
    ## normalize, smooth, and gc-correct
    t.cov <- .runMosdepthTile(t.bam, gt)
    n.cov <- .runMosdepthTile(n.bam, gt)
    ## normalize by sequencing depth
    mx.cov <- .normCoverage(t.cov, n.cov, gt, opts)
    mx.lr <- .covarageRatio(mx.cov, gt, opts)
    ## 
    mcols(gt) <- cbind(mcols(gt), cbind(mx.cov, mx.lr))
    gt <- .fixOutliers(gt)
    return(gt)
}

.addBaf <- function(gt, snp, opts) {
    hits <- findOverlaps(snp, gt, maxgap = opts$shoulder-1)
    hits <- hits[gt[subjectHits(hits)]$target] # prefer assignment to target
    hits <- hits[!duplicated(queryHits(hits))] # if snp close to two targets pick one
    bad=ifelse(snp[queryHits(hits)]$t.AF<0.5,
               round(   snp[queryHits(hits)]$t.AF  * snp[queryHits(hits)]$t.DP),
               round((1-snp[queryHits(hits)]$t.AF) * snp[queryHits(hits)]$t.DP))
    tmp <- data.table(
        idx=subjectHits(hits),
        bad=bad,
        depth=snp[queryHits(hits)]$t.DP
    )
    setkey(tmp, idx)
    tmp <- tmp[J(1:length(gt))]
    tmp <- tmp[,.(baf=sum(bad)/sum(depth)),by=idx]
    gt$baf <- tmp$baf
    gt$baf <- .smooth.outliers.gr(gt, "baf")
    return(gt)
}

.annotateTiles <- function(genome.tile, seqi, skip.chr) {
    
    bl1.fn <- system.file("extdata/sv-blacklist-10x-hg38-ucsc.bed.gz", package="gscars")
    bl2.fn <- system.file("extdata/hg38.blacklist.bed.gz", package="gscars")
    cyto.fn <- system.file("extdata/hg38.cytoband.bed.gz", package="gscars")
    strict.fn <- system.file("extdata/1000G-strict-unmask-hg38.bed.gz", package="gscars")

    ## blacklist
    bl.1 <- .robust.import(bl1.fn, seqi)
    bl.2 <- .robust.import(bl2.fn, seqi)
    bl <- reduce(c(granges(bl.1), granges(bl.2)))
    tmp <- findOverlaps(genome.tile, bl)
    tmp <- data.table(
        tile=queryHits(tmp),
        blacklist=width(pintersect(genome.tile[queryHits(tmp)], bl[subjectHits(tmp)]))
    )
    setkey(tmp, tile)
    tmp <- tmp[J(seq_along(genome.tile))]
    tmp[is.na(blacklist), blacklist:=0]
    tmp <- tmp[,.(blacklist=sum(blacklist)),by=tile]
    genome.tile$blacklist <- tmp$blacklist / width(genome.tile)

    ## masking
    strict <- .robust.import(strict.fn, seqi)
    tmp <- findOverlaps(genome.tile, strict)
    tmp <- data.table(
        tile=queryHits(tmp),
        unmasked=width(pintersect(genome.tile[queryHits(tmp)], strict[subjectHits(tmp)]))
    )
    setkey(tmp, tile)
    tmp <- tmp[J(seq_along(genome.tile))]
    tmp[is.na(unmasked), unmasked:=0]
    tmp <- tmp[,.(unmasked=sum(unmasked)),by=tile]
    genome.tile$unmasked <- tmp$unmasked / width(genome.tile)
    
    ## GC content
    tmp <- getSeq(BSgenome.Hsapiens.UCSC.hg38, genome.tile)
    genome.tile$gc <- letterFrequency(tmp, "GC", as.prob=TRUE)[,1]

    ## cytobands
    cyto <- .robust.import(cyto.fn, seqi, skip.chr=skip.chr)
    tmp <- findOverlaps(genome.tile, cyto, select="first")
    genome.tile$cytoband <- cyto[tmp]$name
    genome.tile$arm <- paste0(seqnames(genome.tile), str_sub(genome.tile$cytoband, 1, 1))

    ## 
    genome.tile <- sort(genome.tile)
    return(genome.tile)
}

.getGenomeTiles <- function(tgt.fn, opts) {
    ## tile genome
    seqi <- keepStandardChromosomes(Seqinfo(genome="hg38"))
    genome.tile <- tileGenome(dropSeqlevels(seqi, opts$skip.chr), tilewidth=opts$tilewidth, cut.last.tile.in.chrom=TRUE)
    genome.tile$target <- TRUE
    genome.tile <- .annotateTiles(genome.tile, seqi, opts$skip.chr)
    return(genome.tile)
}

.getTargetTiles <- function(tgt.fn, opts) {
    ## tile targets
    seqi <- keepStandardChromosomes(Seqinfo(genome="hg38"))
    tgt <- .robust.import(tgt.fn, seqi=seqi, skip.chr = opts$skip.chr)
    tgt$name <- NULL
    tgt$score <- NULL
    tgt$target <- TRUE
    gap <- dropSeqlevels(gaps(tgt), opts$skip.chr, pruning.mode="coarse")
    gap <- gap[strand(gap)=="*"]
    gap <- gap[width(gap)>2*opts$shoulder+opts$min.gap]
    gap <- gap-opts$shoulder
    gap$target <- FALSE
    target.tile <- sort(c(tgt, gap))
    target.tile <- .annotateTiles(target.tile, seqi, opts$skip.chr)
    return(target.tile)
}
