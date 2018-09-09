.annotateTiles <- function(genome.tile, seqi, skip.chr) {
    
    bl1.fn <- system.file("extdata/sv-blacklist-10x-hg38-ucsc.bed.gz", package="cnvex")
    bl2.fn <- system.file("extdata/hg38.blacklist.bed.gz", package="cnvex")
    cyto.fn <- system.file("extdata/hg38.cytoband.bed.gz", package="cnvex")
    strict.fn <- system.file("extdata/1000G-strict-unmask-hg38.bed.gz", package="cnvex")

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

    ## order arms
    genome.tile <- sort(genome.tile)
    genome.tile$arm <- paste0(seqnames(genome.tile), str_sub(genome.tile$cytoband, 1, 1))
    genome.tile$arm <- factor(genome.tile$arm, unique(genome.tile$arm), ordered=TRUE)
    return(genome.tile)
}

.getGenomeTiles <- function(tgt.fn, opts) {
    ## tile genome
    seqi <- keepStandardChromosomes(Seqinfo(genome="hg38"))
    genome.tile <- tileGenome(dropSeqlevels(seqi, opts$skip.chr), tilewidth=opts$tile.width, cut.last.tile.in.chrom=TRUE)
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
    gap <- gap[width(gap)>2*opts$tile.shoulder+opts$tile.min.gap]
    gap <- gap-opts$tile.shoulder
    gap$target <- FALSE
    target.tile <- sort(c(tgt, gap))
    target.tile <- .annotateTiles(target.tile, seqi, opts$skip.chr)
    return(target.tile)
}

getTile <- function(opts) {
    if (opts$target=="genome") {
        target.fun <- .getGenomeTiles
    } else {
        target.fun <- .getTargetTiles
    }
    tile <- target.fun(opts$target.fn, opts)
    return(tile)
}
