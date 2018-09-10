OPTS <- list(
    #### basic
    cores=4,
    sample=c(1,2),
    target="onco1500-v4",
    target.fn=system.file("extdata/onco1500-v4-targets-hg38-ucsc.bed.gz", package="cnvex"),
    gene.fn=system.file("extdata/gene/motr.v2-prot.gff3", package="cnvex"),
    caller="vardict",
    chr.names="ucsc",
    skip.chr="chrM",

    ## tiling
    tile.width=NA_integer_,
    tile.min.gap=50000,
    tile.shoulder=300,

    ## log-ratio
    lr.smooth="outlier",
    lr.smooth.window=13,
    lr.loc.tileflank=NA_integer_,

    #### segmentation
    seg.strategy="joint",
    seg.method="RBS",
    seg.sd.prune=TRUE,
    seg.sd.lr.penalty=1,
    seg.sd.baf.penalty=1,
    seg.local.sd=FALSE,
    seg.min.tile=3,

    #### GC-content
    gc.adjust.trend=TRUE,
    gc.adjust.offset=TRUE,
    gc.adjust.span.on=0.15,
    gc.adjust.span.off=0.30,
    gc.adjust.on=c(0.3, 0.7),
    gc.adjust.off=c(0.0, 1.0),

    #### optimization
    grid.n=64,
    p.lo=0.05,
    p.hi=0.95,
    D.lo=1,
    D.hi=6,
    max.C=9,
    max.sC=20,
    max.len.per.probe=1e6,
    res=0.1
)
