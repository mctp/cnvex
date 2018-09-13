OPTS <- list(
    #### basic
    cores=4,
    sample=c(1,2),
    target="genome",
    target.fn=NA_character_,
    gene.fn=system.file("extdata/gene/motr.v2-prot.gff3", package="cnvex"),
    caller="vardict",
    chr.names="ucsc",
    skip.chr="chrM",

    ## tiling
    tile.width=10000,
    tile.min.gap=NA_integer_,
    tile.shoulder=0,

    ## log-ratio
    lr.smooth="outlier",
    lr.smooth.window=21,
    lr.loc.tileflank=5,

    #### BAF
    baf.min.t.dp=12,
    baf.het.range=c(0.3, 0.7),
    baf.qual.min=0.25,

    #### segmentation
    seg.strategy="joint",
    seg.method="CBS",
    seg.sd.prune=TRUE,
    seg.sd.lr.penalty=1,
    seg.sd.baf.penalty=1,
    seg.len.prune=TRUE,
    seg.len.min=2,
    seg.cbs.baf=list(alpha=0.01, trim=0.025, min.width=2),
    seg.cbs.lr=list(alpha=0.01, trim=0.025, min.width=2),
    seg.rbs.selection="Lebarbier",
    
    #### GC-content
    gc.adjust.trend=TRUE,
    gc.adjust.offset=TRUE,
    gc.adjust.span.on=0.5,
    gc.adjust.span.off=NA_real_,
    gc.adjust.on=c(0.3, 0.7),
    gc.adjust.off=c(NA_real_, NA_real_),

    #### optimization
    opt.local.sd=FALSE,
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
