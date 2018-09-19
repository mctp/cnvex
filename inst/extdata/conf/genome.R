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
    baf.max.eff.dp=300,
    baf.het.range=c(0.3, 0.7),
    baf.min.het.dp=6,
    baf.min.target.dp=24,
    baf.min.genome.dp=12,

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
    opt.grid.n=64,
    opt.p.lo=0.05,
    opt.p.hi=0.95,
    opt.D.lo=1,
    opt.D.hi=6,
    opt.max.C=9,
    opt.max.sC=20,
    opt.max.len.per.probe=1e6,
    opt.res=0.1

)
