OPTS <- list(
    ## basic
    cores=4,
    sample=c(1,2),
    target="nextera-v1.2",
    target.fn=system.file("extdata/nextera-v1.2-targets-hg38-ucsc.bed.gz", package="cnvex"),
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
    lr.smooth.window=21,
    lr.loc.tileflank=5,

    ## BAF
    baf.max.eff.dp=300,
    baf.het.range=c(0.3, 0.7),
    baf.min.het.dp=6,
    baf.min.target.dp=24,
    baf.min.genome.dp=12,

    ## segmentation
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

    ## GC-content
    gc.adjust.trend=TRUE,
    gc.adjust.offset=TRUE,
    gc.adjust.span.on=0.15,
    gc.adjust.span.off=0.30,
    gc.adjust.on=c(0.3, 0.7),
    gc.adjust.off=c(0.0, 1.0),

    ## pool
    pool.method="pca",
    pool.lo.cov=0.25,
    pool.hi.cov=4,
    pool.hi.zero=0.20,
    pool.sd.out=3,
    pool.n.comp=20,

    ## optimization
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
