getOpts <- function(target, opts=list()) {
    if (target=="onco1500-v3") {
        opts.default <- list(
            target=target,
            segment="joint",
            caller="vardict",
            chr.names="ucsc",
            skip.chr="chrM",
            min.gap=50000,
            shoulder=300,
            tilewidth=NA_integer_,
            min.seg=2,
            K=30,
            sd.penalty=1,
            local.sd=FALSE,
            grid.n=64,
            p.lo=0.05,
            p.hi=0.95,
            D.lo=1,
            D.hi=6,
            max.C=9,
            max.sC=20,
            max.len.per.probe=1e6,
            res=0.1,
            gc.adjust.trend=FALSE,
            gc.adjust.offset=FALSE
        )
    } else if (target=="agilent-v4") {
        opts.default <- list(
            target=target,
            segment="joint",
            caller="vardict",
            chr.names="ucsc",
            skip.chr="chrM",
            min.gap=50000,
            shoulder=300,
            tilewidth=NA_integer_,
            min.seg=3,
            K=30,
            sd.penalty=1,
            local.sd=FALSE,
            grid.n=64,
            p.lo=0.05,
            p.hi=0.95,
            D.lo=1,
            D.hi=6,
            max.C=9,
            max.sC=20,
            max.len.per.probe=1e6,
            res=0.1,
            gc.adjust.trend=FALSE,
            gc.adjust.offset=FALSE
        )
    } else if (target=="nextera-v1.2") {
        opts.default <- list(
            target=target,
            segment="joint",
            caller="vardict",
            chr.names="ucsc",
            skip.chr="chrM",
            min.gap=50000,
            shoulder=300,
            tilewidth=NA_integer_,
            min.seg=3,
            K=30,
            sd.penalty=1,
            local.sd=FALSE,
            grid.n=64,
            p.lo=0.05,
            p.hi=0.95,
            D.lo=1,
            D.hi=6,
            max.C=9,
            max.sC=20,
            max.len.per.probe=1e6,
            res=0.1,
            gc.adjust.trend=TRUE,
            gc.adjust.offset=TRUE
        )        
    } else if (target=="genome") {
        opts.default <- list(
            target=target,
            segment="joint",
            caller="vardict",
            chr.names="ucsc",
            skip.chr="chrM",
            min.gap=NA_integer_,
            shoulder=0,
            tilewidth=10000,
            min.seg=20,
            K=30,
            sd.penalty=1,
            local.sd=FALSE,
            grid.n=64,
            p.lo=0.05,
            p.hi=0.95,
            D.lo=1,
            D.hi=6,
            max.C=9,
            max.sC=20,
            max.len.per.probe=1e6,
            res=0.1,
            gc.adjust.trend=FALSE,
            gc.adjust.offset=FALSE
        )
    } else {
        stop("Invalid target.")
    }
    opts <- merge.list(opts, opts.default)
    return(opts)
}
