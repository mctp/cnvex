#' @export
addCoverage <- function(t.bam, n.bam, cnv, opts) {
    ## normalize by sequencing depth
    t.cov <- .runMosdepthTile(t.bam, cnv$tile, opts$cores)
    t.cov <- .normCoverage(t.cov, cnv$tile, opts)
    if (!is.null(n.bam)) {
        n.cov <- .runMosdepthTile(n.bam, cnv$tile, opts$cores)
        n.cov <- .normCoverage(n.cov, cnv$tile, opts)
    } else {
        n.cov <- NA_real_
    }
    mcols(cnv$tile) <- cbind(mcols(cnv$tile), cbind(t.cov, n.cov))
    return(cnv)
}

#' @export
addFakeLogRatio <- function(cnv, pool, opts) {
    p.cov <- .poolFakeCoverage(cnv, pool, opts)
    cnv$tile$lr.raw <- .rawLogRatio(cnv$tile$t.cov, p.cov, opts)
    return(cnv)
}

#' @export
addPcaLogRatio <- function(cnv, pool, opts) {
    cnv$tile$lr.raw <- .poolPcaDenoise(cnv, pool, opts)
    return(cnv)
}

#' @export
addIcaLogRatio <- function(cnv, pool, opts) {
    cnv$tile$lr.raw <- .poolIcaDenoise(cnv, pool, opts)
    return(cnv)
}

#' @export
addPairLogRatio <- function(cnv, opts) {
    cnv$tile$lr.raw <- .rawLogRatio(cnv$tile$t.cov, cnv$tile$n.cov, opts)
    return(cnv)
}

#' @export
addGcLogRatio <- function(cnv, opts) {
    cnv$tile$lr.gc <- .gcLogRatio(cnv$tile$lr.raw, cnv$tile, opts)
    return(cnv)
}

#' @export
addSmoothLogRatio <- function(cnv, opts) {
    cnv$tile$lr.smooth <- cnv$tile$lr.gc
    if (opts$lr.smooth=="hybrid") {
        ## hybrid smooth only on targeted
        cnv$tile$lr.smooth[cnv$tile$target] <-
            .smoothLogRatio(cnv$tile[cnv$tile$target]$lr.smooth, cnv$tile[cnv$tile$target], opts)
    } else if (opts$lr.smooth=="outlier") {
        cnv$tile$lr.smooth <- .smoothOutliers(cnv$tile$lr.smooth, cnv$tile, opts)
    }
    return(cnv)
}

#' @export
addBaf <- function(cnv, opts) {
    cnv$tile <- .getBaf(cnv$tile, cnv$var, opts)
    return(cnv)
}

#' @export
addJointSegment <- function(cnv, opts) {
    ## estimated standard-deviation
    sd.lr <- estimateSd(cnv$tile$lr.smooth) * opts$seg.sd.lr.penalty
    sd.baf <- estimateSd(cnv$tile$baf) * opts$seg.sd.baf.penalty
    ## create segmentation
    cnv$tile <- .addJointSeg(cnv$tile, sd.lr, sd.baf, opts$seg.method,
                             opts$tile.width, opts$seg.len.min, opts$seg.cbs.lr,
                             opts$seg.cbs.baf, opts$seg.rbs.selection, opts$seg.sd.prune, opts$seg.len.prune)
    cnv$seg <- unname(unlist(range(split(cnv$tile, cnv$tile$seg))))
    ## seg -> var
    cnv$var$seg <- findOverlaps(cnv$var, cnv$seg, select="first", maxgap = opts$tile.shoulder-1)
    ## seg -> gene
    hit.prom <- as.data.table(findOverlaps(promoters(cnv$gene, 0, 1), cnv$seg))
    hit.gene <- as.data.table(findOverlaps(cnv$gene, cnv$seg))
    hit.gene <- hit.gene[!(queryHits %in% hit.prom$queryHits)]
    hit.gene <- hit.gene[,.SD[1],by=queryHits]
    hit <- rbind(hit.prom, hit.gene)
    cnv$gene$seg <- NA_real_
    cnv$gene$seg[hit$queryHits] <- hit$subjectHits
    return(cnv)
}

#' @export
addCopy <- function(cnv, opts) {
    if (is.null(opts$sex)) {
        sex <- .detect.sex(cnv$var, cnv$tile)
    } else {
        sex <- opts$sex
    }
    copy <- rep(2L, length(cnv$tile))
    if (sex=="male") {
        x.copy <- 1L
        y.copy <- 1L
    } else if (sex=="female"){
        x.copy <- 2L
        y.copy <- 0L
    } else {
        stop("wrong sex.")
    }
    copy[as.logical(seqnames(gt) %in% c("chrX"))] <- x.copy
    copy[as.logical(seqnames(gt) %in% c("chrY"))] <- y.copy
    cnv$tile$nC <- copy
    return(cnv)
}

#' @export
addGene <- function(cnv, opts) {
    gtf <- import(opts$gene.fn)
    gtf.gene <- gtf[gtf$type=="gene",c("source", "type", "gene_id", "gene_name")]
    cnv$gene <- gtf.gene
    return(cnv)
}

#' @export
importCNVEX <- function(vcf, t.bam, n.bam, opts) {
    var <- importVcf(vcf, opts)
    tile <- getTile(opts)
    cnv <- list(var=var, tile=tile)
    cnv <- addGene(cnv, opts)
    cnv <- addCoverage(t.bam, n.bam, cnv, opts)
    return(cnv)
}

#' @export
addLogRatio <- function(cnv, pool, opts) {
    if (is.null(pool)) {
        cnv <- addPairLogRatio(cnv, opts)
    } else {
        if (pool$method=="fake") {
            cnv <- addFakeLogRatio(cnv, pool, opts)
        }
        else if (pool$method=="pca") {
            cnv <- addPcaLogRatio(cnv, pool, opts)
        }
        else if (pool$method=="ica") {
            cnv <- addIcaLogRatio(cnv, pool, opts)
        }
    }
    cnv <- addGcLogRatio(cnv, opts)
    cnv <- addSmoothLogRatio(cnv, opts)
    return(cnv)
}

#' @export
addSegment <- function(cnv, opts) {
    cnv <- addBaf(cnv, opts)
    cnv <- addJointSegment(cnv, opts)
    return(cnv)
}

#' @export
resolveConfig <- function(config) {
    if (!file.exists(config)) {
        config <- system.file(sprintf("extdata/conf/%s.R", config), package="cnvex")
    }
    return(config)
}

#' @export
getOpts <- function(conf.fn, opts=list()) {
    ENV = new.env(parent = .BaseNamespaceEnv)
    source(conf.fn, local=ENV)
    opts <- merge.list(opts, ENV$OPTS)
    return(opts)
}

#' @export
createPool <- function(cnv.fns, opts) {
    pool.data <- .importPoolData(cnv.fns, opts)
    pool <- .createPool(pool.data, opts)
    return(pool)
}
