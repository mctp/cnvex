.detect.sex <- function(var, tile) {
    ## placeholder implementation of sex detection
    chrx.sel <- tile$target & (seqnames(tile) %in% c("chrX"))
    chrx.cov <- median(tile[chrx.sel]$n.cov)
    chry.sel <- tile$target & (seqnames(tile) %in% c("chrY"))
    chry.cov <- median(tile[chry.sel]$n.cov)
    sex <- if (chrx.cov/chry.cov > 1.5) "female" else "male"
    return(sex)
}

.getCopy <- function(gt, sex) {
    copy <- rep(2L, length(gt))
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
    return(copy)
}

.addCopy <- function(gr, sex) {
    gr$nC <- .getCopy(gr, sex)
    return(gr)
}

CNV <- function(opts, tn.vcf, t.bam, n.bam) {

    if (opts$target=="onco1500-v3") {
        target.fun <- .getTargetTiles
        target.fn <- system.file(sprintf("extdata/onco1500-v3-targets-hg38-%s.bed", opts$chr.names), package="gscars")
        snp.filter.fun <- filterTargetGermlineHets
    } else if (opts$target=="agilent-v4") {
        target.fun <- .getTargetTiles
        target.fn <- system.file(sprintf("extdata/agilent-v4-targets-hg38-%s.bed", opts$chr.names), package="gscars")
        snp.filter.fun <- filterTargetGermlineHets
    } else if(opts$target=="nextera-v1.2") {
        target.fun <- .getTargetTiles
        target.fn <- system.file(sprintf("extdata/nextera-v1.2-targets-hg38-%s.bed", opts$chr.names), package="gscars")
        snp.filter.fun <- filterTargetGermlineHets
    } else if (opts$target=="genome") {
        target.fun <- .getGenomeTiles
        target.fn <- NA_character_
        snp.filter.fun <- filterGenomeGermlineHets
    } else {
        stop("Invalid target.")
    }
    
    ## import all data
    var <- importVcf(tn.vcf, opts)
    tile <- target.fun(target.fn, opts)
    tile <- .addCoverage(tile, t.bam, n.bam, opts)
    
    ## detect sex
    if (is.null(opts$sex)) {
        sex <- .detect.sex(var, tile)
    } else {
        sex <- opts$sex
    }
    tile <- .addCopy(tile, sex)
    var <- .addCopy(var, sex)

    ## add BAF
    snp <- snp.filter.fun(var, tile, opts)
    tile <- .addBaf(tile, snp, opts)
    
    cnv <- list(snp=snp, tile=tile, sex=sex)
    return(cnv)
}
