.runMosdepth <- function(bam, by, cores=4L) {
    if (is.character(by) && file.exists(by)) {
        out.col <- "V5"
    } else if (!is.na(suppressWarnings(as.integer(by))) ) {
        out.col <- "V4"
    } else {
        stop("by should be a file or window-size")
    }
    prefix <- tempfile("mosdepth_")
    ret <- system2("mosdepth", sprintf("-F 772 -n -t%s -b %s %s %s", cores, by, prefix, bam))
    out.fn <- list.files(dirname(prefix), paste0(basename(prefix), ".regions.bed.gz$"),
                         full.names = TRUE)
    if (!ret & file.exists(out.fn)) {
        out <- fread(paste("zcat", out.fn))[[out.col]]
    } else {
        stop(sprintf("mosdepth run failed: %s", ret))
    }
    out.fns <- list.files(dirname(prefix), paste0(basename(prefix)), full.names = TRUE)
    rets <- sapply(out.fns, unlink)
    if (any(rets)) {
        stop("could not remove temp. files")
    }
    return(out)
}

.runMosdepthTile <- function(bam, gt, cores) {
    bed <- tempfile("mosdepth_", fileext=".bed")
    export(granges(gt), bed)
    out <- .runMosdepth(bam, bed, cores)
    unlink(bed)
    return(out)
}

.runCBS <- function(y, ...) {
    n <- length(y)
    chrom <- rep(1, n)
    maploc <- 1:n
    genomdat <- y
    cna <- DNAcopy::CNA(genomdat, chrom, maploc)
    capture.output(
        res <- DNAcopy::segment(cna, ...)
    )
    bkp <- res$output$loc.end[-length(res$output$loc.end)]
    if (length(bkp)==0) {
        bkp <- integer()
    }
    return(bkp)
}
