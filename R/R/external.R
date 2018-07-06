.runMosdepth <- function(bam, by) {
    if (is.character(by) && file.exists(by)) {
        out.col <- "V5"
    } else if (!is.na(suppressWarnings(as.integer(by))) ) {
        out.col <- "V4"
    } else {
        stop("by should be a file or window-size")
    }
    prefix <- tempfile("mosdepth_")
    ret <- system2("tools/bin/mosdepth", sprintf("-F 772 -n -t4 -b %s %s %s", by, prefix, bam))
    out.fn <- list.files(dirname(prefix), paste0(basename(prefix), ".regions.bed.gz$"), full.names = TRUE)
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

.runMosdepthTile <- function(bam, gt) {
    bed <- tempfile("mosdepth_", fileext=".bed")
    export(granges(gt), bed)
    out <- .runMosdepth(bam, bed)
    unlink(bed)
    return(out)
}
