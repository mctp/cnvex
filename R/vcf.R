.importVcf <- function(vcf, param, opts) {
    seqi <- keepStandardChromosomes(Seqinfo(genome="hg38"))
    seqi <- dropSeqlevels(seqi, opts$skip.chr)
    var <- readVcf(vcf, param=param)
    var <- var[,opts$sample]
    shared.levels <- intersect(seqlevels(var), seqlevels(seqi))
    var <- keepSeqlevels(var, shared.levels, pruning.mode="coarse")
    seqlevels(var) <- seqlevels(seqi)
    ## masked regions 1000G
    pilot.fn <- system.file("extdata/1000G-pilot-unmask-hg38.bed.gz", package="cnvex")
    pilot <- .robust.import(pilot.fn, seqi)
    strict.fn <- system.file("extdata/1000G-strict-unmask-hg38.bed.gz", package="cnvex")
    strict <- .robust.import(strict.fn, seqi)
    tmp <- DataFrame(Number=c("1","1"),
                     Type=c("Integer", "Integer"),
                     Description=c("loose mask", "strict mask"),
                     row.names=c("mask.loose", "mask.strict"))
    info(header(var)) <- rbind(info(header(var)), tmp)
    info(var)$mask.loose <- as.integer(!(var %over% pilot))
    info(var)$mask.strict <- as.integer(!(var %over% strict))
    qual(var) <- ifelse(is.na(qual(var)), 0, qual(var))
    ## ignore multi-allelics
    var <- var[
        elementNROWS(fixed(var)$ALT) == 1
    ]
    return(var)
}

.importVcfVarDict <- function(vcf, opts) {
    param <-ScanVcfParam(
        info=c("SOMATIC", "TYPE"),
        geno=c("GT", "AF", "DP", "QUAL")
    )
    if (!is.null(opts$which)) {
        vcfWhich(param) <- which
    }
    .importVcf(vcf, param, opts)
}

.importVarDict <- function(vcf, opts) {
    var <- .importVcfVarDict(vcf, opts)
    var.gr <- rowRanges(var)
    var.gr$mask.loose <- info(var)$mask.loose
    var.gr$mask.strict <- info(var)$mask.strict
    var.gr$SOMATIC <- info(var)$SOMATIC
    var.gr$TYPE <- info(var)$TYPE
    var.gr$t.GT <- geno(var)$GT[,1]
    var.gr$n.GT <- geno(var)$GT[,2]
    var.gr$t.AF <- geno(var)$AF[,1]
    var.gr$n.AF <- geno(var)$AF[,2]
    var.gr$t.DP <- geno(var)$DP[,1]
    var.gr$n.DP <- geno(var)$DP[,2]
    return(var.gr)
}

.importVcfHC <- function(vcf, opts) {
    param <-ScanVcfParam(
        info=c("AF"),
        geno=c("GT", "AD", "DP", "GQ")
    )
    if (!is.null(opts$which)) {
        vcfWhich(param) <- which
    }
    .importVcf(vcf, param, opts)
}

.importHC <- function(vcf, opts) {
    var <- .importVcfHC(vcf, opts)
    var.gr <- rowRanges(var)
    var.gr$mask.loose <- info(var)$mask.loose
    var.gr$mask.strict <- info(var)$mask.strict
    var.gr$SOMATIC <- FALSE
    is.snv <- width(fixed(var)$REF)==1 & width(unlist(fixed(var)$ALT))==1
    var.gr$TYPE <- ifelse(is.snv, "SNV", "Other")
    t.AD <- sapply(geno(var)$AD[,1], "[", 2)
    n.AD <- sapply(geno(var)$AD[,2], "[", 2)
    var.gr$t.GT <- geno(var)$GT[,1]
    var.gr$n.GT <- geno(var)$GT[,2]
    var.gr$t.AF <- t.AD / geno(var)$DP[,1]
    var.gr$n.AF <- n.AD / geno(var)$DP[,2]
    var.gr$t.DP <- geno(var)$DP[,1]
    var.gr$n.DP <- geno(var)$DP[,2]
    return(var.gr)
}

importVcf <- function(vcf, opts) {
    ## import
    if (opts$caller=="vardict") {
        var <- unlist(GRangesList(unname(lapply(names(vcf), function(name) {
            fn <- vcf[[name]]
            v <- .importVarDict(fn, opts)
            mcols(v)$SOURCE <- name
            return(v)
        }))))
    } else if (opts$caller=="hc") {
        var <- unlist(GRangesList(unname(lapply(names(vcf), function(name) {
            fn <- vcf[[name]]
            v <- .importHC(fn, opts)
            mcols(v)$SOURCE <- name
            return(v)
        }))))
    } else {
        stop("Variant caller not supported.")
    }
    return(var)
}
