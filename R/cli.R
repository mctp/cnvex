#' @export
importCNVEX <- function(opts, tn.vcf, t.bam, n.bam) {
    tile <- getTile(opts)
    var <- importVcf(tn.vcf, tile, opts)
    tile <- rawLogRatio(t.bam, n.bam, tile, opts)
    cnv <- list(var=var, tile=tile)
    cnv <- getGene(cnv, opts)
    return(cnv)
}

#' @export
adjustCNVEX <- function(opts, cnv) {
    cnv <- gcLogRatio(cnv, opts)
    cnv <- smoothLogRatio(cnv, opts)
    return(cnv)
}

#' @export
basic <- function() {

    option_list = list(
        optparse::make_option(c("-c", "--config"), type="character",
                            default="genome",
                            help="Preset file"),
        optparse::make_option(c("-p", "--cores"), type="integer",
                            default=4L,
                            help="Number of cores"),
        optparse::make_option(c("-t", "--tumor"), type="character",
                              default=NULL,
                              help="preset"),
        optparse::make_option(c("-n", "--normal"), type="character",
                              default=NULL,
                              help="preset"),
        optparse::make_option(c("-v", "--vcf"), type="character",
                              default=NULL,
                              help="Somatic VCF file"),
        optparse::make_option(c("-o", "--out"), type="character",
                              default="cnvex.rds",
                              help="Output file")
        
    )
    parser = optparse::OptionParser(
      "Rscript -e 'cnvex::basic()' [options]",
      description=c("Import data to create CNVEX object.\n"),
      epilogue=c(
          "Written by Marcin Cieslik (mcieslik@med.umich.edu)",
          "Michigan Center for Translational Pathology (c) 2018\n"),
      option_list=option_list
      )

    args <- optparse::parse_args(parser, positional_arguments=FALSE)
    opts <- getOpts(args$config, list(cores=args$cores))
    cnvex <- importCNVEX(opts, args$vcf, args$tumor, args$normal)
    saveRDS(cnvex, args$out)
}

#' @export
tocsv <- function() {
    
    option_list = list(
        
        optparse::make_option(c("-c", "--config"), type="character",
                            default="genome",
                            help="Preset file"),
        optparse::make_option(c("-p", "--cores"), type="integer",
                            default=4L,
                            help="Number of cores"),
        optparse::make_option(c("-t", "--type"), type="character",
                              default="genelr",
                              help="Output type"),
        optparse::make_option(c("-i", "--inp"), type="character",
                              default=NULL,
                              help="Input CNVEX file"),
        optparse::make_option(c("-o", "--out"), type="character",
                              default=NULL,
                              help="Output CSV file")
        
    )
    parser = optparse::OptionParser(
      "Rscript -e 'cnvex::tocsv()' [options]",
      description=c("Export CNVEX object as a CSV file.\n"),
      epilogue=c(
          "Written by Marcin Cieslik (mcieslik@med.umich.edu)",
          "Michigan Center for Translational Pathology (c) 2018\n"),
      option_list=option_list
      )

    args <- optparse::parse_args(parser, positional_arguments=FALSE)
    opts <- getOpts(args$config, opts=list(cores=args$cores))
    cnv <- readRDS(args$inp)
    if (args$type == "genelr") {
        csv <- .dtRound(.geneLogRatio(cnv, opts))
    } else {
        optparse::print_help(parser)
        write("Output type not supported.\n", stderr())
        quit("no", 1)
    }
    if (is.null(args$out)) {
        args$out <- str_replace(args$inp, ".rds$", sprintf("-%s.csv", args$type))
    }
    fwrite(csv, args$out)
}

#' @export
adjust <- function() {

    option_list = list(
        
        optparse::make_option(c("-c", "--config"), type="character",
                            default="genome",
                            help="Preset file"),
        optparse::make_option(c("-p", "--cores"), type="integer",
                            default=4L,
                            help="Number of cores"),
        optparse::make_option(c("-i", "--inp"), type="character",
                              default=NULL,
                              help="Input CNVEX file"),
        optparse::make_option(c("-o", "--out"), type="character",
                              default=NULL,
                              help="Output CNVEX file")
        
    )
    parser = optparse::OptionParser(
      "Rscript -e 'cnvex::adjust()' [options]",
      description=c("Import data to create CNVEX object.\n"),
      epilogue=c(
          "Written by Marcin Cieslik (mcieslik@med.umich.edu)",
          "Michigan Center for Translational Pathology (c) 2018\n"),
      option_list=option_list
      )

    args <- optparse::parse_args(parser, positional_arguments=FALSE)
    opts <- getOpts(args$config, opts=list(cores=args$cores))
    cnv <- readRDS(args$inp)
    cnv <- adjustCNVEX(opts, cnv)
    if (is.null(args$out)) {
        args$out <- str_replace(args$inp, ".rds$", "-adj.rds")
    }
    saveRDS(cnv, args$out)
}

#' @export
segment <- function() {

    option_list = list(
        
        optparse::make_option(c("-c", "--config"), type="character",
                            default="genome",
                            help="Preset file"),
        optparse::make_option(c("-p", "--cores"), type="integer",
                            default=4L,
                            help="Number of cores"),
        optparse::make_option(c("-i", "--inp"), type="character",
                              default=NULL,
                              help="Input CNVEX file"),
        optparse::make_option(c("-o", "--out"), type="character",
                              default=NULL,
                              help="Output CNVEX file")
        
    )
    parser = optparse::OptionParser(
      "Rscript -e 'cnvex::segment()' [options]",
      description=c("Import data to create CNVEX object.\n"),
      epilogue=c(
          "Written by Marcin Cieslik (mcieslik@med.umich.edu)",
          "Michigan Center for Translational Pathology (c) 2018\n"),
      option_list=option_list
      )

    args <- optparse::parse_args(parser, positional_arguments=FALSE)
    opts <- getOpts(args$config, opts=list(cores=args$cores))
    cnv <- readRDS(args$inp)
    cnv <- getBaf(cnv, opts)
    cnv <- jointSegment(cnv, opts)
    cnv <- getSeg(cnv, opts)
    if (is.null(args$out)) {
        args$out <- str_replace(args$inp, ".rds$", "-seg.rds")
    }
    saveRDS(cnv, args$out)
}

#' @export
plot <- function() {
    suppressPackageStartupMessages(library(ggpubr))
    suppressPackageStartupMessages(library(gridExtra))

    option_list = list(
        
        optparse::make_option(c("-c", "--config"), type="character",
                            default="genome",
                            help="Preset file"),
        optparse::make_option(c("-p", "--cores"), type="integer",
                            default=4L,
                            help="Number of cores"),
        optparse::make_option(c("-t", "--type"), type="character",
                              default="cnv",
                              help="Plot type"),
        optparse::make_option(c("-l", "--lr"), type="character",
                              default="lr.smooth",
                              help="log2ration type [lr, lr.gc, lr.smooth]"),
        optparse::make_option(c("-i", "--inp"), type="character",
                              default=NULL,
                              help="Input CNVEX file"),
        optparse::make_option(c("-r", "--chr"), type="character",
                              default=NULL,
                              help="chromosome"),
        optparse::make_option(c("-o", "--out"), type="character",
                              default=NULL,
                              help="Output file name")
        
    )
    parser = optparse::OptionParser(
      "Rscript -e 'cnvex::plot()' [options]",
      description=c("Plot CNVEX object.\n"),
      epilogue=c(
          "Written by Marcin Cieslik (mcieslik@med.umich.edu)",
          "Michigan Center for Translational Pathology (c) 2018\n"),
      option_list=option_list
      )

    args <- optparse::parse_args(parser, positional_arguments=FALSE)
    opts <- getOpts(args$config, opts=list(cores=args$cores))
    
    cnv <- readRDS(args$inp)
    if (args$type=="gc") {
        plt <- plotGC(cnv)
        if (is.null(args$out)) {
            args$out <- str_replace(args$inp, ".rds$", "-gc.png")
        }
        suppressMessages(ggsave(args$out, plt, width=7, height=7))
    }
    if (args$type=="cnv") {
        plt <- plotCNV(cnv, sel.lr = args$lr, sel.chr = args$chr)
        if (is.null(args$out)) {
            if (is.null(args$chr)) {
                args$out <- str_replace(args$inp, ".rds$", sprintf("-cnv-%s.png", args$lr))
            } else {
                args$out <- str_replace(args$inp, ".rds$", sprintf("-cnv-%s-%s.png", args$lr, args$chr))
            }
        }
        suppressMessages(ggsave(args$out, plt, width=15, height=5))
    }
    if (args$type=="seg") {
        plt <- plotSeg(cnv, sel.chr = args$chr)
        if (is.null(args$out)) {
            if (is.null(args$chr)) {
                args$out <- str_replace(args$inp, ".rds$", "-seg.png")
            } else {
                args$out <- str_replace(args$inp, ".rds$", sprintf("-seg-%s.png", args$chr))
            }
        }
        suppressMessages(ggsave(args$out, plt, width=15, height=5))
    }
}
