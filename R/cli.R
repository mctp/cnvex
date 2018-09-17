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
        optparse::make_option(c("-g", "--gvcf"), type="character",
                              default=NULL,
                              help="Genome Somatic VCF file"),
        optparse::make_option(c("-e", "--evcf"), type="character",
                              default=NULL,
                              help="Exome Somatic VCF file"),
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

    config <- resolveConfig(args$config)
    if (config=="") {
        optparse::print_help(parser)
        write("Config file not found.\n", stderr())
        quit("no", 1)
    }
    opts <- getOpts(config, opts=list(cores=args$cores))

    if (
        is.null(args$tumor) ||
        is.null(args$normal) ||
        is.null(args$vcf) ||
        !file.exists(args$tumor) ||
        !file.exists(args$normal) ||
        !file.exists(args$vcf)
    ) {
        optparse::print_help(parser)
        write("Input file(s) not found.\n", stderr())
        quit("no", 1)
    }
    cnvex <- importCNVEX(c(args$gvcf, args$evcf), args$tumor, args$normal, opts)
    saveRDS(cnvex, args$out)
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
                              help="Output CNVEX file"),
        optparse::make_option(c("-s", "--suffix"), type="character",
                              default=NULL,
                              help="Optional file name suffix")
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

    config <- resolveConfig(args$config)
    if (config=="") {
        optparse::print_help(parser)
        write("Config file not found.\n", stderr())
        quit("no", 1)
    }
    opts <- getOpts(config, opts=list(cores=args$cores))
    if(is.null(args$inp) || !file.exists(args$inp)) {
        optparse::print_help(parser)
        write("Input file not found.\n", stderr())
        quit("no", 1)
    }
    cnv <- readRDS(args$inp)

    cnv <- addCorrections(cnv, opts)
    if (is.null(args$out)) {
        if (is.null(args$suffix)) {
            args$out <- str_replace(args$inp, ".rds$", "-adj.rds")
        } else {
            args$out <- str_replace(args$inp, ".rds$", sprintf("-adj-%s.rds", args$suffix))
        }
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
                              help="Output CNVEX file"),
        optparse::make_option(c("-s", "--suffix"), type="character",
                              default=NULL,
                              help="Optional file name suffix")
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

    config <- resolveConfig(args$config)
    if (config=="") {
        optparse::print_help(parser)
        write("Config file not found.\n", stderr())
        quit("no", 1)
    }
    opts <- getOpts(config, opts=list(cores=args$cores))
    if(is.null(args$inp) || !file.exists(args$inp)) {
        optparse::print_help(parser)
        write("Input file not found.\n", stderr())
        quit("no", 1)
    }
    cnv <- readRDS(args$inp)
    
    cnv <- addBaf(cnv, opts)
    cnv <- addJointSegment(cnv, opts)
    if (is.null(args$out)) {
        if (!is.null(args$suffix)) {
            args$out <- str_replace(args$inp, ".rds$", sprintf("-seg-%s.rds", args$suffix))
        } else {
            args$out <- str_replace(args$inp, ".rds$", "-seg.rds")
        }
    }
    saveRDS(cnv, args$out)
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
                              default="gene",
                              help="Output type"),
        optparse::make_option(c("-i", "--inp"), type="character",
                              default=NULL,
                              help="Input CNVEX file"),
        optparse::make_option(c("-o", "--out"), type="character",
                              default=NULL,
                              help="Output CSV file"),
        optparse::make_option(c("-s", "--suffix"), type="character",
                              default=NULL,
                              help="Optional file name suffix")
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

    if(is.null(args$inp) || !file.exists(args$inp)) {
        optparse::print_help(parser)
        write("Input file not found.\n", stderr())
        quit("no", 1)
    }
    cnv <- readRDS(args$inp)

    if (args$type == "gene") {
        csv <- .dtRound(.geneCSV(cnv, opts))
    }
    else if (args$type == "segment") {
        csv <- .dtRound(.segmentCSV(cnv, opts))
    }
    else {
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
plot <- function() {
    suppressPackageStartupMessages(library(ggpubr))
    suppressPackageStartupMessages(library(ggsci))
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

    config <- resolveConfig(args$config)
    if (config=="") {
        optparse::print_help(parser)
        write("Config file not found.\n", stderr())
        quit("no", 1)
    }
    opts <- getOpts(config, opts=list(cores=args$cores))

    if(is.null(args$inp) || !file.exists(args$inp)) {
        optparse::print_help(parser)
        write("Input file not found.\n", stderr())
        quit("no", 1)
    }
    cnv <- readRDS(args$inp)
    
    if (args$type=="gc") {
        plt <- plotGC(cnv)
        if (is.null(args$out)) {
            args$out <- str_replace(args$inp, ".rds$", "-gc.png")
        }
        dummy <- capture.output({
            pdf(NULL)
            suppressMessages(ggsave(args$out, plt, width=15, height=5))
            dev.off()
        })
    }
    if (args$type=="cnv") {
        plt <- plotCNV(cnv, opts, sel.lr = args$lr, sel.chr = args$chr)
        if (is.null(args$out)) {
            if (is.null(args$chr)) {
                args$out <- str_replace(args$inp, ".rds$", sprintf("-cnv-%s.png", args$lr))
            } else {
                args$out <- str_replace(args$inp, ".rds$", sprintf("-cnv-%s-%s.png", args$lr, args$chr))
            }
        }
        dummy <- capture.output({
            pdf(NULL)
            suppressMessages(ggsave(args$out, plt, width=15, height=5))
            dev.off()
        })
    }
    if (args$type=="seg") {
        plt <- plotSeg(cnv, opts, sel.chr = args$chr)
        if (is.null(args$out)) {
            if (is.null(args$chr)) {
                args$out <- str_replace(args$inp, ".rds$", "-seg.png")
            } else {
                args$out <- str_replace(args$inp, ".rds$", sprintf("-seg-%s.png", args$chr))
            }
        }
        dummy <- capture.output({
            pdf(NULL)
            suppressMessages(ggsave(args$out, plt, width=15, height=5))
            dev.off()
        })
    }
}
