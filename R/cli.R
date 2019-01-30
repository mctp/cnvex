#' @export
load <- function() {

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
                              help="Target Somatic VCF file"),
        optparse::make_option(c("-o", "--out"), type="character",
                              default="cnvex.rds",
                              help="Output file")
        
    )
    parser = optparse::OptionParser(
      "Rscript -e 'cnvex::import()' [options]",
      description=c("Import data to create CNVEX object.\n"),
      epilogue=c(
          "Written by Marcin Cieslik (mcieslik@med.umich.edu)",
          "Michigan Center for Translational Pathology (c) 2018\n"),
      option_list=option_list
      )

    args <- optparse::parse_args(parser, positional_arguments=FALSE)
    args$vcf <- c(genome=args$gvcf, target=args$evcf)
    args$tumor <- str_split(args$tumor, ",|;")[[1]]
    args$normal <- str_split(args$normal, ",|;")[[1]]
    
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
    cnvex <- importCNVEX(args$vcf, args$tumor, args$normal, opts)
    saveRDS(cnvex, args$out)
}

#' @export
pool <- function() {

    option_list = list(
        optparse::make_option(c("-c", "--config"), type="character",
                              default="genome",
                              help="Preset file"),
        optparse::make_option(c("-p", "--cores"), type="integer",
                              default=4L,
                              help="Number of cores"),
        optparse::make_option(c("-o", "--out"), type="character",
                              default="pool.rds",
                              help="Output file")
    )
    parser = optparse::OptionParser(
      "Rscript -e 'cnvex::pool()' [options] files",
      description=c("Create pool of normals CNVEX object.\n"),
      epilogue=c(
          "Written by Marcin Cieslik (mcieslik@med.umich.edu)",
          "Michigan Center for Translational Pathology (c) 2018\n"),
      option_list=option_list
      )

    args <- optparse::parse_args(parser, positional_arguments=TRUE)

    config <- resolveConfig(args$options$config)
    if (config=="") {
        optparse::print_help(parser)
        write("Config file not found.\n", stderr())
        quit("no", 1)
    }
    opts <- getOpts(config, opts=list(cores=args$options$cores))
    if (length(args$args) < 2) {
        optparse::print_help(parser)
        write("Multiple CNVEX input files required.\n", stderr())
        quit("no", 1)
    }
    
    pool <- createPool(args$args, opts)
    saveRDS(pool, args$options$out)
}


#' @export
analyze <- function() {

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
        optparse::make_option(c("-l", "--pool"), type="character",
                              default=NULL,
                              help="Pool of normals"),
        optparse::make_option(c("-o", "--out"), type="character",
                              default=NULL,
                              help="Output CNVEX file"),
        optparse::make_option(c("-s", "--suffix"), type="character",
                              default=NULL,
                              help="Optional file name suffix")
    )
    parser = optparse::OptionParser(
      "Rscript -e 'cnvex::analyze()' [options]",
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
    if (is.null(args$inp) || !file.exists(args$inp)) {
        optparse::print_help(parser)
        write("Input file not found.\n", stderr())
        quit("no", 1)
    }
    cnv <- readRDS(args$inp)

    if (!is.null(args$pool)) {
        if (file.exists(args$pool)) {
            pool <- readRDS(args$pool)
        } else {
            optparse::print_help(parser)
            write("Pool file not found.\n", stderr())
            quit("no", 1)
        }
    } else {
        pool <- NULL
    }
    
    cnv <- addLogRatio(cnv, pool, opts)
    cnv <- addSegment(cnv, opts)
    
    if (is.null(args$out)) {
        if (is.null(args$suffix)) {
            args$out <- str_replace(args$inp, ".rds$", "-seg.rds")
        } else {
            args$out <- str_replace(args$inp, ".rds$", sprintf("-%s-seg.rds", args$suffix))
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
        plt <- plotCNV(cnv, opts, sel.chr = args$chr)
        if (is.null(args$out)) {
            if (is.null(args$chr)) {
                args$out <- str_replace(args$inp, ".rds$", "-cnv.png")
            } else {
                args$out <- str_replace(args$inp, ".rds$", sprintf("-cnv-%s.png", args$chr))
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
