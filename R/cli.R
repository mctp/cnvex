importCNVEX <- function(opts, tn.vcf, t.bam, n.bam) {
    tile <- getTile(opts)
    var <- importVcf(tn.vcf, tile, opts)
    tile <- rawLogRatio(t.bam, n.bam, tile, opts)
    cnv <- list(var=var, tile=tile)
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
          "Written by Marcin Cieslik (mcieslik@med.umich.edu) ",
          "Michigan Center for Translational Pathology (c) 2018\n"),
      option_list=option_list
      )

    args <- optparse::parse_args(parser, positional_arguments=FALSE)
    opts <- getOpts(args$config, list(cores=args$cores))
    cnvex <- importCNVEX(opts, args$vcf, args$tumor, args$normal)
    saveRDS(cnvex, args$out)
}
