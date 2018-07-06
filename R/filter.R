filterGenomeGermlineHets <- function(var, tile, opts) {
    min.qual <- quantile(var$QUAL, 0.20, na.rm=TRUE)
    var <- var[
        ## germline
        !var$SOMATIC &
        ## simple variant
        var$TYPE=="SNV" &
        ## heterozygous
        var$n.GT %in% c("0/1", "1/0") &
        var$t.GT %in% c("0/1", "1/0") &
        ## right range
        var$n.AF > 0.25 & var$n.AF < 0.75 &
        ## high-quality
        (!var$mask.strict | var$QUAL > min.qual)
    ]
    return(var)
}

filterTargetGermlineHets <- function(var, tile, opts) {
    splash <- tile[tile$target] + opts$shoulder
    var <- var[
        var %over% splash &
        ## germline
        !var$SOMATIC &
        ## simple variant
        var$TYPE=="SNV" &
        ## heterozygous
        var$n.GT %in% c("0/1", "1/0") &
        var$t.GT %in% c("0/1", "1/0") &
        ## right range
        var$n.AF > 0.25 & var$n.AF < 0.75 &
        ## high 
        var$n.DP > 16 &
        var$t.DP > 16 &
        ##
        (!var$mask.strict)
    ]
    return(var)
}


filterGermlineCont <- function(var) {
    var <- var[
        ## not likely artifact
        !var$mask.strict &
        ## germline
        !var$SOMATIC &
        ## simple variant
        var$TYPE=="SNV" &
        ## heterozygous
        var$n.GT %in% c("0/1", "1/0") &
        ## contamination level
        var$n.AF<0.2
    ]
    return(var)
}

filterSomatic <- function(var, qual.cut=0.2, min.n.dp=10, max.n.alt=1) {
    min.qual <- quantile(var$QUAL, qual.cut, na.rm=TRUE)
    var <- var[
        ## not likely artifact
        !var$mask.strict &
        ## somatic
        var$SOMATIC & round(var$n.AF*var$n.DP)<=max.n.alt &
        ## enough coverage
        var$n.DP>min.n.dp &
        ## high quality
        var$QUAL > min.qual
    ]
    return(var)
}
