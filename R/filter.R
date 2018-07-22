filterGenomeGermlineHets <- function(var, tile, opts) {
    min.qual <- quantile(var$QUAL, 0.20, na.rm=TRUE)
    var$germline <- 
        ## germline
        !var$SOMATIC &
        ## simple variant
        var$TYPE=="SNV" &
        ## heterozygous
        var$n.GT %in% c("0/1", "1/0") &
        ## right range
        var$n.AF > 0.25 & var$n.AF < 0.75 &
        ## high-quality
        (!var$mask.strict | var$QUAL > min.qual)
    return(var)
}

filterTargetGermlineHets <- function(var, tile, opts) {
    splash <- tile[tile$target] + opts$shoulder
    var$germline <- 
        var %over% splash &
        ## germline
        !var$SOMATIC &
        ## simple variant
        var$TYPE=="SNV" &
        ## heterozygous
        var$n.GT %in% c("0/1", "1/0") &
        ## right range
        var$n.AF > 0.25 & var$n.AF < 0.75 &
        ## high 
        var$n.DP > 16 &
        var$t.DP > 16 &
        ## mappability
        (!var$mask.strict)
    return(var)
}
