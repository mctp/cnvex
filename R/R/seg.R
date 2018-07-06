.smooth.outliers <- function(data, chr, pos) {
    print(pos)
    obj <- CNA(data, as.integer(factor(chr)), pos,
              data.type = "logratio",
              sampleid = "sample")
    adj <- smooth.CNA(obj)$sample
    
}

.smooth.outliers.gr <- function(gr, data.col) {
    obj <- CNA(mcols(gr)[[data.col]], as.integer(seqnames(gr)), floor((start(gr)+end(gr))/2),
              data.type = "logratio",
              sampleid = "sample")
    adj <- smooth.CNA(obj)$sample
    return(adj)
}

.len.penalty <- function(x) {
    1/x
}

.jointSegArm <- function(arm, K, min.seg, sd.lr, sd.baf, only.target=TRUE) {
    ## make sure we have enough points to segment
    max.K <- sum(!is.na(arm$lr) & !is.na(arm$baf))/10
    K <- min(K, max.K)
    if (K>1) {
        ## initial segmentation
        arm.lr <- arm$lr
        arm.baf <- arm$baf
        if (only.target) {
            arm.lr[!arm$target] <- NA_real_
            arm.baf[!arm$target] <- NA_real_
        }
        seg0 <- jointSeg(cbind(arm.lr, arm.baf), K=K)$bestBkp
        ## skip short segments
        bpt0 <- c(0, seg0) + 1
        len0 <- diff(bpt0)
        seg1 <- seg0[len0>=min.seg]
        if (length(seg1)>0) {
            ## corner-
            if ((length(arm)-tail(seg1, 1))<min.seg) {
                seg1 <- head(seg1, -1)
            }
            if (length(seg1)>0) {
                bpt1 <- c(1, seg1 + 1)
                len1 <- diff(c(bpt1, length(arm)+1))
                idx1 <- rep(seq_along(bpt1), len1)
                ## compute breakpoints stats
                lr1 <- split(arm$lr, idx1)
                baf1 <- split(arm$baf, idx1)
                stat1 <- data.table(
                    seg1=seg1,
                    lr.diff =sapply(2:length(bpt1), function(i) .absMedDiff( lr1[[i]],  lr1[[i-1]])),
                    baf.diff=sapply(2:length(bpt1), function(i) .absMedDiff(baf1[[i]], baf1[[i-1]])),
                    min.len =sapply(2:length(bpt1), function(i) min(length(lr1[[i]]), length(lr1[[i-1]])))
                )
                stat1[is.na(lr.diff), lr.diff:=0]
                stat1[is.na(baf.diff), baf.diff:=0]
                stat1[,len.penalty:=.len.penalty(min.len)]
                
                ## filtered segmentation
                seg2 <- stat1[(
                    lr.diff  > sd.lr  + sd.lr  * len.penalty |
                    baf.diff > sd.baf + sd.baf * len.penalty
                ), seg1]
                bpt2 <- c(1, seg2 + 1)
                len2 <- diff(c(bpt2, length(arm)+1))
                idx2 <- rep(seq_along(bpt2), len2)
                arm$seg <- idx2
            } else {
                arm$seg <- 1L
            }
        } else {
            arm$seg <- 1L
        }
    } else {
        arm$seg <- 1L
    }
    return(arm)
}

.addJointSeg <- function(gt, ...) {
    gt <- sort(unname(unlist(endoapply(split(gt, gt$arm), .jointSegArm, ...))))
    ## provide globally unique ids
    tmp <- paste(gt$arm, gt$seg)
    gt$seg <- as.integer(factor(tmp, levels=unique(tmp)))
    return(gt)
}

jointSegment <- function(cnv, opts) {
    ## estimated standard-deviation
    sd.lr <- estimateSd(cnv$tile$lr) * opts$sd.penalty
    sd.baf <- estimateSd(cnv$tile$baf) * opts$sd.penalty
    ## create segmentation
    cnv$tile <- .addJointSeg(cnv$tile, opts$K, opts$min.seg, sd.lr, sd.baf)
    cnv$seg <- unname(unlist(range(split(cnv$tile, cnv$tile$seg))))
    cnv$snp$seg <- findOverlaps(cnv$snp, cnv$seg, select="first", maxgap = opts$shoulder-1)
    return(cnv)
}
