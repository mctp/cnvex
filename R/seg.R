.len.penalty <- function(x) {
    1/x
}

.runJointSeg <- function(arm, method, tile.width, len.min, cbs.lr, cbs.baf, rbs.selection, only.target) {
    arm.lr <- arm$lr.smooth
    arm.baf <- arm$baf
    if (only.target) {
        arm.lr[!arm$target] <- NA_real_
        arm.baf[!arm$target] <- NA_real_
    }
    ## make sure we have enough points to segment
    n.points <- sum(!is.na(arm.lr) & !is.na(arm.baf))
    seg0 <- integer()
    if (n.points >= 2*len.min) {
        if (method %in% c("RBS", "DynamicProgramming")) {
            tile.per.mb <- 1e6 / tile.width
            opt.K <- ceiling(n.points / tile.per.mb)
            seg0 <- suppressWarnings(
                jointSeg(cbind(arm.lr, arm.baf), method=method,
                         modelSelectionMethod=rbs.selection, K=opt.K)$bestBkp
            )
        } else if (method=="CBS") {
            seg0.lr <- .runCBS(arm.lr, cbs.lr)
            seg0.baf <- .runCBS(arm.baf, cbs.baf)
            seg0 <- unique(sort(c(seg0.lr, seg0.baf)))
        } else {
            stop("Joint segmentation method not supported.")
        }
    }
    return(seg0)
}

.mergeBreakpoints <- function(seg0, sd.lr, sd.baf, arm) {
    best1 <- integer()
    ## got at least one breakpoint
    if (length(seg0) > 0) {
        ## compute breakpoint stats
        beg0 <- c(1, seg0+1)
        end0 <- c(seg0, length(arm))
        len0 <- diff(c(0, end0))
        idx0 <- rep(seq_along(beg0), len0)
        lr0 <- split(arm$lr.smooth, idx0)
        baf0 <- split(arm$baf, idx0)
        stat0 <- data.table(
            seg=seg0,
            lr.diff =sapply(2:length(beg0), function(i) .absMedDiff( lr0[[i]],  lr0[[i-1]])),
            baf.diff=sapply(2:length(beg0), function(i) .absMedDiff(baf0[[i]], baf0[[i-1]])),
            min.len =sapply(2:length(beg0), function(i) min(length(lr0[[i]]), length(lr0[[i-1]])))
        )
        stat0[is.na(lr.diff), lr.diff:=0]
        stat0[is.na(baf.diff), baf.diff:=0]
        stat0[,len.penalty:=.len.penalty(min.len)]
        ## pick weakest breakpoint to merge
        stat0 <- stat0[order(lr.diff, baf.diff)]
        seg1 <- stat0[(
             lr.diff <  sd.lr +  sd.lr * len.penalty &
            baf.diff < sd.baf + sd.baf * len.penalty
        ), seg]
        best1 <- head(seg1, 1)
    }
    return(best1)
}

.jointSegArm <- function(arm, sd.lr, sd.baf, method, tile.width, len.min, cbs.lr, cbs.baf, rbs.selection, sd.prune, len.prune) {
    ## initial segmentation
    seg1 <- .runJointSeg(arm, method, tile.width, len.min, cbs.lr, cbs.baf, rbs.selection, TRUE)
    ## iteratively remove spurious breakpoints
    if (sd.prune) {
        repeat({
            best1 <- .mergeBreakpoints(seg1, sd.lr, sd.baf, arm)
            seg1 <- setdiff(seg1, best1)
            if (length(best1) == 0) {break}
        })
    }
    ## skip short segments
    if (len.prune) {
        end1 <- c(seg1, length(arm))
        len1 <- diff(c(0, end1))
        seg1 <- head(end1[len1>=len.min], -1)
    }
    ## append segmentation
    bpt1 <- c(1, seg1 + 1)
    len1 <- diff(c(bpt1, length(arm)+1))
    idx1 <- rep(seq_along(bpt1), len1)
    arm$seg <- idx1
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
    sd.lr <- estimateSd(cnv$tile$lr.smooth) * opts$seg.sd.lr.penalty
    sd.baf <- estimateSd(cnv$tile$baf) * opts$seg.sd.baf.penalty
    ## create segmentation
    cnv$tile <- .addJointSeg(cnv$tile, sd.lr, sd.baf, opts$seg.method,
                             opts$tile.width, opts$seg.len.min, opts$seg.cbs.lr,
                             opts$seg.cbs.baf, opts$seg.rbs.selection, opts$seg.sd.prune, opts$seg.len.prune)
    cnv$seg <- unname(unlist(range(split(cnv$tile, cnv$tile$seg))))
    return(cnv)
}

getSeg <- function(cnv, opts) {
    ## seg -> var
    cnv$var$seg <- findOverlaps(cnv$var, cnv$seg, select="first", maxgap = opts$tile.shoulder-1)
    ## seg -> gene
    hit.prom <- as.data.table(findOverlaps(promoters(cnv$gene, 0, 1), cnv$seg))
    hit.gene <- as.data.table(findOverlaps(cnv$gene, cnv$seg))
    hit.gene <- hit.gene[!(queryHits %in% hit.prom$queryHits)]
    hit.gene <- hit.gene[,.SD[1],by=queryHits]
    hit <- rbind(hit.prom, hit.gene)
    cnv$gene$seg <- NA_real_
    cnv$gene$seg[hit$queryHits] <- hit$subjectHits
    return(cnv)
}
