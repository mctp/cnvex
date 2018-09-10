plotGC <- function(cnv) {
    tmp <- as.data.table(mcols(cnv$tile))
    tmp <- melt(tmp[,.(gc, blacklist, target, lr.raw, lr.gc, lr.off=lr.raw-lr.gc)], id.vars=c("gc", "blacklist", "target"))
    tmp[,delta:=FALSE]
    tmp[variable=="lr.off", delta:=TRUE]
    tmp[variable=="lr.off", variable:="lr.raw"]
    tmp[,target.lab:=ifelse(target, "on-target", "off-target")]
    tmp[,lr.lab:=ifelse(variable=="lr.raw", "raw", "adjusted")]
    plt <- ggplot(tmp) + aes(x=gc, y=value, color=delta) +
        facet_grid(lr.lab~target.lab) +
        geom_point(alpha=0.05, size=0.1) +
        scale_color_manual(values=c("black", "red"), guide=FALSE) +
        scale_x_continuous(labels=scales::percent) +
        coord_cartesian(xlim=c(0.2, 0.8), ylim=c(-3, 3)) +
        ylab("log2(tumor/normal)") +
        xlab("GC [%]") +
        theme_pubr()
    return(plt)
}

plotCNV <- function(cnv, sel.lr="lr.smooth", sel.chr=NULL) {
    if (is.null(sel.chr)) {
        sel.chr <- paste("chr", c(1:22, "X", "Y"), sep="")
    }
    germline <- filterGermlineHets(cnv$tile, cnv$var, opts)
    snp <- cnv$var[germline]
    snp <- snp[seqnames(snp) %in% sel.chr]
    cov <- cnv$tile
    cov <- cov[seqnames(cov) %in% sel.chr]

    snp.dt <- data.table(
        chr=as.character(seqnames(snp)),
        pos=start(snp),
        val=snp$t.AF,
        type="BAF"
    )[!is.na(val)]
    snp.dt[,chr:=factor(chr, sel.chr, ordered=TRUE)]
    
    cov.dt <- data.table(
        chr=as.character(seqnames(cov)),
        pos=floor((start(cov)+end(cov))/2),
        val=mcols(cov)[[sel.lr]],
        type="COV",
        tgt=mcols(cov)[["target"]]
    )[!is.na(val)]
    cov.dt[,chr:=factor(chr, sel.chr, ordered=TRUE)]
    cov.dt <- cov.dt[-c(1,nrow(cov.dt))][order(-tgt)]
    ## fix
    snp.dt <- rbind(snp.dt, cov.dt[,.SD[1],by=chr][,.(chr, pos, val=0.5, type="BAF")])
    cov.dt <- rbind(cov.dt, snp.dt[,.SD[1],by=chr][,.(chr, pos, val=0.0, type="COV", tgt=TRUE)])

    ## COV
    if (length(unique(cov.dt$tgt))==2) {
        cols <- c("red", "black")
    } else {
        cols <- "black"
    }
    cov.plt <- ggplot(cov.dt) + facet_grid(.~chr, scales="free_x") +
        aes(x=pos, y=val, color=tgt) +
        scale_color_manual(values=cols, guide=FALSE) +
        geom_point(size=0.5, alpha=0.5, shape=16) +
        coord_cartesian(ylim=c(min(-3, min(cov.dt$val, na.rm=TRUE)),
                               max (3, max(cov.dt$val, na.rm=TRUE)))) +
        coord_cartesian(ylim=c(-3,3)) +
        geom_hline(yintercept = 0, color="blue") +
        theme_pubr() +
        ylab("log2(tumor/normal)") +
        theme(
            panel.spacing = unit(0, "lines"),
            axis.title.x=element_blank(),
            axis.ticks.x = element_blank(),
            axis.text.x = element_blank()
        )
    
    ## SNP
    if (nrow(snp.dt)>1e5) {
        p.size=0.2
        p.alpha=0.2
    } else {
        p.size=0.5
        p.alpha=0.5
    }
    snp.plt <- ggplot(snp.dt) + facet_grid(.~chr, scales="free_x") +
        aes(x=pos, y=val) +
        geom_point(size=p.size, alpha=p.alpha, shape=16) +
        coord_cartesian(ylim=c(0,1)) + 
        theme_pubr() +
        ylab("BAF") +
        theme(
            panel.spacing = unit(0, "lines"),
            axis.title.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.text.x = element_blank()
        )
    
    plt <- rbind(ggplotGrob(cov.plt), ggplotGrob(snp.plt), size = "last")
    return(plt)
}

plotSeg <- function(cnv, opts, sel.lr="lr.smooth", sel.chr=NULL) {
    if (is.null(sel.chr)) {
        sel.chr <- paste("chr", c(1:22, "X", "Y"), sep="")
    }
    germline <- filterGermlineHets(cnv$tile, cnv$var, opts)
    snp <- cnv$var[germline]
    snp <- snp[seqnames(snp) %in% sel.chr]
    cov <- cnv$tile
    cov <- cov[seqnames(cov) %in% sel.chr]

    snp.dt <- data.table(
        chr=as.character(seqnames(snp)),
        pos=start(snp),
        val=snp$t.AF,
        seg=snp$seg,
        type="BAF"
    )[!is.na(val)]
    snp.dt[,chr:=factor(chr, sel.chr, ordered=TRUE)]
    
    cov.dt <- data.table(
        chr=as.character(seqnames(cov)),
        pos=floor((start(cov)+end(cov))/2),
        val=mcols(cov)[[sel.lr]],
        type="COV",
        seg=cov$seg,
        tgt=mcols(cov)[["target"]]
    )[!is.na(val)]
    cov.dt[,chr:=factor(chr, sel.chr, ordered=TRUE)]
    cov.dt <- cov.dt[-c(1,nrow(cov.dt))][order(-tgt)]
    ## fix
    snp.dt <- rbind(snp.dt, cov.dt[,.SD[1],by=chr][,.(chr, pos, val=0.5, seg=0, type="BAF")])
    cov.dt <- rbind(cov.dt, snp.dt[,.SD[1],by=chr][,.(chr, pos, val=0.0, seg=0, type="COV", tgt=TRUE)])

    ## COV
    cov.plt <- ggplot(cov.dt) + facet_grid(.~chr, scales="free_x") +
        aes(x=pos, y=val, color=factor(ifelse(!tgt,4,as.integer(seg%%3)))) +
        scale_color_npg(guide=FALSE) +        
        geom_point(size=0.75, alpha=0.5, shape=16) +
        coord_cartesian(ylim=c(min(-3, min(cov.dt$val, na.rm=TRUE)),
                               max (3, max(cov.dt$val, na.rm=TRUE)))) +
        coord_cartesian(ylim=c(-3,3)) +
        geom_hline(yintercept = 0, color="blue") +
        theme_pubr() +
        ylab("log2(tumor/normal)") +
        theme(
            panel.spacing = unit(0, "lines"),
            axis.title.x=element_blank(),
            axis.ticks.x = element_blank(),
            axis.text.x = element_blank()
        )
    
    ## SNP
    if (nrow(snp.dt)>1e5) {
        p.size=0.2
        p.alpha=0.5
    } else {
        p.size=0.5
        p.alpha=0.75
    }
    snp.plt <- ggplot(snp.dt) + facet_grid(.~chr, scales="free_x") +
        aes(x=pos, y=val, color=factor(as.integer(seg%%3))) +
        geom_point(size=p.size, alpha=p.alpha, shape=16) +
        coord_cartesian(ylim=c(0,1)) +
        scale_color_npg(guide=FALSE) +
        theme_pubr() +
        ylab("BAF") +
        theme(
            panel.spacing = unit(0, "lines"),
            axis.title.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.text.x = element_blank()
        )
    
    plt <- rbind(ggplotGrob(cov.plt), ggplotGrob(snp.plt), size = "last")
    return(plt)
}

## plotSeg <- function(cnv, sel.chr=NULL, sel.lr="lr") {
##     gr <- unlist(cnv$tile)
##     l2r <- ifelse(is.finite(mcols(gr)[[sel.lr]]), mcols(gr)[[sel.lr]], NA_real_)
##     baf <- mcols(gr)[["baf"]]
##     tmp <- data.table(
##         chr=as.character(seqnames(gr)),
##         pos=floor((start(gr)+end(gr))/2),
##         l2r=l2r,
##         baf=baf,
##         seg=mcols(gr)[["seg"]],
##         tgt=mcols(gr)[["target"]]
##     )
##     print(tmp)
##     if (!is.null(sel.chr)) {
##         tmp <- tmp[chr %in% sel.chr]
##     }
##     plt.l2r <- ggplot(tmp) + facet_grid(chr~.) +
##         aes(x=pos, y=l2r, color=) +
##         geom_point(size=0.5) +
##         coord_cartesian(ylim=c(min(-3, min(tmp$l2r, na.rm=TRUE)),
##                                max (3, max(tmp$l2r, na.rm=TRUE)))) +
##         
##         theme_pubr()
    
##     plt.baf <- ggplot(tmp) + facet_grid(chr~.) +
##         aes(x=pos, y=baf,
##             color=
##             ) +
##         geom_point(size=0.5) +
##         coord_cartesian(ylim=c(0,0.5))+ 
##         theme_pubr()
    
##     plt <- grid.arrange(plt.l2r, plt.baf, ncol=1)
##     return(plt)
## }

plotGrid <- function(grid, cand) {
    midp <- min(cand$L1)
    delta <- (max(cand$L1) - midp)
    lowp <- midp - delta
    tmp <- copy(grid)
    tmp[L1<lowp, L1:=lowp]
    plt <- ggplot(tmp) +
        aes(x=p0, y=D0, fill=L1) +
        geom_tile() +
        scale_fill_gradient2(low="blue", mid="white", high="red",
                             name="likelihood", midpoint=midp
                             ) +  
        geom_point(data=cand, size=2, color="black") +
        theme_pubr(legend="right")
    return(plt)
}

plotFine <- function(grid, fine) {
    plt <- ggplot() +
        geom_tile(aes(x=p0, y=D0, fill=ifelse(abs(L1)>max(L1), -max(L1), L1)), data=grid) +
        scale_fill_gradient2(low="blue", mid="white", high="red",
                             name="likelihood"
                             ) +
        geom_segment(aes(x=p0, y=D0, xend=pi, yend=Di), fine, size=1.2) +
        coord_cartesian(ylim=c(1,6), xlim=c(0.05,0.95)) +
        theme_pubr(legend="right")
}


plotStsC <- function(sts, sel.chr=NULL) {
    gr <- unlist(sts$tile)
    tmp <- data.table(
        chr=as.character(seqnames(gr)),
        pos=floor((start(gr)+end(gr))/2),
        l2r=mcols(gr)[["lr"]],
        baf=mcols(gr)[["baf"]],
        C=mcols(gr)[["C"]]
    )
    if (!is.null(sel.chr)) {
        tmp <- tmp[chr %in% sel.chr]
    }
    plt.l2r <- ggplot(tmp) + facet_grid(chr~.) +
        aes(x=pos, y=l2r, color=factor(C)) +
        geom_point(size=0.5) +
        scale_color_npg() +
        theme_pubr()
    
    plt.baf <- ggplot(tmp) + facet_grid(chr~.) +
        aes(x=pos, y=baf, color=factor(C))+
        geom_point(size=0.5) +
        scale_color_npg(guide=FALSE) +
        theme_pubr()
    
    plt <- grid.arrange(plt.l2r, plt.baf, ncol=1)
    return(plt)
}

