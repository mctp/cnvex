plotGC <- function(sts) {
    tbl <- as.data.table(mcols(cnv$tile))
    ## off-target
    tmp <- tbl[(!target)]
    tmp <- melt(tmp[,.(gc, blacklist, lr.raw, lr, lr.off=lr.raw-lr)], id.vars=c("gc", "blacklist"))
    tmp[,delta:=FALSE]
    tmp[variable=="lr.off", delta:=TRUE]
    tmp[variable=="lr.off", variable:="lr.raw"]    
    plt <- ggplot(tmp) + aes(x=gc, y=value, color=delta) + facet_grid(variable~.) +
        geom_point(alpha=0.1) + 
        scale_color_manual(values=c("black", "red")) +
        coord_cartesian(xlim=c(0.2, 0.8)) +
        theme_pubr()

    tmp <- tbl[(target & !is.na(baf))]
    tmp <- melt(tmp[,.(gc, blacklist, lr.raw, lr, lr.off=lr.raw-lr)], id.vars=c("gc", "blacklist"))
    tmp[,delta:=FALSE]
    tmp[variable=="lr.off", delta:=TRUE]
    tmp[variable=="lr.off", variable:="lr.raw"]    
    plt <- ggplot(tmp) + aes(x=gc, y=value, color=delta) + facet_grid(variable~.) +
        geom_point(alpha=0.1) +
        scale_color_manual(values=c("black", "red")) +
        coord_cartesian(xlim=c(0.2, 0.8)) +
        theme_pubr()
    ggsave("blah.png", plt)


    tmp <- tbl[(target & !is.na(baf))]
    plt <- ggplot(tmp) + aes(x=lr.raw, y=lr, color=gc) + 
        geom_point(alpha=0.1) +
        geom_smooth(method="lm")+
        scale_color_gradient2(midpoint=0.49) +
        theme_pubr()
    ggsave("blah.png", plt)    

    
}

plotSeg <- function(sts, sel.chr=NULL, sel.lr="lr") {
    gr <- unlist(sts$tile)
    l2r <- ifelse(is.finite(mcols(gr)[[sel.lr]]), mcols(gr)[[sel.lr]], NA_real_)
    baf <- mcols(gr)[["baf"]]
    tmp <- data.table(
        chr=as.character(seqnames(gr)),
        pos=floor((start(gr)+end(gr))/2),
        l2r=l2r,
        baf=baf,
        seg=mcols(gr)[["seg"]],
        tgt=mcols(gr)[["target"]]
    )
    print(tmp)
    if (!is.null(sel.chr)) {
        tmp <- tmp[chr %in% sel.chr]
    }
    plt.l2r <- ggplot(tmp) + facet_grid(chr~.) +
        aes(x=pos, y=l2r, color=factor(ifelse(!tgt,4,as.integer(seg%%3)))) +
        geom_point(size=0.5) +
        coord_cartesian(ylim=c(min(-3, min(tmp$l2r, na.rm=TRUE)),
                               max (3, max(tmp$l2r, na.rm=TRUE)))) +
        scale_color_npg(guide=FALSE) +
        theme_pubr()
    
    plt.baf <- ggplot(tmp) + facet_grid(chr~.) +
        aes(x=pos, y=baf,
            color=factor(as.integer(seg%%3))
            ) +
        geom_point(size=0.5) +
        coord_cartesian(ylim=c(0,0.5))+ 
       scale_color_npg(guide=FALSE) +
        theme_pubr()
    
    plt <- grid.arrange(plt.l2r, plt.baf, ncol=1)
    return(plt)
}

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

