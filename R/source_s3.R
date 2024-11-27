print.grafpop <- function(x, ...) {

    cat("\nPredicted reference population counts:")
    tab <- table(x$table[, "Refpop", drop=TRUE], exclude=NULL)
    tmp <- gp_getReturnRefData(x)
    if (!length(tmp)) {
      nms <- names(tab)
      ref <- gp_getOrder()$refpop
      tmp <- ref %in% nms
      ref <- ref[tmp]  
      print(tab[ref])
    } else {
      print(tab) 
    }
    invisible(x)
}

plot.grafpop <- function(x, legend.pos="right", showRefData=TRUE,
                         refObj=NULL, ...){

    plotGrafOut(x, legend.pos=legend.pos, 
            showRefData=showRefData, refObj=refObj)
    invisible(NULL)
}

plotGrafOut <- function(obj, xvar="GD1_x", yvar="GD2_y", 
                        legend.pos="right", showRefData=TRUE, refObj=NULL) {

    interactive <- FALSE
    triangle    <- 1
    if ((xvar != "GD1_x") || (yvar != "GD2_y")) triangle <- 0 

    hpflag   <- obj$objects$hpflag
    refcols  <- obj$objects$ref.pops
    vpops    <- obj$objects$vertex.pops
    chr.col  <- obj$objects$chr.col
    prophage <- obj$objects$prophage 
    data     <- obj$table
    tri      <- pgo_getVertex(obj)
    train_results <- getTrainResults(retobj=obj, refobj=refObj, prophage=prophage)
    X <- Y <- Refpop <- SampleID <- Refpop_n <- .data <- NULL

    # Set up the data
    if (hpflag) {
        tmp <- pgo_getData(data, train_results, xvar, yvar, legend.pos)
    } else {
        lab3vec <- paste0(vpops, "_percent")
        lab3    <- paste0(vpops, collapse=", ")
        lab3    <- paste0(lab3," Ancestry")
        tmp <- pgo_getData(data, train_results, xvar, yvar, legend.pos,
            lab3vec=lab3vec, lab3=lab3)  
    }
    nb.cols <- length(refcols)
    df  <- tmp$data
    train_results <- tmp$trn

    # Define the colors 
    mycolors <- colorRampPalette(brewer.pal(min(c(9, nb.cols)), "Set1"))(nb.cols)

    p1 <- pgo_ggplot(df, X, Y, xvar, yvar, legend.pos, mycolors)  
    if (showRefData) {
        # ellipse of training data
        p1 <- p1 + stat_ellipse(data=train_results, aes(fill=Refpop_n), 
            geom="polygon", alpha=0.25) + guides(fill="none")
    }
    if (triangle) p1 <- pgo_triangle(p1, tri)
    p1 <- p1 + geom_point(aes(color=Refpop_n), size=1)
    print(p1)
    NULL
}

pgo_ggplot <- function(df, X, Y, xvar, yvar, legend.pos, mycolors) {

    p1  <- ggplot(df, aes(X, Y)) +
        xlab(xvar) + ylab(yvar) +
        theme_classic() +       
        theme(panel.background = element_blank()) +
        theme(panel.border = element_rect(fill = NA, color="black")) +
        theme(legend.title = element_blank()) +
        theme(legend.text = element_text(size=10, family="serif")) +
        theme(strip.text = element_text(face="bold",family="serif",size=10)) +
        theme(axis.text=element_text(face="bold",family="serif", 
                                        color="black")) +
        theme(axis.title = element_text(face="bold",family="serif", 
                                        color="black")) +
        theme(plot.title = element_text(face = "bold",family="serif",
                                        color="black", hjust=0.5, size=12)) +
        theme(plot.subtitle=element_text(color="black", family="serif",
                                        hjust = 0.5,size=12)) +
        theme(legend.position=legend.pos) + 
        theme(legend.justification="right") + 
        scale_fill_manual(values = mycolors)+
        scale_color_manual(values = mycolors)+
        scale_x_continuous(xvar,breaks=pretty_breaks(n=10)) +
        scale_y_continuous(yvar,breaks=pretty_breaks(n=10)) 
 
    p1
}

pgo_triangle <- function(p1, tri) {

    p1 <- p1 + annotate("segment", x=tri$F.x, xend=tri$E.x, y=tri$F.y, 
                    yend=tri$E.y, color="black") +
    annotate("segment", x=tri$E.x, xend=tri$A.x, y=tri$E.y, yend=tri$A.y, 
                    color="black") +
    annotate("segment", x=tri$F.x, xend=tri$A.x, y=tri$F.y, yend=tri$A.y, 
                    color="black")
    p1
}

pgo_refactor <- function(vec) {

    ref  <- gp_getOrder()$refpop
    ref  <- paste0(ref, "\n")
    n    <- length(ref)
    levs <- NULL

    for (i in seq_len(n)) {
        str <- ref[i]    
        len <- nchar(str)
        tmp <- substr(vec, 1, len) == str
        tmp[is.na(tmp)] <- FALSE
        if (any(tmp)) {
            v2   <- vec[tmp]
            levs <- c(levs, v2[1])  
        }
    }
    if (length(levs)) {
        ret <- factor(vec, levels=levs)
    } else {
        ret <- factor(vec)
    }
    ret
}

pgo_getData <- function(data, train_results, xvar, yvar, legend.pos,
    lab3vec=c("A_percent", "F_percent", "E_percent"),
    lab3="African, European, Asian Ancestry") {

    tmp                  <- train_results
    train_results[, "X"] <- tmp[, xvar, drop=TRUE]
    train_results[, "Y"] <- tmp[, yvar, drop=TRUE]

    # Set up data frame
    df                 <- data
    df[, "SampleID"]   <- df[, "Sample", drop=TRUE]
    tmp                <- df
    df[, "X"]          <- tmp[, xvar, drop=TRUE]
    df[, "Y"]          <- tmp[, yvar, drop=TRUE]
    df[, "legend.pos"] <- legend.pos
    df[, "xvar"]       <- xvar
    df[, "yvar"]       <- yvar

    # Frequency counts, add to train_results also
    vec  <- df[, "Refpop", drop=TRUE]
    vec2 <- vec
    tvec <- train_results[, "Refpop", drop=TRUE]
    pops <- unique(tvec) 
    for (pop in pops) {
        tmp1      <- vec == pop
        m         <- sum(tmp1)
        if (m) vec[tmp1]  <- paste0(vec[tmp1], "\n n=", m)
        tmp       <- tvec == pop
        n         <- sum(tmp) 
        tvec[tmp] <- paste0(tvec[tmp], "\n n=", m)
        if (m) vec2[tmp1]  <- paste0(vec2[tmp1], "\n n=", n)
    }
    df[, "Refpop_n"]            <- pgo_refactor(vec)
    df[, "Refpop_n2"]           <- pgo_refactor(vec2)
    train_results[, "Refpop_n"] <- pgo_refactor(tvec)

    # hovering text
    df[, "label1"] <- paste0(df[, "Refpop", drop=TRUE], ", ", 
                        df[, "Sample", drop=TRUE])
    df[, "label2"] <- paste0(df[, "Nearest_neighbor", drop=TRUE], ", ", 
                        round(df[, "Separation_percent", drop=TRUE], 0))
    df[, "label3"] <- paste0(round(df[, lab3vec[1], drop=TRUE], 0), "%, ",
                        round(df[, lab3vec[2], drop=TRUE], 0), "%, ",
                        round(df[, lab3vec[3], drop=TRUE], 0), "%")
    df[, "Refpop, SampleID"]    <- df[, "label1", drop=TRUE]
    df[, "Nearest, Separation"] <- df[, "label2", drop=TRUE]
    df[, lab3]                  <- df[, "label3", drop=TRUE]

    list(data=df, trn=train_results)
}

pgo_getVertex <- function(obj) {

  def <- !length(gp_getReturnRefData(obj))
  if (def) {
    ret <- pgo_getVertex_def(obj)
  } else {
    vt  <- obj$vertex
    ret <- list()
    nms <- c("E", "F", "A")
    for (i in 1:3) {
      nm  <- nms[i]
      tmp <- vt[[i]]
      ret[[paste0(nm, ".x")]] <- tmp$x
      ret[[paste0(nm, ".y")]] <- tmp$y
    }
  }
  ret
}

pgo_getVertex_def <- function(obj) {

    vt       <- obj$vertex
    tmp      <- vt$Africa
    F.x      <- tmp$x
    F.y      <- tmp$y
    tmp      <- vt$Asia
    A.x      <- tmp$x
    A.y      <- tmp$y
    tmp      <- vt$Europe
    E.x      <- tmp$x
    E.y      <- tmp$y

    list(F.x=F.x, F.y=F.y, A.x=A.x, A.y=A.y, E.x=E.x, E.y=E.y)
}

grafGenPlot <- function(obj, which=1, legend.pos=NULL, 
                        ylim=NULL, showRefData=TRUE,
                        jitter=0, refObj=NULL) {

    check_grafpop(obj)
    which <- check_numVec(which, "which", valid=seq_len(5), def=1)
    if (length(ylim)) check_numVec(ylim, "ylim", len=2)
    legend.pos <- check_legendPos(legend.pos)
    check.logical(showRefData, "showRefData")

    ylim0 <- ylim
    lpos  <- legend.pos
    if (1 %in% which) {
        if (!length(lpos)) legend.pos <- "right"
        plotGrafOut(obj, xvar="GD1_x", yvar="GD2_y", legend.pos=legend.pos, 
            showRefData=showRefData, refObj=refObj) 
    }
    if (2 %in% which) {
        if (!length(lpos)) legend.pos <- "right"
        plotGrafOut(obj, xvar="GD1_x", yvar="GD3_z", legend.pos=legend.pos, 
            showRefData=showRefData, refObj=refObj) 
    }
    if (3 %in% which) {
        if (!length(lpos)) legend.pos <- "right"
        plotGrafOut(obj, xvar="GD2_y", yvar="GD3_z", legend.pos=legend.pos, 
            showRefData=showRefData, refObj=refObj) 
    }
    if (4 %in% which) {
        plotRefDist(obj, ylim=ylim0, showTrainData=showRefData, refObj=refObj)
    }
    if (5 %in% which) {
        plotRefPerc(obj, ylim=ylim0, showTrainData=showRefData, jitter=jitter,
            refObj=refObj)
    }

    invisible(NULL)
}

plotRefDist <- function(obj, ylim=NULL, showTrainData=TRUE, refObj=NULL) {

    legend.pos <- "top"
    prophage   <- obj$objects$prophage 
    train_results <- getTrainResults(retobj=obj, refobj=refObj, 
        prophage=prophage)

    # Get the reference pops
    trn_ref <- train_results[, "Refpop", drop=TRUE]
    pops    <- sort(unique(trn_ref))
    npop    <- length(pops)

    tst     <- obj$table
    tst_ref <- tst[, "Refpop", drop=TRUE]
    tst     <- as.matrix(tst[, pops, drop=FALSE])

    if (showTrainData) trn <- as.matrix(train_results[, pops, drop=FALSE])
    if (!length(ylim)) ylim <- plotRefDist_ylim(pops, showTrainData, trn_ref,
                                trn, tst_ref, tst)

    main <- "Reference Population Genetic Distances"
    plot(seq_len(npop), rep(0, npop), type="n", ylim=ylim, ylab="Distance",
        xlab="", xaxt="n", main=main)
    plotRefDist_points(pops, showTrainData, trn_ref, trn, tst_ref, tst)

    labs <- gsub("hpgp", "", pops, fixed=TRUE)
    axis(1, at=seq_len(npop), labels=labs, las=2, cex.axis=0.65, font=2)
    if (showTrainData) {
        legend(legend.pos,c("User Data","Training Data"),col=c("red","blue"),
            pch=c(20, 3), horiz=TRUE, bty="n")
    }
    NULL
}

plotRefDist_ylim <- function(pops, showTrainData, trn_ref, trn, 
                                tst_ref, tst) {

    npop <- length(pops)
    ymin <- Inf
    ymax <- -Inf
    for (i in seq_len(npop)) {
        pop  <- pops[i]
        tmp  <- tst_ref %in% pop
        if (any(tmp)) {
            y    <- tst[tmp, pop, drop=TRUE]
            ymin <- min(c(ymin, y))
            ymax <- max(c(ymax, y))
        }
        if (showTrainData) {
            tmp  <- trn_ref %in% pop
            y    <- trn[tmp, pop, drop=TRUE]
            ymin <- min(c(ymin, y))
            ymax <- max(c(ymax, y))
        }
    }   
    if (showTrainData) ymax <- ymax + (ymax - ymin)/10
    ylim <- c(ymin, ymax)

    ylim
}

plotRefDist_points <- function(pops, showTrainData, trn_ref, trn,
                                tst_ref, tst) {

    x.eps <- 0.05
    if (!showTrainData) x.eps <- 0
    npop  <- length(pops)
    for (i in seq_len(npop)) {
        pop <- pops[i]

        # Plot train first
        if (showTrainData) {
            tmp <- trn_ref %in% pop
            y   <- trn[tmp, pop, drop=TRUE]
            x   <- rep(i-x.eps, length(y))
            points(x, y, col="blue", pch=3, cex=0.75) 
        }    

        # Plot test
        tmp <- tst_ref %in% pop
        if (any(tmp)) { 
            y   <- tst[tmp, pop, drop=TRUE]
            x   <- rep(i+x.eps, length(y))
            points(x, y, col="red", pch=20, cex=0.75)
        } 
    }
    NULL
}

plotRefPerc <- function(obj, ylim=NULL, showTrainData=TRUE, jitter=0, refObj=NULL) {

    legend.pos <- "top"  
    clrs       <- c("red", "green", "blue")
    if (!length(ylim)) ylim <- c(0, 103)
    jeps     <- abs(jitter)
    objects  <- obj$objects
    prophage <- objects$prophage
    hpflag   <- objects$hpflag
    vpops    <- objects$vertex.pops
    rpops    <- objects$ref.pops
    user.ref <- objects$user.ref

    train_results <- getTrainResults(retobj=obj, refobj=refObj, prophage=prophage)
    if (hpflag && !prophage) {
        pvars <- c("F_percent", "A_percent", "E_percent")
        leg   <- c("Africa", "Asia", "Europe")
        NREF  <- 9
    } else if (prophage && !user.ref) {
        pvars <- c("F_percent", "A_percent", "E_percent")
        leg   <- c("Africa1", "EastAsia", "SWEurope")
        NREF  <- 4
    } else {
        pvars <- paste0(vpops, "_percent")
        leg   <- vpops
        NREF  <- length(rpops)
    }
    npvars <- length(pvars)

    # Get the reference pops
    trn_ref <- train_results[, "Refpop", drop=TRUE]
    pops    <- sort(unique(trn_ref))
    npop    <- length(pops)
    tst     <- obj$table
    tst_ref <- tst[, "Refpop", drop=TRUE]
    tst     <- as.matrix(tst[, pvars, drop=FALSE])

    if (showTrainData) trn <- as.matrix(train_results[, pvars, drop=FALSE])

    # Get the x-axis values to plot, max 6 lines per ref pop
    plotx <- seq_len(6)
    dashx <- NULL
    a     <- 0
    b     <- 7
    labx  <- 7/2
    for (i in 2:NREF) {
        m     <- max(plotx)
        add   <- (m+2):(m+7)
        plotx <- c(plotx, add)  
        a     <- a + 7
        dashx <- c(dashx, a)
        b     <- b + 14
        labx  <- c(labx, b/2)
    }
    xlim <- c(1, max(plotx))
    
    main <- "Vertex Population Ancestry Percents"
    plot(seq_len(npop), rep(0, npop), type="n", ylim=ylim, ylab="Percent",
        xlim=xlim, xlab="", xaxt="n", main=main)
    plotRefPerc_loop(pops, pvars, showTrainData, trn_ref, trn, 
                    plotx, jeps, clrs, tst_ref, tst) 
    labs <- gsub("hpgp", "", pops, fixed=TRUE)
    axis(1, at=labx, labels=labs, las=2, cex.axis=0.65, font=2)

    # Add dashed lines
    plotRefPerc_addLines(dashx)

    # Get legend
    plotRefPerc_addLegend(showTrainData, clrs, legend.pos, leg=leg)

    NULL
}

plotRefPerc_addLines <- function(dashx) {

    for (i in seq_len(length(dashx))) {
        x <- dashx[i]
        segments(x0=x, y0=0, x1=x, y1=100, lty=2, col="grey")
    }
    NULL
}

plotRefPerc_points <- function(ref, pop, dat, pv, plotxind, jeps, color, pch) {

    tmp <- ref %in% pop
    if (any(tmp)) {
        y   <- dat[tmp, pv, drop=TRUE] 
        ny  <- length(y)
        x   <- rep(plotxind, ny)
        if (jeps) x <- x + runif(ny, min=-jeps, max=jeps)
        points(x, y, col=color, pch=pch, cex=0.75) 
    }
    NULL
}

plotRefPerc_loop <- function(pops, pvars, showTrainData, trn_ref, trn, 
                            plotx, jeps, clrs, tst_ref, tst) {

    npop    <- length(pops)
    npvars  <- length(pvars)
    ind     <- 1
    pch.trn <- 3
    pch.tst <- 20
    for (i in seq_len(npop)) {
        pop <- pops[i]

        for (j in seq_len(npvars)) {
            pv <- pvars[j]

            # Plot train first
            if (showTrainData) plotRefPerc_points(trn_ref, pop, trn, pv, 
                                         plotx[ind], jeps, clrs[j], pch.trn)
            ind <- ind + 1

            # Plot test
            plotRefPerc_points(tst_ref,pop,tst,pv,plotx[ind],jeps,clrs[j], pch.tst)
            ind <- ind + 1
        }
    }
    NULL
}

plotRefPerc_addLegend <- function(showTrainData, clrs, legend.pos,
    leg=c("Africa", "Asia", "Europe")) {

    pch <- rep(20, 3)
    if (showTrainData) {
        leg  <- c(leg, "User Data", "Training Data")
        clrs <- c(clrs, "black", "black")
        pch  <- c(pch, 20, 3)
    }
    legend(legend.pos, leg, col=clrs, pch=pch, horiz=TRUE, bty="n",
            cex=0.75, text.font=2)
    NULL
}
