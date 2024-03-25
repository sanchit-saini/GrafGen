grafGen <- function(genoFile, print=1) {

    checkGenoFile(genoFile)
    ret <- gp_main(genoFile, print) 

    ret
}

gp_main <- function(genoFile, print) {

    DEBUG <- 0

    # Check working directory
    check_out.dir(getwd(), name="getwd()")

    outfile    <- gp_getTempFile()
    ancSnpFile <- gp_getTempFile(prefix="anc_")

    # Save data in text file
    write.table(gp_getDefAncData(), file=ancSnpFile, quote=FALSE,
                sep="\t", col.names=TRUE, row.names=FALSE)

    type <- gp_getFileType(genoFile) 
    if (type == 0) {
        # Remove ".bed"
        genoFile <- substr(genoFile, 1, nchar(genoFile)-4)
    }
    cargs   <- c(genoFile, outfile, ancSnpFile)
    iargs   <- c(type, 100, DEBUG, print)
    tmp <- try(.C("C_main", cargs, as.integer(iargs), 
                PACKAGE="GrafGen"), silent=TRUE)
    file.remove(ancSnpFile) 
    err <- inherits(tmp, "try-error")
    if (err) {
        if (file.exists(outfile)) file.remove(outfile)  
        stop("calling C_main")
    }
    ret             <- readGrafOut(outfile)
    tmp             <- ret$table
    tmp[, "Refpop"] <- rules.mindist(tmp)
    tmp             <- gp_nearest_neighbor(tmp)
    ret$table       <- tmp

    if (file.exists(outfile)) file.remove(outfile)  

    # Re-order columns
    ret <- gp_setReturn(ret)

    class(ret) <- "grafpop"

    ret
}

gp_setReturn <- function(obj) {

    tmp <- gp_getOrder()
    vt  <- tmp$vertex
    per <- tmp$percent
    ref <- tmp$refpop 

    vertex     <- obj$vertex
    vertex     <- vertex[vt]
    obj$vertex <- vertex
    tab        <- obj$table
    ord        <- c("Sample", "N_SNPs", "GD1_x", "GD2_y", "GD3_z",
        per, ref, "Refpop", "Nearest_neighbor", "Separation_percent")
    obj$table  <- tab[, ord, drop=FALSE]
    obj
}

gp_getOrder <- function() {

    vertex  <- c("Africa", "Europe", "Asia")
    percent <- c("F_percent", "E_percent", "A_percent")
    refpop  <- c("hpgpAfrica", "hpgpAfrica-distant", "hpgpAfroamerica", 
        "hpgpEuroamerica", "hpgpMediterranea", "hpgpEurope",
        "hpgpEurasia", "hpgpAsia", "hpgpAklavik86-like")
    list(vertex=vertex, percent=percent, refpop=refpop)
}

gp_getRefPops <- function() {

    ret <- c("hpgpAfrica", "hpgpAfrica-distant", "hpgpAfroamerica",
            "hpgpAklavik86-like", "hpgpAsia", "hpgpEurasia", 
            "hpgpEuroamerica", "hpgpEurope", "hpgpMediterranea")
    ret
}

gp_nearest_neighbor <- function(x) {

    nr    <- nrow(x)
    ref   <- x[, "Refpop", drop=TRUE]
    uref  <- gp_getRefPops()
    mat   <- as.matrix(x[, uref, drop=FALSE]) 
    nref  <- length(uref)
    mind  <- rep(NA, nr)
    for (i in seq_len(nref)) {
        tmp         <- ref %in% uref[i]
        mind[tmp]   <- mat[tmp, i]
        mat[tmp, i] <- Inf
    }
    mind2 <- mat[, 1, drop=TRUE]
    ref2  <- rep(uref[1], nr)
    for (i in 2:nref) {
        vec <- mat[, i, drop=TRUE]
        tmp <- vec < mind2
        tmp[is.na(tmp)] <- FALSE
        if (any(tmp)) {
            mind2[tmp] <- vec[tmp]
            ref2[tmp]  <- uref[i]
        }
    }
    sep       <- abs(mind - mind2)/mind
    sep       <- sep*100
    sep       <- round(sep, digits=2)
    new1      <- "Nearest_neighbor"
    new2      <- "Separation_percent"
    x[, new1] <- ref2
    x[, new2] <- sep

    x
}

gp_getDefAncData <- function() {

    HpyloriData <- NULL
    dir <- system.file("data", package="GrafGen", mustWork=TRUE)
    f   <- file.path(dir, "HpyloriData.rda")
    load(f)
    if (is.null(HpyloriData)) stop("with HpyloriData")
    HpyloriData <- as.data.frame(HpyloriData)
    HpyloriData <- HpyloriData[, -c(1, 3)]
    HpyloriData
}

gp_getTempFile <- function(prefix="tmp_") {

    dir <- getwd()
    vec <- c(letters, as.character(0:9))
    vec <- sample(vec, 12)
    ret <- paste0(dir, "/", prefix, paste0(vec, collapse=""))
    ret 
}

gp_getFileType <- function(genoFile) {

    #define BED_FILE 0
    #define VCF_FILE 1
    f   <- tolower(basename(genoFile))
    f   <- gsub(".gz", "", f, fixed=TRUE)
    nc  <- nchar(f)
    ext <- substr(f, nc-3, nc)
    if (ext == ".bed") {
        ret <- 0
    } else if (ext == ".vcf") {
        ret <- 1
    } else {
        stop("with file extension")
    }
    ret 
}

readGrafOut.vertex <- function(f) {

    x <- scan(f, what="char", nlines=3, skip=3, sep="\n", quiet=TRUE)
    x <- gsub("\t", "", x, fixed=TRUE)
    x <- gsub("#", "", x, fixed=TRUE)
    x <- gsub(":", "", x, fixed=TRUE)
    x <- gsub("\\s+", " ",trimws(x))

    nms <- c("Africa", "Asia", "Europe")
    names(nms) <- c("F", "A", "E")
    ret <- list()
    for (i in seq_len(length(x))) {
        vec <- trimws(unlist(strsplit(x[i], " ", fixed=TRUE)))
        anc <- vec[1]
        x1  <- as.numeric(vec[2])
        y1  <- as.numeric(vec[3])
        nm  <- nms[anc]
        ret[[nm]] <- list(x=x1, y=y1)
    }
    ret
}

readGrafOut <- function(f) {

    vertex <- readGrafOut.vertex(f)
    x      <- read.table(f, skip=7, sep="\t", header=1, as.is=TRUE,
                comment.char="", check.names=FALSE)
    for (i in 2:ncol(x)) x[, i] <- as.numeric(x[, i])

    list(table=x, vertex=vertex)
}

gp_readAncSnpFile <- function(flist, nrows=-1) {
    ret <- read.table(flist$file, header=1, sep=flist$sep, as.is=TRUE,
                    comment.char="", check.names=FALSE, nrows=nrows)
    ret
}

getTrainResults <- function() {

    grafGen_reference_dataframe <- NULL
    dir <- system.file("data", package="GrafGen", mustWork=TRUE)
    f   <- file.path(dir, "grafGen_reference_dataframe.rda")
    load(f)
    grafGen_reference_dataframe
}
