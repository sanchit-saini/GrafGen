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
    err <- "try-error" %in% class(tmp)
    if (err) {
        if (file.exists(outfile)) file.remove(outfile)  
        stop("ERROR calling C_main")
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
    if (is.null(HpyloriData)) stop("ERROR with HpyloriData")
    HpyloriData <- as.data.frame(HpyloriData)
    HpyloriData <- HpyloriData[, -c(1, 3)]
    HpyloriData
}

gp_getRandStr <- function() {

    vec <- c(letters, as.character(0:9))
    vec <- sample(vec, 6)
    ret <- paste0(vec, collapse="")
    ret
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
        stop("ERROR")
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

gp_setUpAncSnpFile <- function(flist, op, out.file) {

    pos.v      <- flist$position
    ref.v      <- flist$ref.allele
    alt.v      <- flist$alt.allele
    vertex.pop <- flist$vertex.pop
    ref.pop    <- flist$ref.pop

    if (isFile(flist$file)) {
        x <- gp_readAncSnpFile(flist, nrows=-1)
    } else {
        x <- flist$file
    }
    nr0 <- nrow(x)

    # Get correct order
    vv <- c(pos.v, ref.v, alt.v, vertex.pop, ref.pop)
    x  <- x[, vv, drop=FALSE]

    # Check for missing
    x <- checkForMissing(x, ref.v, alt.v, pos.v, vertex.pop, ref.pop)

    # Positions must be integers
    x <- checkPositions(x, pos.v) 
    if (nrow(x) < op$min.snp) stop("ERROR: too few SNPs for analysis")

    # Allele freqs cannot be too large or too small
    x <- checkAlleleFreqs(x, op, vertex.pop, ref.pop)

    # Save file
    if (length(out.file)) {
        write.table(x, file=out.file, sep="\t", quote=FALSE, 
                row.names=FALSE, col.names=FALSE)
    }
    x
}

checkForMissing <- function(x, ref.v, alt.v, pos.v, vertex.pop, ref.pop) {

    nr0        <- nrow(x)
    x[, ref.v] <- toupper(trimws(x[, ref.v, drop=TRUE]))  
    x[, alt.v] <- toupper(trimws(x[, ref.v, drop=TRUE])) 
    tmp        <- (nchar(x[, ref.v, drop=TRUE]) > 0) &
                    (nchar(x[, alt.v, drop=TRUE]) > 0)
    if (!all(tmp)) x <- x[tmp, , drop=FALSE]
    nr <- nrow(x)
    if (nr < nr0) {
        message(paste0(nr-nr0, " rows removed due to invalid alleles"))
    }
    nr0 <- nr
    tmp <- rep(TRUE, nr)
    vv  <- c(pos.v, ref.v, alt.v, vertex.pop, ref.pop)
    for (v in vv) {
        x[, v] <- as.numeric(x[, v, drop=TRUE])
        tmp    <- tmp & is.finite(x[, v, drop=TRUE])
    }
    if (!all(tmp)) x <- x[tmp, , drop=FALSE]
    nr <- nrow(x)
    if (nr < nr0) {
        message(paste0(nr-nr0, " rows removed due to non-finite data"))
    }
    x
}

checkPositions <- function(x, pos.v) {

    vec <- x[, pos.v, drop=TRUE]
    tmp <- vec != floor(vec)
    m   <- sum(tmp)
    if (m) {
        x <- x[!tmp, , drop=FALSE]
        message(paste0(m, " rows removed due to non-integer positions"))
    }
    x
}

checkAlleleFreqs <- function(x, op, vertex.pop, ref.pop) {

    min.af <- op$min.af
    max.af <- op$max.af
    vv     <- c(vertex.pop, ref.pop)
    for (v in vv) {
        vec <- x[, v, drop=TRUE]
        tmp <- (vec < min.af)
        if (any(tmp)) vec[tmp] <- min.af
        tmp <- (vec > max.af)
        if (any(tmp)) vec[tmp] <- max.af
        x[, v] <- vec
    }
    x
}

getTrainResults <- function() {

    grafGen_reference_dataframe <- NULL
    dir <- system.file("data", package="GrafGen", mustWork=TRUE)
    f   <- file.path(dir, "grafGen_reference_dataframe.rda")
    load(f)
    grafGen_reference_dataframe
}
