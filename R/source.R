grafGen <- function(genoFile, referenceData=NULL, print=1) {

    checkGenoFile(genoFile)
    tmp <- check_refdata(referenceData)
    obj <- list(print=print, chr.col=tmp, prophage=0)
    ret <- gp_main(genoFile, referenceData, obj) 

    ret
}

grafGen_prophage <- function(genoFile, referenceData=NULL, print=1) {

    checkGenoFile(genoFile)
    tmp <- check_refdata(referenceData)
    obj <- list(print=print, chr.col=tmp, prophage=1)
    ret <- gp_main(genoFile, referenceData, obj) 

    ret
}

gp_main <- function(genoFile, ancSnpData, obj) {

    DEBUG       <- 0
    obj$min.snp <- 100
    chr.col     <- obj$chr.col
    prophage    <- obj$prophage
    user.ref    <- !is.null(ancSnpData)

    # Get ref data for prophages if not passed in
    if (prophage && !length(ancSnpData)) {
        ancSnpData <- gp_getDefAncData(prophage=1)
    } 

    # Check working directory
    check_out.dir(getwd(), name="getwd()")

    outfile    <- gp_getTempFile()
    ancSnpFile <- gp_getTempFile(prefix="anc_")
    ancFlag    <- length(ancSnpData)

    # Set up and write out ancestry snps
    tmp        <- gp_setUpAncSnpFile(ancSnpData, obj, ancSnpFile)
    ancSnpData <- tmp$data
    reference  <- tmp[["reference", exact=TRUE]]
    hpflag     <- !length(reference)
    tmp        <- NULL
    gc()
    nancsnp    <- nrow(ancSnpData)
    if (hpflag) {
        nrefpop <- 9
    } else {
        nrefpop <- length(reference$ref.pop)
    }
    if (nrefpop < 3) stop("with nrefpop")

    # If data has a chr column, then get the map of chrs to ints
    if (chr.col) {
      chr.map <- as.character(unique(ancSnpData[, 1, drop=TRUE])) # should already be char
      map.len <- length(chr.map)
    } else {
      chr.map <- ""
      map.len <- 0
    }

    # Save vertex and ref pop names. chr column is added to ancSpData
    vpops <- getVertexPopNames(colnames(ancSnpData), TRUE, ref=reference)
    rpops <- getRefPopNames(colnames(ancSnpData), TRUE, ref=reference)

    ancSnpData <- NULL
    gc()

    type <- gp_getFileType(genoFile) 
    if (type == 0) {
        # Remove ".bed"
        genoFile <- substr(genoFile, 1, nchar(genoFile)-4)
    } else {
        # VCF file, check for "CHROM column
        check_valid_vcf(genoFile) 
    }
    cargs   <- c(genoFile, outfile, ancSnpFile)
    iargs   <- c(type, 100, DEBUG, obj$print, nancsnp, nrefpop, !chr.col, map.len, obj$prophage)
    tmp     <- try(.C("C_main", cargs, as.integer(iargs), chr.map, 
                    PACKAGE="GrafGen"), silent=TRUE)
    file.remove(ancSnpFile) 
    err <- inherits(tmp, "try-error")
    if (err) {
        if (file.exists(outfile)) file.remove(outfile)  
        stop("calling C_main")
    }
    ret             <- readGrafOut(outfile, reference)
    tmp             <- ret$table
    tmp[, "Refpop"] <- rules.mindist(tmp, reference=reference)
    tmp             <- gp_nearest_neighbor(tmp, reference=reference)
    ret$table       <- tmp

    if (file.exists(outfile)) file.remove(outfile)  

    # Re-order columns
    ret <- gp_setReturn(ret, reference=reference)

    # Save other important info
    if (prophage && !user.ref) hpflag <- TRUE
    lst <- list(hpflag=hpflag, prophage=prophage, chr.col=chr.col, 
                vertex.pops=vpops, ref.pops=rpops, user.ref=user.ref)
    ret$objects <- lst
  
    class(ret) <- "grafpop"

    ret
}

gp_setReturn <- function(obj, reference=NULL) {

    tmp <- gp_getOrder(reference=reference)
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

    # For non-Hpylori, add reference ancestry snps
    if (length(reference)) obj$reference <- reference
    obj
}

gp_getReturnRefData <- function(obj) {

  ret <- NULL
  if (is.list(obj)) ret <- obj[["reference", exact=TRUE]]
  ret
}

gp_getOrder <- function(reference=NULL) {

    n <- length(reference)
    if (!n) {
        vertex  <- c("Africa", "Europe", "Asia")
        percent <- c("F_percent", "E_percent", "A_percent")
        refpop  <- c("hpgpAfrica", "hpgpAfrica-distant", "hpgpAfroamerica", 
        "hpgpEuroamerica", "hpgpMediterranea", "hpgpEurope",
        "hpgpEurasia", "hpgpAsia", "hpgpAklavik86-like")
    } else {
        vertex  <- reference$vertex.pop
        percent <- paste0(vertex, "_percent")
        refpop  <- reference$ref.pop
    }
    list(vertex=vertex, percent=percent, refpop=refpop)
}

gp_getRefPops <- function(ancSnpCols=NULL) {

    if (!length(ancSnpCols)) {
        ret <- c("hpgpAfrica", "hpgpAfrica-distant", "hpgpAfroamerica",
            "hpgpAklavik86-like", "hpgpAsia", "hpgpEurasia", 
            "hpgpEuroamerica", "hpgpEurope", "hpgpMediterranea")
    } else {
        ret <- getRefPopNames(ancSnpCols) 
    }
    ret
}

gp_nearest_neighbor <- function(x, reference=NULL) {

    nr    <- nrow(x)
    ref   <- x[, "Refpop", drop=TRUE]
    #uref  <- gp_getRefPops(ancSnpCols=ancSnpCols)
    uref  <- getRefPopNames(NULL, NULL, ref=reference)
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

gp_getDefAncData <- function(prophage=0) {

    HpyloriData           <- NULL
    Hpylori_prophage_data <- NULL

    dir <- system.file("data", package="GrafGen", mustWork=TRUE)
    if (!prophage) {
        f   <- file.path(dir, "HpyloriData.rda")
        load(f)
        if (is.null(HpyloriData)) stop("with HpyloriData")
        HpyloriData <- as.data.frame(HpyloriData)
        HpyloriData <- HpyloriData[, -c(1, 3)]
    } else {
        f   <- file.path(dir, "Hpylori_prophage_data.rda")
        load(f)
        if (is.null(Hpylori_prophage_data)) stop("with Hpylori_prophage_data")
        HpyloriData <- Hpylori_prophage_data
    }
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

readGrafOut.vertex <- function(f, vpops=NULL) {

    x <- scan(f, what="char", nlines=3, skip=3, sep="\n", quiet=TRUE)
    x <- gsub("\t", "", x, fixed=TRUE)
    x <- gsub("#", "", x, fixed=TRUE)
    x <- gsub(":", "", x, fixed=TRUE)
    x <- gsub("\\s+", " ",trimws(x))

    n   <- length(vpops)
    if (!n) {
      # Default HP
      nms        <- c("Africa", "Asia", "Europe")
      names(nms) <- c("F", "A", "E")
    } else {
      # Order written out in file is E, F, A
      nms        <- vpops
      names(nms) <- c("E", "F", "A")
    }
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

readGrafOut <- function(f, reference) {

    outcolnames <- getOutFileCols(reference)
    vpops       <- NULL
    if (length(reference)) {
      vpops <- reference$vertex.pop
    }
    vertex <- readGrafOut.vertex(f, vpops=vpops)
    x      <- read.table(f, skip=7, sep="\t", header=1, as.is=TRUE,
                comment.char="", check.names=FALSE)
    if (length(outcolnames) == ncol(x)) colnames(x) <- outcolnames
    for (i in 2:ncol(x)) x[, i] <- as.numeric(x[, i])

    list(table=x, vertex=vertex)
}

gp_readAncSnpFile <- function(flist, nrows=-1) {
    ret <- read.table(flist$file, header=1, sep=flist$sep, as.is=TRUE,
                    comment.char="", check.names=FALSE, nrows=nrows)
    ret
}

getTrainResults <- function(retobj=NULL, refobj=NULL, prophage=0) {

    if (!length(retobj)) {
        def <- TRUE
    } else {
        def <- !retobj$objects$user.ref
    }
    if (def) {
        if (!prophage) {
            grafGen_reference_dataframe <- NULL
            dir <- system.file("data", package="GrafGen", mustWork=TRUE)
            f   <- file.path(dir, "grafGen_reference_dataframe.rda")
            load(f)
            ret <- grafGen_reference_dataframe
        } else {
            prophage_results <- NULL
            dir <- system.file("data", package="GrafGen", mustWork=TRUE)
            f   <- file.path(dir, "prophage_results.rda")
            load(f)
            ret <- prophage_results$table
        }
    } else {
        if (length(refobj)) {
            ret <- refobj$table
        } else {
            ret <- retobj$table
        }
    }
    ret
}


