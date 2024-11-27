gp_setUpAncSnpFile <- function(x, op, out.file) {

    chr.col  <- op$chr.col
    prophage <- op$prophage
    flist    <- gp_getAncFileCols(x, chr.col)
    if (!length(x)) {
        # Assume HP analysis
        if (chr.col) stop("INTERNAL CODING ERROR")
        x         <- gp_getDefAncData()
        hpflag    <- 1
        reference <- NULL
    } else {
        hpflag <- 0
    }

    if (!hpflag) {
      chr.v      <- flist[["chr", exact=TRUE]]
      pos.v      <- flist$position
      ref.v      <- flist$ref.allele
      alt.v      <- flist$alt.allele
      vertex.pop <- flist$vertex.pop
      ref.pop    <- flist$ref.pop
      nr0        <- nrow(x)

      # Get correct order
      vv          <- c(chr.v, pos.v, ref.v, alt.v, vertex.pop, ref.pop)
      x           <- x[, vv, drop=FALSE]
    
      # Check for missing
      x <- checkForMissing(x, chr.v, ref.v, alt.v, pos.v, vertex.pop, ref.pop)

      # Positions must be integers
      x <- checkPositions(x, pos.v) 
      if (nrow(x) < op$min.snp) stop("ERROR: too few SNPs for analysis")

      # Allele freqs cannot be too large or too small
      x <- checkAlleleFreqs(x, op, vertex.pop, ref.pop)

      reference <- list(refcols=colnames(x), chr.col=chr.col,
          vertex.pop=vertex.pop, ref.pop=ref.pop)
    }

    # Add a chr column if needed
    new <- "...chr..." 
    if (!chr.col) {
      cx       <- colnames(x)
      x[, new] <- 0
      cx       <- c(new, cx)
      x        <- x[, cx, drop=FALSE]
    }
    cx <- colnames(x)
    cx[1] <- new
    colnames(x) <- cx

    # Get correct order
    if (hpflag) {
        tmp <- c(new, flist$position, flist$ref.allele, flist$alt.allele, 
                flist$vertex.pop, flist$ref.pop)
        x   <- x[, tmp, drop=FALSE] 
    }

    # Format positions
    pos.v      <- flist$position
    x[, pos.v] <- trimws(formatC(as.numeric(x[, pos.v, drop=TRUE]), digits=30))

    # Save file
    if (length(out.file)) {
        write.table(x, file=out.file, sep="\t", quote=FALSE, 
                row.names=FALSE, col.names=TRUE)
    }
    list(data=x, reference=reference)
}

gp_getDefAncFileCols <- function(prophage=0) {

    if (!prophage) {
        # order is E, F, A
        vpop <- c("vt_European", "vt_African", "vt_Asian")

        # order is NEur, Afr, EAsia, AfrAm, EurAm, SEur, Eurasia, Akl, ZAF
        rpop <- c("rf_Europe", "rf_Africa", "rf_Asia", "rf_Afroamerica",
              "rf_Euroamerica", "rf_Mediterranea",  "rf_Eurasia",
              "rf_Aklavik86.like", "rf_Africa.distant")
    } else {
        vpop <- c("v_Africa1", "v_EastAsia", "v_SWEurope")
        rpop <- c("Africa1", "EastAsia", "SWEurope", "NEurope")
    }

    list(position="pos", ref.allele="REF", alt.allele="ALT", vertex.pop=vpop, ref.pop=rpop)
} 

gp_getAncFileCols <- function(x, chr.col, prophage=0) {

    if (!length(x)) {
        ret <- gp_getDefAncFileCols(prophage=prophage)
    } else {
        # Assume x has order <chr>, pos, refAllele, altAllele, vpopAFs, refPopAFs
        cx   <- colnames(x)
        vpop <- getVertexPopNames(cx, chr.col)
        rpop <- getRefPopNames(cx, chr.col)
        if (!chr.col) {
            ret <- list(position=cx[1], ref.allele=cx[2], alt.allele=cx[3], 
                vertex.pop=vpop, ref.pop=rpop)
        } else {
            ret <- list(chr=cx[1], position=cx[2], ref.allele=cx[3], 
                alt.allele=cx[4], vertex.pop=vpop, ref.pop=rpop)
        }
    }
    ret
} 

checkForMissing <- function(x, chr.v, ref.v, alt.v, pos.v, vertex.pop, ref.pop, prophage=0) {

    nr0 <- nrow(x)

    if (length(chr.v)) {
      x[, chr.v] <- trimws(x[, chr.v, drop=TRUE])
      vec        <- x[, chr.v, drop=TRUE]
      tmp        <- is.na(vec) | (nchar(vec) < 1)
      tmp[is.na(tmp)] <- TRUE 
      if (any(tmp)) x <- x[!tmp, , drop=FALSE]
      nr <- nrow(x)
      if (nr < nr0) {
          message(paste0(nr0-nr, " rows removed due to invalid chromosomes"))
      }
      nr0 <- nr
    }

    x[, ref.v] <- toupper(trimws(x[, ref.v, drop=TRUE]))  
    x[, alt.v] <- toupper(trimws(x[, alt.v, drop=TRUE])) 
    # Single char for ref allele and alt for non-prophage
    tmp <- (nchar(x[, ref.v, drop=TRUE]) == 1) 
    if (!prophage) tmp <- tmp & (nchar(x[, alt.v, drop=TRUE]) == 1)
    if (!all(tmp)) x <- x[tmp, , drop=FALSE]
    nr <- nrow(x)
    if (nr < nr0) {
        message(paste0(nr0-nr, " rows removed due to invalid alleles"))
    }
    nr0 <- nr
    tmp <- rep(TRUE, nr)
    vv  <- c(pos.v, vertex.pop, ref.pop)
    for (v in vv) {
        x[, v] <- as.numeric(x[, v, drop=TRUE])
        tmp    <- tmp & is.finite(x[, v, drop=TRUE])
    }
    if (!all(tmp)) x <- x[tmp, , drop=FALSE]
    nr <- nrow(x)
    if (nr < nr0) {
        message(paste0(nr0-nr, " rows removed due to non-finite data"))
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

    #min.af <- op[["min.af", exact=TRUE]]
    #if (!length(min.af)) return(x)
    #if (min.af <= 0) return(x)

    #max.af <- 1 - min.af
    vv     <- c(vertex.pop, ref.pop)
    for (v in vv) {
        vec <- x[, v, drop=TRUE]
        tmp <- (vec <= 0) | (vec >= 1)
        if (any(tmp)) {
            stop("Allele frequencies must be in (0, 1)")
        } 
        #tmp <- (vec < min.af)
        #if (any(tmp)) vec[tmp] <- min.af
        #tmp <- (vec > max.af)
        #if (any(tmp)) vec[tmp] <- max.af
        #x[, v] <- vec
    }
    x
}

getOutFileDefCols <- function() {

    # Copied from util.h
    REFPOP0 <- "hpgpEurope"
    REFPOP1 <- "hpgpAfrica"
    REFPOP2 <- "hpgpAsia"
    REFPOP3 <- "hpgpAfroamerica"
    REFPOP4 <- "hpgpEuroamerica"
    REFPOP5 <- "hpgpMediterranea"
    REFPOP6 <- "hpgpEurasia"
    REFPOP7 <- "hpgpAklavik86-like"
    REFPOP8 <- "hpgpAfrica-distant"

    ret <- c("Sample", "N_SNPs", "GD1_x", "GD2_y", "GD3_z",
             "E_percent", "F_percent", "A_percent",
             REFPOP0, REFPOP1, REFPOP2, REFPOP3, 
             REFPOP4, REFPOP5, REFPOP6, REFPOP7, REFPOP8)
    ret
}

getVertexPopNames <- function(cx, chr.col, ref=NULL) {
    if (length(ref)) {
        ret <- ref$vertex.pop
    } else if (length(cx)) {
        if (!chr.col) { 
            ret <- cx[4:6]
        } else {
            ret <- cx[5:7]
        }
    } else {
        ret <- c("Europe", "Africa", "Asia")
    }
    ret
}
getRefPopNames <- function(cx, chr.col, ref=NULL) {

    if (length(ref)) {
        ret <- ref$ref.pop
    } else {
        n <- length(cx)
        if (n) {
            if (!chr.col) { 
                ret <- cx[7:n]
            } else {
                ret <- cx[8:n]
            }
        } else {
            ret <- gp_getRefPops()   
        }
    }
    ret
}

getOutFileCols <- function(reference) {

    # cx col names of ancSnpData or NULL for default cols
    n   <- length(reference)
    if (n) {
        ret <- c("Sample", "N_SNPs", "GD1_x", "GD2_y", "GD3_z",
                 paste0(reference$vertex.pop, "_percent"),
                 reference$ref.pop)
    } else {
        ret <- getOutFileDefCols()
    }
    ret
}

