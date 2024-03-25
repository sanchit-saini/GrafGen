checkGenoFile <- function(x, nm="genoFile", ext=NULL) {

    if (!isFile(x)) stop(nm, " must be a valid file")  
    if (!length(ext)) ext <- c(".bed", ".vcf", ".vcf.gz")
    ok  <- checkFileExtensions(x, checkForExt=ext, case=0) 
    if (!ok) {
        str <- paste0("'", ext, "'")
        str <- paste0(str, collapse=", ")
        stop(nm, " must have an extension ", str)  
    }

    # If bed file, check for bim and fam
    if (isBedFile(x)) {
        len   <- nchar(x)
        bfile <- substr(x, 1, len-4)
        bim   <- paste0(bfile, ".bim")
        if (!file.exists(bim)) stop(bim, " not found")
        fam   <- paste0(bfile, ".fam")
        if (!file.exists(fam)) stop(fam, " not found")
    }
    NULL
}

isBedFile <- function(x) {

    ret <- FALSE
    if (isFile(x)) {
        len <- nchar(x)
        ext <- tolower(substr(x, len-3, len))
        ret <- ext == ".bed"
    }
    ret
}

# Function to check that an object is a string
isString <- function(obj) {

    if ((length(obj) == 1) && is.character(obj)) {
        ret <- TRUE
    } else {
        ret <- FALSE
    }

    ret

} # END: isString

isFile <- function(x) {
    isString(x) && file.exists(x)
}

checkFileExtensions <- function(filevec, checkForExt=c(".xlsx", ".rda"), 
                                case=0) {

    tmp <- substr(checkForExt, 1, 1) != "."
    if (any(tmp)) checkForExt <- paste0(".", checkForExt[tmp])
    lens   <- unique(nchar(checkForExt))
    ok     <- rep(FALSE, length(filevec))
    flen   <- nchar(filevec)
    if (!case) {
        checkForExt <- tolower(checkForExt)
        filevec     <- tolower(filevec)
    }
    for (i in seq_len(length(lens))) {
        len <- lens[i]
        tmp <- substr(filevec, flen-len+1, flen) %in% checkForExt
        ok  <- ok | tmp
    }
    ok
}

checkForSep <- function(x) {

    sep <- .Platform$file.sep
    x   <- trimws(x)
    n   <- nchar(x)
    if (!n) {
        x <- "./"
    } else {
        if (substr(x, n, n) != sep) x <- paste0(x, sep)
    }
    x
}

check_out.dir <- function(x, name="out.dir") {

    if (!isString(x)) stop(" with ", name)
    x <- trimws(x)
    if (!dir.exists(x)) stop(name, " does not exist")
    if (file.access(x, mode=2)) {
        stop(name, " needs write permission")
    }
    x <- gsub("\\", "/", x, fixed=TRUE)
    x <- checkForSep(x)
    x
}

checkRequiredListNames <- function(x, req, name) {

    if (!length(x)) stop(name, " has length 0")
    if (!is.list(x)) stop(name, " must be a list")

    tmp  <- !(req %in% names(x))
    miss <- req[tmp]
    if (length(miss)) {
        miss <- paste0("'", miss, "'")
        tmp  <- paste0(miss, collapse=", ")
        stop(tmp, " not found in ", name)  
    }
    NULL

} # END: checkRequiredListNames

check_grafpop <- function(x, nm="obj") {

    if (!inherits(x, "grafpop")) {
        stop(nm, " must be of class grafpop")
    }
    NULL
}

check_numVec <- function(x, nm, len=NA, valid=NULL, def=NULL) {

    n <- length(x)
    if (!n) {
        if (is.na(len)) return(def)
        if (n != len) {
            stop(nm, " must have length ", len)
        }
    }
    if (!is.numeric(x)) stop(nm, " must be numeric")
    if (!is.vector(x)) stop(nm, " must be a vector")
    if (is.finite(len) && (n != len)) {
        stop(nm, " must have length ", len)
    }
    if (length(valid)) {
        tmp <- !(x %in% valid)
        if (any(tmp)) {
            stop(nm, " contains invalid values")
        }
    }

    x
}

check_legendPos <- function(x, nm="legend.pos") {

    if (!length(x)) return(NULL)
    if (!isString(x)) {
        stop(nm, " must be a string")
    }
    x
}

check.logical <- function(x, name, len=1) {

    err <- 0
    n   <- length(x)
    if (!n && !len) return(NULL)

    if (!n || (n > 1)) {
        err <- 1
    } else {
        if ((x != 0) && (x != 1)) err <- 1
    }
    if (err) {
        stop(name, " must be logical")
    }
    NULL
}

check_metadata <- function(x, idv, nm="metadata") {

    if (!length(x)) return(NULL)
    if (!is.data.frame(x)) stop(nm, " must be a data frame")
    if (!nrow(x)) stop(nm, " contains no rows")
    if (!ncol(x)) stop(nm, " contains no columns")
    cx <- colnames(x)
    if (!length(idv)) idv <- cx[1]
    if (!isString(idv)) stop("id must be a column in ", nm)
    if (!(idv %in% cx)) stop("id must be a column in ", nm)
    NULL
}

check_variable <- function(var, nm, dat) {

    if (!length(dat)) return(NULL)
    if (!length(var)) return(NULL)
    if (!isString(var)) stop("'", nm, "' must be a column in the data")
    cx <- colnames(dat)
    if (!(var %in% cx)) stop("'", nm, "' must be a column in the data")
    NULL
}
