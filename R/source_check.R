checkGenoFile <- function(x, nm="genoFile", ext=NULL) {

    if (!isFile(x)) {
        msg <- paste0("ERROR: ", nm, " must be a valid file")
        stop(msg)  
    }
    if (!length(ext)) ext <- c(".bed", ".vcf", ".vcf.gz")
    ok  <- checkFileExtensions(x, checkForExt=ext, case=0) 
    if (!ok) {
        str <- paste0("'", ext, "'")
        str <- paste0(str, collapse=", ")
        msg <- paste0("ERROR: ", nm, " must have an extension ", str)
        stop(msg)  
    }

    # If bed file, check for bim and fam
    if (isBedFile(x)) {
        len   <- nchar(x)
        bfile <- substr(x, 1, len-4)
        bim   <- paste0(bfile, ".bim")
        if (!file.exists(bim)) stop(paste0("ERROR: file ", bim, " not found"))
        fam   <- paste0(bfile, ".fam")
        if (!file.exists(fam)) stop(paste0("ERROR: file ", fam, " not found"))
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

checkAncFileList <- function(x, nm="ancSnpFile") {

    req <- c("file", "position", "ref.allele", "alt.allele",
            "vertex.pop", "ref.pop")
    checkRequiredListNames(x, req, nm)

    # Check file
    f <- x$file
    if (!isFile(f) && !is.data.frame(f)) {
        msg <- paste0("ERROR: ", nm, "$file must be a file or data frame")
        stop(msg) 
    }
    if (isFile(f)) {
        f  <- try(gp_readAncSnpFile(x, nrows=2))
        if ("try-error" %in% class(f)) {
            msg <- paste0("ERROR reading ", x$file)
            stop(msg)
        } 
    }

    # See that is has required columns 
    cnames <- colnames(f)
    nms    <- req[-1]
    for (y in nms) {
        var <- x[[y, exact=TRUE]]
        if (!length(var)) {
            msg <- paste0("ERROR: ", nm, "$", y, " must be specified")
            stop(msg)   
        }  
        if (!(var %in% cnames)) {
            msg <- paste0("ERROR: ", nm, "$", y, "=", var, " not found in data")
            stop(msg)   
       }
    }
    NULL
}

check_options <- function(x, nm="options") {

    valid <- c("DEBUG", "min.snp", "min.af", "max.af")
    def   <- list(0, 100, 1e-4, 1-1e-4)

    if (!length(x)) x <- list()
    if (!isList(x)) stop(paste0("ERROR: ", nm, " must be a list"))
    checkOptionListNames(x, valid, nm)
    x <- default.list(x, def)
    x
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

checkFiles <- function(filevec, name="filevec") {

    if (!length(filevec)) stop("ERROR 0")
    if (!is.character(filevec)) stop("ERROR 1")
    filevec <- trimws(filevec)
    tmp     <- !file.exists(filevec)
    if (any(tmp)) {
        stop("ERROR: file does not exist")
    }
    filevec
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

check_out.file <- function(x, nm="out.file", valid.ext=NULL, delIfExists=1) {

    if (!isString(x)) stop("ERROR 0")
    x <- trimws(x)
    if (length(valid.ext)) {
        tmp <- checkFileExtensions(tolower(x), checkForExt=valid.ext)
        if (!tmp) {
            msg <- "ERROR: not valid file extention"
            stop(msg) 
        } 
    }
    dir <- dirname(x)
    check_out.dir(dir, name=nm) 
    if (delIfExists && file.exists(x)) file.remove(x)

    x 

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

check.dir <- function(x, name) {

    if (!isString(x)) stop("ERROR 0")
    x <- trimws(x)
    if (!dir.exists(x)) stop("ERROR 1")
    x <- checkForSep(x)
    x
}

check_out.dir <- function(x, name="out.dir") {

    if (!isString(x)) stop(paste0("ERROR: with ", name))
    x <- trimws(x)
    if (!dir.exists(x)) stop(paste0("ERROR: ", name, " does not exist"))
    if (file.access(x, mode=2)) {
        stop(paste0("ERROR: ", name, " needs write permission"))
    }
    x <- gsub("\\", "/", x, fixed=TRUE)
    x <- checkForSep(x)
    x
}

isList <- function(x) {

    ret <- is.list(x) && ("list" %in% class(x))
    ret
}

checkOptionListNames <- function(op, valid, name) {
    if (!length(op)) return(NULL)

    # Names cannot be ""
    nms <- trimws(names(op))
    tmp <- nchar(nms) < 1
    if (any(tmp)) {
        stop("ERROR: list must have names") 
    }
    if (length(valid)) {
        tmp <- !(nms %in% valid)
        if (any(tmp)) {
            nms <- paste0("'", nms[tmp], "'")
            err <- paste0(nms, collapse=",")
            stop(paste0("ERROR: the option(s) ", err, " are not valid"))
        }
    }  

    NULL

} # END: checkOptionListNames

checkRequiredListNames <- function(x, req, name) {

    if (!length(x)) stop(paste0("ERROR: ", name, " has length 0"))
    if (!is.list(x)) stop(paste0("ERROR: ", name, " must be a list"))

    tmp  <- !(req %in% names(x))
    miss <- req[tmp]
    if (length(miss)) {
        miss <- paste0("'", miss, "'")
        tmp  <- paste0(miss, collapse=", ")
        msg  <- paste0("ERROR: ", tmp, " not found in ", name)
        stop(msg)  
    }
    NULL

} # END: checkRequiredListNames

check.list <- function(x, name, valid) {

    if (!isList(x)) stop(paste0("ERROR: ", name, " must be a list"))
    checkOptionListNames(x, valid, name) 
    NULL 

} # END: check.list

check_numVec <- function(x, nm, finite=1, int=1) {

    n <- length(x)
    if (!n) stop(paste0("ERROR: ", nm, " has length 0"))
    if (!is.numeric(x)) stop(paste0("ERROR: ", nm, " must be numeric"))
    if (!is.vector(x)) stop(paste0("ERROR: ",nm," must be a numeric vector"))
    if (finite) {
        tmp <- !is.finite(x)
        if (any(tmp)) stop(paste0("ERROR: ",nm," contains non-finite values"))
    }
    if (int) {
        tmp <- x != floor(x) 
        if (any(tmp)) stop(paste0("ERROR: ",nm," contains non-integer values"))
    }

    NULL
}

check_grafpop <- function(x, nm="obj") {

    if (!("grafpop" %in% class(x))) {
        msg <- paste0("ERROR: ", nm, " must be of class grafpop")
        stop(msg)
    }
    NULL
}

check_numVec <- function(x, nm, len=NA, valid=NULL, def=NULL) {

    n <- length(x)
    if (!n) {
        if (is.na(len)) return(def)
        if (n != len) {
            msg <- paste0("ERROR ", nm, " must have length ", len)
            stop(msg)
        }
    }
    if (!is.numeric(x)) stop(paste0("ERROR ", nm, " must be numeric"))
    if (!is.vector(x)) stop(paste0("ERROR ", nm, " must be a vector"))
    if (is.finite(len) && (n != len)) {
        stop(paste0("ERROR ", nm, " must have length ", len))
    }
    if (length(valid)) {
        tmp <- !(x %in% valid)
        if (any(tmp)) {
            msg <- paste0("ERROR ", nm, " contains invalid values")
            stop(msg)
        }
    }

    x
}

check_legendPos <- function(x, nm="legend.pos") {

    if (!length(x)) return(NULL)
    if (!isString(x)) {
        msg <- paste0("ERROR: ", nm, " must be a string")
        stop(msg)
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
        msg <- paste0("ERROR: ", name, " must be logical")
        stop(msg)
    }
    NULL
}

check_metadata <- function(x, idv, nm="metadata") {

    if (!length(x)) return(NULL)
    if (!is.data.frame(x)) stop(paste0("ERROR: ", nm, " must be a data frame"))
    if (!nrow(x)) stop(paste0("ERROR: ", nm, " contains no rows"))
    if (!ncol(x)) stop(paste0("ERROR: ", nm, " contains no columns"))
    cx <- colnames(x)
    if (!length(idv)) idv <- cx[1]
    if (!isString(idv)) stop(paste0("ERROR: id must be a column in ", nm))
    if (!(idv %in% cx)) stop(paste0("ERROR: id must be a column in ", nm))
    NULL
}

check_variable <- function(var, nm, dat) {

    if (!length(dat)) return(NULL)
    if (!length(var)) return(NULL)
    if (!isString(var)) stop(paste0("ERROR: '", nm, "' must be a column in the data"))
    cx <- colnames(dat)
    if (!(var %in% cx)) stop(paste0("ERROR: '", nm, "' must be a column in the data"))
    NULL
}

# Function to assign a default value to an element in a list
default.list <- function(inList, names, default, error=NULL,
                        checkList=NULL) {

    n1 <- length(names)
    n2 <- length(default)
    if (n1 != n2) stop("ERROR: in calling default.list")

    if (is.null(error)) {
        error <- rep(0, times=n1)
    } else if (n1 != length(error)) {
        stop("ERROR: in calling default.list")
    }

    if (!is.null(checkList)) {
        if (n1 != length(checkList)) stop("ERROR: in calling default.list")
        checkFlag <- 1
    } else {
        checkFlag <- 0
    } 

    if (is.null(inList)) inList <- list()

    listNames <- names(inList)
    for (i in seq_len(n1)) {
        if (!(names[i] %in% listNames)) {
            if (!error[i]) {
                inList[[names[i]]] <- default[[i]]
            } else {
                temp <- paste("ERROR: the name ", names[i], 
                " was not found", sep="")
                stop(temp)
            }
        }
    }

    inList

} # END: default.list

