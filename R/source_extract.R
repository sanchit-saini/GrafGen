extractPositions <- function(vcf.file, out.file, positions=NULL) {

    checkGenoFile(vcf.file, nm="vcf.file", ext=c(".vcf", ".vcf.gz"))
    out.file <- check_out.file(out.file)
    if (length(positions)) {
        check_numVec(positions, "positions")
    } else {
        positions <- getPosFromAncSnpFile()
    }
    extractByLocation(vcf.file, list(pos=positions), out.file)

    NULL
}

getPosFromAncSnpFile <- function() {

    x   <- gp_getDefAncData()
    ret <- as.numeric(x[, 1])
    ret
}

extractByLocation <- function(fileVec, locList, outfile, op=NULL) {

    n <- length(fileVec)
    if (length(locList) != n) {
        stop("ERROR: fileVec and locList must have the same length")
    }
    tmp  <- !file.exists(fileVec)
    miss <- fileVec[tmp]
    if (length(miss)) {
        stop("ERROR: the file does not exist")
    }
    op <- default.list(op, 
            c("print", "delete", "locCol", "isRange"), 
            list(0, 1, 2, 0))
    python.prog <- op[["python.prog", exact=TRUE]]
    tmpfile <- paste(outfile, "_", gp_getRandStr(), sep="")
    if (is.null(python.prog)) {
        python.prog <- system.file("exec", "extractByLoc.py", 
                    package="GrafGen", mustWork=TRUE)
    }
    if (is.null(python.prog)) stop("ERROR: python.prog is NULL")
    x <- rep("", n)
    for (i in seq_len(n)) {
        vec <- as.numeric(locList[[i]])
        tmp <- !is.finite(vec)
        if (any(tmp)) stop("ERROR: non-finite values in locList")
        vec  <- trimws(formatC(as.numeric(vec), digits=30))
        x[i] <- paste(c(fileVec[i], vec), collapse=",", sep="")
    }
    write(x, file=tmpfile, ncolumns=1)
    str <- paste("python ", python.prog, " files=", tmpfile,
                " out=", outfile, " locCol=", op$locCol, 
                " print=", op$print, " range=", op$isRange, 
                sep="")
    callOS(str)
    if (op$delete) file.remove(tmpfile)

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

# Function to call the operating system
callOS <- function(command, intern=FALSE) {

    # Determine the platform
    os      <- .Platform$OS.type
    winFlag <- (os == "windows")

    if (winFlag) {
        ret <- shell(command, intern=intern)
    } else {
        ret <- system(command, intern=intern)
    }
    ret

} # END: callOS
