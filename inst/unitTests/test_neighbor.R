test_neighbor <- function() {

    dat <- GrafGen:::getTrainResults()
    tmp <- dat[, "Refpop", drop=TRUE] %in% "hpgpAfrica"
    x   <- dat[tmp, , drop=FALSE]
    x   <- GrafGen:::gp_nearest_neighbor(x)

    ref0 <- rep("hpgpAfroamerica", nrow(x))
    ref1 <- x[, "Nearest_neighbor", drop=TRUE]
    checkEquals(ref0, ref1)

}
