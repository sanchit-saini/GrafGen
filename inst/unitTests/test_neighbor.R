test_neighbor <- function() {

    dir <- system.file("data", package="GrafGen", mustWork=TRUE)
    f   <- file.path(dir, "train_results.rda")
    load(f)

    tmp <- train_results[, "Refpop", drop=TRUE] %in% "hpgpAfrica"
    x   <- train_results[tmp, , drop=FALSE]
    x   <- GrafGen:::gp_nearest_neighbor(x)

    ref0 <- rep("hpgpAfroamerica", nrow(x))
    ref1 <- x[, "Nearest_neighbor", drop=TRUE]
    checkEquals(ref0, ref1)

}
