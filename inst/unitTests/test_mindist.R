test_mindist <- function() {

    dir <- system.file("data", package="GrafGen", mustWork=TRUE)
    f   <- file.path(dir, "train_results.rda")
    load(f)

    ref0 <- train_results[, "Refpop", drop=TRUE]
    x    <- train_results
    ref1 <- GrafGen:::rules.mindist(x)

    checkEquals(ref0, ref1)

}
