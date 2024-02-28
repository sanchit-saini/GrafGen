test_mindist <- function() {

    x    <- GrafGen:::getTrainResults()
    ref0 <- x[, "Refpop", drop=TRUE]
    ref1 <- GrafGen:::rules.mindist(x)

    checkEquals(ref0, ref1)

}
