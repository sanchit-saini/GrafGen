
\name{grafGen_reference_results}
\alias{grafGen_reference_results}

\title{
    Reference data results
}
\description{
    The returned object from \code{\link{grafGen}} in the 
    analysis of the reference data.
}
\details{ An object of class "grafpop" containing the 
    \code{\link{grafGen}} results for each 
    of the 1011 genomes in the reference data. 
    The full set of reference data can be found at
    https://github.com/wheelerb/GrafGen/tree/reference/data .
}
\value{ An object of class "grafpop".}
\examples{
    data(grafGen_reference_results, package="GrafGen")
    grafGen_reference_results
}
\keyword{data}
\seealso{ \code{\link{grafGen}} }
