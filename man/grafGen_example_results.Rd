
\name{grafGen_example_results}
\alias{grafGen_example_results}

\title{
    Example results
}
\description{
    The returned object from \code{\link{grafGen}} in the
    analysis of a subset of the reference data.
}
\details{ An object of class "grafpop" containing the 
    \code{\link{grafGen}} results for a subset  
    of 206 genomes and 35528 SNPs in the reference data. 
    This subset of the reference data is included in the package
    (\code{/extdata/data.vcf.gz}). 
}
\value{ An object of class "grafpop".}
\examples{
    data(grafGen_example_results, package="GrafGen")
    grafGen_example_results
}
\keyword{data}
\seealso{ \code{\link{example_metadata}} }
