
\name{grafGen_reference_dataframe}
\alias{grafGen_reference_dataframe}

\title{
    Reference data results for plots
}
\description{
    A data frame of the reference data results used in
    creating plots.
}
\details{ The data frame contains the results for each 
    of the 1011 genomes in the reference data used
    in training the model along with some additional 
    columns. 
}
\value{ A data frame}
\examples{

    data(grafGen_reference_dataframe, package="GrafGen")

    # Display a few rows
    grafGen_reference_dataframe[seq_len(5), ]
}
\keyword{data}
\seealso{ \code{\link{grafGen_reference_results}} }
