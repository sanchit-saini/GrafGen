
\name{train_results}
\alias{train_results}

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
\examples{

    data(train_results, package="GrafGen")

    # Display a few rows
    train_results[seq_len(5), ]
}
\keyword{data}
\seealso{ \code{\link{plot.grafpop}} }
