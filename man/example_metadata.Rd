
\name{example_metadata}
\alias{example_metadata}

\title{
    Metadata for examples
}
\description{
    A data frame containing metadata
    used in examples.
}
\details{ The data frame contains the sample id,
    type (i.e. source country),
    and country abbreviation for the 206 genomes in 
    \code{\link{grafGen_example_results}}. 
}
\value{ A data frame}
\examples{

    data(example_metadata, package="GrafGen")

    # Display a few rows
    example_metadata[seq_len(5), ]
}
\keyword{data}
\seealso{ \code{\link{grafGen_example_results}} }
