
\name{createApp}
\alias{createApp}

\title{R Shiny App}
\usage{
    createApp(obj, metadata=NULL, id=NULL) 
}
\arguments{
    \item{obj}{Return object from \code{\link{grafGen}}.}
    \item{metadata}{NULL or data frame containing meta data
        for the plot. This data frame must contain an id 
        variable.}
    \item{id}{Name of the id column in \code{metadata}.
        If NULL, then the first column will be used.}
}
\description{
    To return an R Shiny app for the reference data.
}
\details{
    This R function returns an R Shiny app that can be launched by
    calling \code{\link[shiny]{runApp}}.
}

\value{A list containing an R Shiny app and data
    frames needed to run the app.}

\examples{
    library(GrafGen)
    data(grafGen_results, package="GrafGen")
    tmp <- createApp(grafGen_results)
    train_results <- tmp$train_results
    test_results  <- tmp$test_results
    test_metadata <- tmp$test_metadata
    if (interactive()) {
        shiny::runApp(tmp$app)  
    }
}

\seealso{ \code{\link{grafGen}} }

