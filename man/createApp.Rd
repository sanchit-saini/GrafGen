
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
    To return an R Shiny app for the user's data.
}
\details{
    This R function returns an R Shiny app that can be launched by
    calling \code{\link[shiny]{runApp}}. The app allows the user 
    to view and filter the plot using up to two variables.
}

\value{A list containing an R Shiny app and data
    frames needed to run the app.}

\examples{
    library(GrafGen)
    data(grafGen_example_results, package="GrafGen")
    data(example_metadata, package="GrafGen")
    tmp <- createApp(grafGen_example_results, metadata=example_metadata, 
            id="Sample")
    reference_results <- tmp$reference_results
    user_results      <- tmp$user_results
    user_metadata     <- tmp$user_metadata
    if (interactive()) {
        shiny::runApp(tmp$app)  
    }
}

\seealso{ \code{\link{grafGen}} }

