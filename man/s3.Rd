
\name{S3 methods}
\alias{plot.grafpop}
\alias{print.grafpop}

\title{Plot and Print}
\usage{
    \method{plot}{grafpop}(x, legend.pos="right", showRefData=TRUE,
        ...) 
    \method{print}{grafpop}(x, ...)
}

\description{Plot or print an object of class "grafpop".}

\arguments{
    \item{x}{An object of class "grafpop" returned from
            \code{\link{grafGen}}.}
    \item{legend.pos}{The position of the legend. 
            The default is "topleft".
            See \code{\link[graphics]{legend}}.}
    \item{showRefData}{TRUE or FALSE to display the 95 percent
            confidence ellipses for the reference
            data results.}
    %\item{interactive}{TRUE or FALSE to display an interactive
    %        plot.}
    \item{...}{Additional arguments.}
}
\details{Printing an object of class "grafpop" will display
        the frequency counts of the predicted
        reference populations.      
}
\value{NULL}
\seealso{ \code{\link{grafGen}}, \code{\link{grafGenPlot}} }

\examples{
    data(grafGen_example_results, package="GrafGen")
    obj <- grafGen_example_results
    print(obj)
    plot(obj)
}
