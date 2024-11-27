
\name{S3 methods}
\alias{plot.grafpop}
\alias{print.grafpop}

\title{Plot and Print}
\usage{
    \method{plot}{grafpop}(x, legend.pos="right", showRefData=TRUE,
        refObj=NULL, ...) 
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
    \item{refObj}{The returned object from \code{\link{grafGen}} for
        the reference data. This should only be passed in for
        non-H. pylori genomes and if \code{x} is not
        the result of the reference data.}
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
