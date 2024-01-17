
\name{interactiveTestPlot}
\alias{interactiveTestPlot}

\title{Interactive plot of user test data}
\usage{
    interactiveTestPlot(obj, metadata=NULL, id=NULL, country=NULL, group=NULL) 
}

\description{Create an interactive plot of user data}

\arguments{
    \item{obj}{Return object from \code{\link{grafGen}}.}
    \item{metadata}{NULL or data frame containing meta data
        for the plot. This data frame must contain an id 
        variable.}
    \item{id}{Name of the id column in \code{metadata}.
        If NULL, then the first column will be used.}
    \item{country}{Name of the country variable in \code{metadata}.
        If NULL, then NA values will appear when hovering over points.}
    \item{group}{Name of the group variable in \code{metadata}.
        This variable should be categorical. If NULL, then
        \code{group} will be set to \code{country}.}
}

\details{This plot will all show the results of all samples in the user's data.
    Hovering over a point in the plot will display three lines of information. 
    Line 1 contains the group, country and id of that sample.
    Line 2 contains the sample's assigned reference population, next
    nearest reference population, and separation
    to the next nearest reference population defined as 
    \code{100*abs(d1 - d2)/d1}, where \code{d1} and \code{d2} are the
    genetic distances to the sample's assigned reference population and
    next nearest reference population respectively.
    Line 3 contains the percent Asian, African, and European for that sample.
    The legend shows the country abbreviations for all samples, and clicking a
    country will add or remove those samples from the plot.
}

\value{NULL}

\seealso{ \code{\link{grafGenPlot}} }

\examples{
    if (interactive()) {
        data(grafGen_results, package="GrafGen")
        data(train_results, package="GrafGen")
        interactiveTestPlot(grafGen_results, metadata=train_results, 
            id="Sample", country="Country")
    }
}
