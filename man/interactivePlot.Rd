
\name{interactivePlot}
\alias{interactivePlot}

\title{Interactive plot of user's data}
\usage{
    interactivePlot(obj, metadata=NULL, id=NULL, type=NULL, group=NULL) 
}

\description{Create an interactive plot of user data}

\arguments{
    \item{obj}{Return object from \code{\link{grafGen}}.}
    \item{metadata}{NULL or data frame containing meta data
        for the plot. This data frame must contain an id 
        variable.}
    \item{id}{Name of the id column in \code{metadata}.
        If NULL, then the first column will be used.}
    \item{type}{Name of the type variable in \code{metadata}, which is
        rendered as different colors. This variable should be
        categorical.
        If NULL, then NA values will appear when hovering over points.}
    \item{group}{Name of the group variable in \code{metadata}, which is
        rendered as different symbols.
        This variable should be categorical. If NULL, then
        \code{group} will be set to \code{type}.}
}

\details{This plot will all show the results of all samples in the user's data.
    Hovering over a point in the plot will display three lines of information. 
    Line 1 contains the group, type and id of that sample.
    Line 2 contains the sample's assigned reference population, next
    nearest reference population, and separation
    to the next nearest reference population defined as 
    \code{100*abs(d1 - d2)/d1}, where \code{d1} and \code{d2} are the
    genetic distances to the sample's assigned reference population and
    next nearest reference population respectively.
    Line 3 contains the percent African, European and Asian ancestry for that sample.
    The legend shows the types for all samples, and clicking a
    type will add or remove those samples from the plot.

    Note that printing the returned object from \code{\link{grafGen}} with the command
    \code{print(obj)} will display the frequency counts for each reference population.
}

\value{NULL}

\seealso{ \code{\link{grafGenPlot}} }

\examples{
    if (interactive()) {
        data(grafGen_example_results, package="GrafGen")
        data(example_metadata, package="GrafGen")
        interactivePlot(grafGen_example_results, metadata=example_metadata, 
            id="Sample", type="Country")
    }
}
