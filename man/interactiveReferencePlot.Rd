
\name{interactiveReferencePlot}
\alias{interactiveReferencePlot}

\title{Interactive plot of the reference data}
\usage{interactiveReferencePlot() 
}

\description{Create an interactive plot of the reference data}

%\arguments{}
\details{This plot will all show the results of all samples in the reference data.
    Hovering over a point in the plot will display three lines of information. 
    Line 1 contains the type (i.e., the source country) and id of that sample.
    Line 2 contains the sample's assigned reference population, next
    nearest reference population, and separation
    to the next nearest reference population defined as 
    \code{100*abs(d1 - d2)/d1}, where \code{d1} and \code{d2} are the
    genetic distances to the sample's assigned reference population and
    next nearest reference population respectively.
    Line 3 contains the percent African, European and Asian ancestry for that sample.
    The legend shows the abbreviated names of the source countries for all samples, 
    and clicking a
    country will add or remove those samples from the plot.
}

\value{NULL}

\seealso{ \code{\link{grafGenPlot}} }

\examples{
    if (interactive()) {
        interactiveReferencePlot()
    }
}
