
\name{grafGenPlot}
\alias{grafGenPlot}

\title{Plot results}
\usage{grafGenPlot(obj, which=1, legend.pos=NULL, 
        ylim=NULL, showRefData=TRUE,
        jitter=0) 
}

\description{Plot results}

\arguments{
    \item{obj}{An object of class "grafpop" returned from
            \code{\link{grafGen}}.}
    \item{which}{A vector of integers in \code{1, 2, 3, 4, 5} to
            determine which plots are produced where:
    \tabular{ll}{
        which \tab plot \cr
        1     \tab GD1_x vs GD2_y \cr
        2     \tab GD1_x vs GD3_z \cr
        3     \tab GD2_y vs GD3_z \cr
        4     \tab Distance from each strain to its predicted
                    population \cr
        5     \tab Ancestry percents for each strain within
                    each reference population
    }
    } 
    \item{legend.pos}{The position of the legend. 
                    See \code{\link[graphics]{legend}}.}
    \item{ylim}{NULL or the limits of the y-axis. 
                    See \code{\link[graphics]{plot}}.}
    \item{showRefData}{TRUE or FALSE to display the 95 percent
                    confidence ellipses for the 
                    reference data results.}
    \item{jitter}{Numeric value for the amount of jitter to add
                for the plot \code{which = 5}.
                Values less than 0.5 work well.
                The default is 0.}
    %\item{interactive}{TRUE or FALSE to display an interactive
    %                plot.}
}
\details{The option \code{legend.pos} is only available for 
    \code{which = 1-3},
    option \code{ylim} is only available for \code{which = 4-5},
    and option \code{jitter} is only available for \code{which = 5}.}
\value{NULL}
\seealso{ \code{\link{grafGen}} }

\examples{
    data(grafGen_example_results, package="GrafGen")
    grafGenPlot(grafGen_example_results)
}
