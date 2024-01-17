
\name{extractPositions}
\alias{extractPositions}

\title{Extract SNPs by position}
\usage{
    extractPositions(vcf.file, out.file, positions=NULL) 
}
\arguments{
\item{vcf.file}{The complete path to the VCF file.}
\item{out.file}{The complete path to the output file. 
                If this file ends in ".gz", then the
                file will be compressed. If this file is 
                going to be used with \code{\link{grafGen}},
                then it must have a .vcf or .vcf.gz extension.}
\item{positions}{NULL or a vector of positions to extract.
                If NULL, then the positions of the 
                ancestry SNPs for H. pylori will be
                extracted.}
}
\description{
    To extract SNPs from a VCF file by position.
}
\details{
    This R function calls a python program stored in
    /exec/extractByLoc.py.
}

\value{NULL, the extracted genotypes are stored in \code{out.file}.}

\examples{
    dir      <- system.file("extdata", package="GrafGen", mustWork=TRUE)
    vcf.file <- file.path(dir, .Platform$file.sep, "data.vcf.gz")
    out.file <- paste0(getwd(), .Platform$file.sep, "extractData.vcf")
    \donttest{
    extractPositions(vcf.file, out.file)
    }
}

\seealso{ \code{\link{grafGen}} }

