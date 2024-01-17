
\name{HpyloriData}
\alias{HpyloriData}

\title{
    H. pylori data
}
\description{ 
    Reference allele frequencies and poistions of H. pylori training SNPs
}
\details{
    A \code{\link[GenomicRanges]{GPos}} class object containing
    the vertex and reference population 
    reference allele frequencies for the
    set of SNPs used in the analysis for H. pylori.
    The SNPs were created using Tomb26695 as the reference genome.
    The set of SNPs was selected using a MAF threshold of 0.01.
    The total sample size was from a set of 1011 H. pylori strains.
}
\examples{
    # Load data and view the first few rows
    data(HpyloriData, package="GrafGen")
    HpyloriData
}
\keyword{data}
