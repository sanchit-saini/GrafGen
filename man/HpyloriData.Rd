
\name{HpyloriData}
\alias{HpyloriData}

\title{
    H. pylori data
}
\description{ 
    SNP positions and allele frequencies for the reference data
}
\details{
    A \code{\link[GenomicRanges]{GPos}} class object containing
    the vertex and reference population 
    allele frequencies for the
    set of 143705 SNPs used in the analysis for H. pylori.
    The SNPs were created using 
    26695 (NCBI GenBank Accession NC_000915.1) as the 
    reference genome.
    The set of SNPs was selected using a MAF threshold of 0.01.
    The total sample size was from a set of 1011 H. pylori strains.
}
\value{ An object of class \code{\link[GenomicRanges]{GPos}}.}
\examples{
    # Load data and view the first few rows
    data(HpyloriData, package="GrafGen")
    HpyloriData
}
\keyword{data}
