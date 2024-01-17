
\name{grafGen}
\alias{grafGen}

\title{Ancestry for  H. pylori strains}
\usage{
    grafGen(genoFile, print=1)
}
\arguments{
\item{genoFile}{The complete path to the input genotype file.
    This file can only be a PLINK binary file (.bed)
    or a VCF file (.vcf, .vcf.gz). If it is a .bed file,
    then the corresponding .bim and .fam files
    must also exist.
    If a VCF file, then the format should be genotypes:\cr
    \code{##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">}.
    }
\item{print}{0 or 1 to print information as the program runs.}
}
\description{
    To determine the ancestry of H. pylori strains.
}
\details{
    See the references for complete details of the algorithm.

    This function is more efficient if the input genotype file only 
    contains the set (or subset) of SNPs defined in \code{\link{HpyloriData}}.
    The SNPs can be extracted by calling \code{\link{extractPositions}} if
    the genotype file is a VCF file. For a binary PLINK file, the PLINK
    software must be used to extract the SNPs. 
}

\value{A list of class "grafpop" containing a data frame (\code{table})
    that includes the ancestry percents 
    (\code{A_percent, F_percent, E_percent}) for
    Asian, African, and European respectively,  
    normalized genetic distance scores (\code{GD1_x, GD2_y, GD3_z}),
    the predicted reference population (\code{Refpop}), and
    the genetic distances to each reference population 
    (names of the reference populations) as defined by equation 2 in
    Jin (2019).  
    The list also includes the list \code{vertex} which
    gives the x-y coordinates of the vertex populations.
}

\examples{
    dir  <- system.file("extdata", package="GrafGen", mustWork=TRUE)
    file <- file.path(dir, .Platform$file.sep, "data.vcf.gz")
    grafGen(file)
}

\seealso{ \code{\link{extractPositions}} }
