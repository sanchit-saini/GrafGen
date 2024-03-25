
\name{grafGen}
\alias{grafGen}

\title{Reference population for  H. pylori strains}
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
    The SNPs can be extracted by utilizing the VCFtools software if
    the genotype file is a VCF file. For a binary PLINK file, the PLINK
    software can be used to extract the SNPs. 
}

\value{A list of class "grafpop" containing a data frame (\code{table})
    that includes the ancestry percents 
    (\code{F_percent, E_percent, A_percent}) for
    African, European and Asian respectively,  
    normalized genetic distance scores (\code{GD1_x, GD2_y, GD3_z}),
    the predicted reference population (\code{Refpop}), 
    next nearest reference population (\code{Nearest_neighbor}),
    separation to the next nearest reference population 
    (\code{Separation_percent}) defined as \code{100*abs(d1 - d2)/d1},
    where \code{d1} and \code{d2} are the genetic distances to the 
    sample's assigned reference population and next nearest reference
    population respectively, and
    the genetic distances to each reference population 
    (\code{hpgpAfrica}, \code{hpgpAfrica-distant},
    \code{hpgpAfroamerica}, \code{hpgpEuroamerica}, \code{hpgpMediterranea}, 
    \code{hpgpEurope}, \code{hpgpEurasia}, \code{hpgpAsia}, and 
    \code{hpgpAklavik86-like}) as defined by equation 2 in
    Jin (2019).
    The returned object also includes the list \code{vertex} which
    gives the x-y coordinates of the vertex populations.
}

\examples{
    dir  <- system.file("extdata", package="GrafGen", mustWork=TRUE)
    file <- file.path(dir, "data.vcf.gz")
    grafGen(file)
}

\seealso{ \code{\link{HpyloriData}} }
