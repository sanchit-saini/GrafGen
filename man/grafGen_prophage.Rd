
\name{grafGen_prophage}
\alias{grafGen_prophage}

\title{GrafGen for Prophages}
\usage{
    grafGen_prophage(genoFile, referenceData=NULL, print=1)
}
\arguments{
\item{genoFile}{The complete path to the input genotype file.
    A VCF file (.vcf, .vcf.gz) with genotypes:\cr
    \code{##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">}.
    }
\item{referenceData}{(Only to use a different set of training loci)
    Data frame with at least 9 columns where for
    a single chromosome haploid genome,\cr
    Column 1 - SNP position \cr
    Column 2 - SNP non-variant allele\cr
    Column 3 - SNP variant alleles\cr
    Column 4 - Vertex population 1 allele frequency\cr
    Column 5 - Vertex population 2 allele frequency\cr
    Column 6 - Vertex population 3 allele frequency\cr
    Column 7 - Reference population 1 allele frequency\cr
    Column 8 - Reference population 2 allele frequency\cr
    Column 9 - Reference population 3 allele frequency\cr
    Columns 10 through the last column are for other reference
    population allele frequencies. The order of the first three
    reference populations must match the order of the three
    vertex populations.\cr
    For multi-chromosome haploid genomes, the first column must 
    be a \bold{character} column for the chomosome id, and the other
    columns must be as above starting on column 2.
    }
\item{print}{0 or 1 to print information as the program runs.}
}
\description{
    To determine the ancestry of prophage genomes.
}
\details{
    See the references for complete details of the algorithm.

    For non-H. pylori prophage data, the data frame \code{referenceData} must be passed in.
    All the allele frequencies in this file must be greater than zero and
    less than one.
    The column names defining the vertex and reference populations 
    in \code{referenceData} will appear in the output and plots.
}

\value{A list of class "grafpop" containing a data frame (\code{table})
    that includes the ancestry percents,  
    normalized genetic distance scores (\code{GD1_x, GD2_y, GD3_z}),
    the predicted reference population (\code{Refpop}), 
    next nearest reference population (\code{Nearest_neighbor}),
    separation to the next nearest reference population 
    (\code{Separation_percent}) defined as \code{100*abs(d1 - d2)/d1},
    where \code{d1} and \code{d2} are the genetic distances to the 
    sample's assigned reference population and next nearest reference
    population respectively, and
    the genetic distances to each reference population 
    as defined by equation 2 in Jin (2019).
    For H. pylori prophages, the ancestry percents are the columns 
    \code{F_percent, A_percent, E_percent} for percent
    African, East Asian, and southwest European respectively, and the 
    genetic distances to each reference population are 
    \code{Africa1}, \code{EastAsia},
    \code{SWEurope}, and \code{NEurope}.

    The returned object also includes the list \code{vertex} which
    gives the x-y coordinates of the vertex populations.
}

\examples{
    dir  <- system.file("extdata", package="GrafGen", mustWork=TRUE)
    file <- file.path(dir, "data_prophage.vcf.gz")
    grafGen_prophage(file)
}

\seealso{ \code{\link{grafGen}} }
