
\name{GrafGen-package}
\alias{GrafGen}
\alias{GrafGen-package}

\docType{package}
\title{
    Classify H. pylori genomes
}
\description{
    To classify H. pylori genomes according to genetic distance 
    from nine reference populations.
}
\details{
    This package was modified from the GrafPop software\cr
    \bold{(https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/Software.cgi)}\cr
    to be 
    applied on H. pylori genomes. The three vertex populations are
    "Africa", "Europe" and "Asia".
    The nine reference populations are "hpgpAfrica", "hpgpAfrica-distant",
    "hpgpAfroamerica", "hpgpEuroamerica", "hpgpMediterranea", "hpgpEurope",
    "hpgpEurasia", "hpgpAsia", and "hpgpAklavik86-like". 
    The training data is based on The Helicobacter pylori Genome 
    Project (HpGP), see \cr
    \bold{https://www.ncbi.nlm.nih.gov/bioproject/?term=HpGP} or \cr
    \bold{https://zenodo.org/records/10048320}.

    To use this package, the user must have a file of genotypes for 
    H. pylori strains. The genotype file can be a binary PLINK file
    in SNP-major format, or a VCF file of genotypes. If a PLINK file,
    then the corresponding \code{bim} and \code{fam} files must also
    be present. If a VCF file, then the format should be genotypes:\cr
    \code{##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">}.

    Ideally, the
    genotype file will contain all the SNPs with positions given
    in the \code{\link{HpyloriData}} data frame, where the positions
    are based on the reference genome 
    26695 (NCBI GenBank Accession NC_000915.1). However, the
    software has been shown to work well with only using a much smaller
    fraction of SNPs in \code{\link{HpyloriData}}. The main function
    in this package is \code{\link{grafGen}}.
}
\references{
    Jin Y, Schaffer AA, Feolo M, Holmes JB and Kattman BL (2019). 
    GRAF-pop: A Fast Distance-based Method to Infer Subject Ancestry 
    from Multiple Genotype Datasets without Principal Components Analysis.
    G3: Genes | Genomes | Genetics. DOI: 10.1534/g3.118.200925. 

    Thorell K, Munoz-Ramirez ZY, Wang D, Sandoval-Motta S, 
    Boscolo Agostini R, Ghirotto S, Torres RC, HpGP Research Network, 
    Falush D, Camargo MC and Rabkin CS (2023). 
    New insights into Helicobacter pylori population structure from 
    analysis of a worldwide collection of complete genomes: the 
    H. pylori genome project.  
    Nature Communications. DOI: 10.1038/s41467-023-43562-y.
}
\author{William Wheeler, Difei Wang, Isaac Zhao, Yumi Jin, Charles Rabkin}
\keyword{package}
