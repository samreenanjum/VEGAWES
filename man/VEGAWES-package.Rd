\name{VEGAWES-package}
\alias{VEGAWES-package}
\alias{VEGAWES}
\docType{package}
\title{An R package for copy number data segmentation on WES data}
\description{VEGAWES (Variational Estimator for Genomic Aberrations on Whole Exome Sequencing data) is an algorithm that adapts a very popular variational model (Mumford and Shah) and used in image segmentation so that chromosomal aberrant regions can be efficiently detected.
}
\details{
\tabular{ll}{
Package: \tab VEGAWES\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
License: \tab GPL-2\cr
}
}

\references{ Morganella S. et al. (2010). VEGA: Variational segmentation for copy number detection, Bioinformatics.}
\author{Samreen Anjum, Michele Ceccarelli\cr
Maintainer: Samreen Anjum <samreen.anjum@gmail.com>}

\examples{
\dontrun{
# Run VEGAWES segmentation on all chromosomes of the WES data and save the results in a .txt tab delimited file
vegawes(CNVData, chromosomes=c(1:22), alpha = 0.001, beta=0.5);
}
}
