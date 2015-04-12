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

\references{ 
Anjum, S., Ceccarelli, M. (2015). VEGAWES: Variational segmentation on Whole Exome Sequencing, Submitted. \cr
Morganella, S., Cerulo, L., Viglietto, G., & Ceccarelli, M. (2010). VEGA: Variational segmentation for copy number detection. Bioinformatics, 26(24), 3020-3027.
}
\author{Samreen Anjum, Michele Ceccarelli\cr
Maintainer: Samreen Anjum <samreen.anjum@gmail.com>}

\examples{
\dontrun{
# Start VEGAWES segmentation on chromosomes 8 of the data with alpha = 0.001, beta = 0.7 and save the results in chr8.txt
runVEGAWES(inputfile, chr.list=c(8), alpha = 0.001, beta=0.7);
}
}
