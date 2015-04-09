\name{vegawes}
\alias{vegawes}

\title{
Variational Segmentation on WES Data
}
\description{
This function computes the variational segmentation on WES data based on the Mumford Shah variational model.
}
\usage{
vegawes(CNVdata, chromosomes, beta = 0.5, alpha = 0.001, min_region_size = 2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{CNVdata}{This argument is matrix containing the data that have to be analyzed. This matrix must have 4 columns:\cr 
  	- the first row indicates the chromosome;\cr
		- the second row indicates the start bp of the probe;\cr
		- the third row indicates the end bp of the probe;\cr
		- the fourth row reports the measured Log R Ratio;\cr

}

  \item{chromosomes}{This is a vector containing the chromosome sthat have to be analyzed. By using c(1:22, "X", "Y") the whole genome will be segmented.
  }
  \item{alpha}{This argument is the weight given the distance between exomes factor in the segmentation algorithm.}
  \item{beta}{This argument is used for the stop condition definition.}
  \item{min_region_size}{This argument specifies the minimum size for the segmented regions.}
}


\value{
  \item{segmentation}{This is a data frame containing the computed segmentation.}
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Samreen Anjum, Michele Ceccarelli, Sandro Morganella
}


\examples{
\dontrun{
vegawes(CNVData, chromosomes=c(1:22), alpha = 0.001, beta = 0.7)
}
}

