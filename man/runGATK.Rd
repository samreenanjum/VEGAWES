\name{runGATK}
\alias{runGATK}

\title{
Compute Average Read Count using GATK on the BAM files
}
\description{
This function runs the GATK command with the DepthOfCoverage utility to compute the average read coverage of the normal and tumor BAM files and saves it in the output folder named after the sample name. Assumption: Java exists in the current path
}
\usage{
runGATK(params)
}

\arguments{
  \item{params}{
This argument is the output object from the function 'initializeParams' containing the path to the GATK jar, the exon list, the reference genome, the sample name as well the input BAM files. 
}
}


\author{
Samreen Anjum, Michele Ceccarelli
}

\examples{
\dontrun{
#params is the output of the function 'initializeParams'
runGATK(params)
}
}

