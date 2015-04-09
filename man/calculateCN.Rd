\name{calculateCN}
\alias{calculateCN}

\title{
Copy Number Classication
}
\description{
This function computes the copy number value for the given LRR values based on the following classification: <0.36 for amplication, >0.24 for deletion, normal otherwise.
}
\usage{
calculateCN(logR)
}

\arguments{
  \item{logR}{
This argument is vector containing all the logR values computed by the segmentation algorithm
}
}

\value{
Returns 2 for normal, 3 for amplification, and 1 for deletion 
}

\author{
Samreen Anjum, Michele Cecceralli
}


\examples{
\dontrun{
#calculate copy number classification based on the logR values
calculate(logR)

}
}