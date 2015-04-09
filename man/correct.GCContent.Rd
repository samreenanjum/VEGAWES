\name{correct.GCContent}
\alias{correct.GCContent}

\title{
GC Bias correction of exon average read coverage
}
\description{
This function adjusts the average read coverage values based on the GC Content. The approach is based on the median normalization approach described in Yoon et al.
}
\usage{
correct.GCContent(average.coverage, GCContent, step)
}

\arguments{
  \item{average.coverage}{
This argument is the original average read coverage values
}
  \item{GCContent}{
This argument is the percentage of GCContent in an exon
}
  \item{step}{
The argument is the interval of GC percentage values chosen to compute the median to adjust each exon
}
}

\value{
 \item{average.coverage.normalized}{The adjusted values of average read coverage based on the GC normalization}
}
\references{
Magi, Alberto, et al. "EXCAVATOR: detecting copy number variants from whole-exome sequencing data." Genome biology 14.10 (2013): R120.\cr
Yoon, Seungtai, et al. "Sensitive and accurate detection of copy number variants using read depth of coverage." Genome research 19.9 (2009): 1586-1592.

}
\author{
Samreen Anjum, Michele Ceccarelli
}

\examples{
\dontrun{
average.coverage.normalized = correct.GCContent(average.coverage, GCContent, step=10)
}
}
