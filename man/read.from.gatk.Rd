\name{read.from.gatk}
\alias{read.from.gatk}

\title{
Read GATK files
}
\description{
This function reads in the GATK file and creates an object containing the required information.
}
\usage{
read.from.gatk(gatk.file)
}

\arguments{
  \item{gatk.file}{The output file from the GATK computation
}
}

\value{
\item{probe}{The probe}
\item{chr}{The chromosome number}
\item{probe_start}{Start coordinate of the probe in the genome }
\item{probe_end}{End coordinate of the probe in the genome }
\item{targeted.base}{The length of the probe }
\item{coverage}{The total read coverage of the probe}
\item{average.coverage}{The average read coverage of the probe}
          
}

\author{
Samreen Anjum, Michele Ceccarelli
}


\examples{
\dontrun{
normal_RCData = read.from.gatk(normal.gatk.file)
}
}
