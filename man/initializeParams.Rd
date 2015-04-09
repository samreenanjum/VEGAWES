\name{initializeParams}
\alias{initializeParams}

\title{
Initialization of filenames and parameters
}
\description{
This function reads in the "parameters.txt" file, creates an output folder named after the sample name and returns an object that contains all the required paths and parameters.

}
\usage{
initializeParams(inputfile)
}

\arguments{
  \item{inputfile}{
This argument is the "parameters.txt" file the contains the required filepaths and samplename.
}
}

\value{
params The object containing the path information
 \item{working.folder}{ The current working directory containing the reference genome and where the output will be saved}
 \item{reference.file}{Path to the reference genome - Fasta file}
 \item{path.to.GATK}{Path to the GATK jar}
 \item{interval.list.file}{Path to the file containing the exon interval list}
 \item{sample.name}{The name of the sample to create the output folder directory and save the GATK and segmentation results, eg: '06-0125'}
 \item{normal.inputfile}{Path to the input normal BAM file}
 \item{tumor.inputfile}{Path to the input tumor BAM file}
 \item{GCContent.folder}{Path to the folder containing the GCContent .Rdata}

}

\author{
Samreen Anjum, Michele Ceccarelli
}


\examples{
\dontrun{
params = initializeParams("parameters.txt")
}
}
