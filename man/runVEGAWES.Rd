\name{runVEGAWES}
\alias{runVEGAWES}
\title{Function to compute variational copy number segmentation on WES data}

\description{
 This is the main function that computes segmentation on the WES data and is an extended version of 'vega'. It reads in an inputfile containing the filepaths to WES BAM files and calls GATK to compute average read counts. It then normalizes the data and removes GC Content Bias. Finally, it computes segmentation on the normalized data and saves it as a tab-delimited file. In order to use this function, users have to modify the parameters.txt file to include the paths to their files and make the necessary path changes to other required data. 
  

}
\usage{
runVEGAWES(inputfile, chr.list = c(1:22), alpha = 0.001, beta = 0.7)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{inputfile}{This argument is path to a file that has 8 lines and consists the following information:\cr
    - Path to the working directory\cr
    - Reference genome file path\cr
    - Path to the GATK jar\cr
    - List of exons to run GATK (chr, probe_start, probe_end)\cr
    - Name of the sample - to create an output folder to save results\cr
    - Normal input BAM file path\cr
    - Tumor input BAM file path\cr
    - Path to GC Content Data folder\cr

}

  \item{chr.list}{This is a vector containing the chromosomes that have to be analyzed. Setting chr.list= c(1:22) the whole genome will be segmented.
  }
  \item{alpha}{This is the weight parameter used for the distance factor in the segmentation. Settling alpha = 0 will compute the original vega segmentation }
  \item{beta}{This argument is used for the stop condition definition.}
}



\author{
Samreen Anjum, Michele Ceccarelli
}

\examples{
\dontrun{
# Start VEGAWES segmentation on all chromosomes and save it in chr1txt, ..., chr22.txt
runVEGAWES(inputfile, chr.list=c(1:22));

# Start VEGAWES segmentation on chromosomes 8 of the data with alpha = 0.001, beta = 0.7 and save the results in chr8.txt
runVEGAWES(inputfile, chr.list=c(8), alpha = 0.001, beta=0.7);
}
}
