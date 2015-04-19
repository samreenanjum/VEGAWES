# VEGAWES package

This is an R package that implements VEGAWES, a variational model based segmentation algorithm, for copy number segmentation on Whole Exome Sequencing data. 
The package and the source code can be downloaded from https://github.com/samreenanjum/VEGAWES. 

##### To install the package on R:

install.packages("https://github.com/samreenanjum/VEGAWES/VEGAWES_1.0.tar.gz")

##### Load the library:

library(VEGAWES)

##### Input files required for running this pipeline: 

* Paired Tumor-Normal BAM files

* Reference Genome: ftp://ftp.ncbi.nlm.nih.gov/sra/reports/Assembly/GRCh37-HG19_Broad_variant/Homo_sapiens_assembly19.fasta

* GATK jar

In addition, the file included in the package - "parameters.txt" must be modified to include all the required paths. 

##### Run the pipeline

To run the pipeline, the user needs to run the runVEGAWES function giving atleast two parameters - "parameters.txt" and the list of chromosomes to be analyzed. This function first runs GATK to compute average read counts from the BAM files, normalizes the counts based on GC Content, computes log-ratio values for each exome, and finally runs VEGGAWES segmentation on the required chromosomes. The average read counts and other output files are stored in a new output directory named after the sample-name. The function can be called in the following way:

> runVEGAWES( "parameters.txt", chr.list = c(1:22) )

The segmentation results for each chromosome are stored in the output directory as "Segmentation.\<chromosome\>.txt". 


Please refer to the VEGAWES.pdf under /inst/doc for further details and descriptions.
