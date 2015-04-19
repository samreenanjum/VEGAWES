# VEGAWES package

This is an R package that implements VEGAWES, a segmentation algorithm based on Mumford and Shah variational model, for copy number segmentation on Whole Exome Sequencing data. 
The package and the source code can be downloaded [here](https://github.com/samreenanjum/VEGAWES). 

##### To install the package on R:

    library(devtools)
    install_github("samreenanjum/VEGAWES")

##### Load the library:

    library(VEGAWES)

##### Dependencies:
* Hg19 Reference Genome: (ftp://ftp.ncbi.nlm.nih.gov/sra/reports/Assembly/GRCh37-HG19_Broad_variant/Homo_sapiens_assembly19.fasta)
* Genome Analysis Tool Kit (GATK) jar


##### Input files required: 

* Paired Tumor-Normal BAM files
* Exome Interval List (Included in the package as "targets.interval\_list" under /inst/extdata)
* GC Content Data (Included in the package under /data)

In addition, "parameters.txt" must be modified to include all the required paths. See Examples below and VEGAWES.pdf under /inst/doc for further details.

The file "parameters.txt" contains the following information:

    * Path to the working directory
    * Reference genome file path
    * Path to the GATK jar
    * List of exons filepath to run GATK
    * Name of the sample (to create an output folder under the working directory to save GATK results and the segmentation results)
    * Normal input BAM file path
    * Tumor input BAM file path
    * Path to the GC Content Data folder

Example of "parameters.txt":

    /home/samreen/VEGAWES
    /home/samreen/VEGAWES/Homo\_sapiens\_assembly19.fasta
    /home/samreen/tools/GATK/GenomeAnalysisTK.jar
    /home/samreen/VEGAWES/inst/extdata/targets.interval_list
    06-0125
    /home/samreen/data/TCGA/GBM/TCGA-06-0125/6f138046-a805-489a-a335-6686173d6505/C484.TCGA-06-0125-10A-01D-1490-08.6.bam
    /home/samreen/data/TCGA/GBM/TCGA-06-0125/ced14238-6593-4648-b3bb-ca7711f71346/C484.TCGA-06-0125-01A-01D-1490-08.6.bam
    /home/samreen/VEGAWES/data/GCC



##### Run the pipeline

To run the pipeline, the user needs to run the `runVEGAWES` function giving atleast two parameters - "parameters.txt" and the list of chromosomes to be analyzed. In addition, Java must be in the current path to run GATK.

This function first runs GATK to compute average read counts from the BAM files, normalizes the counts based on GC Content, computes log-ratio values for each exome, and finally runs VEGAWES segmentation on the required chromosomes. The average read counts and other output files are stored in a new output directory named after the sample-name. 

The function can be called in the following way:

    runVEGAWES( "parameters.txt", chr.list = c(1:22) )

The segmentation results for each chromosome are stored in the output directory as "Segmentation.\<chromosome\>.txt". 





Please refer to the VEGAWES.pdf under /inst/doc for further details and descriptions.
