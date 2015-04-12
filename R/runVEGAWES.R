source("vegawes.functions.R")

#################################################################################
##                                                                             ##
## Runs the VEGAWES segmentation on WES data                                    ##
## Assumptions: Java exists in the current path                                ##
##                                                                             ##
## Input:  inputfile - File containing all required parameters                 ##
##         chr.list  - Chromosomes on which segmentation is performed          ##
##         alpha     - Weight given to the distance factor in the segmentation ##
##         beta      - Argument used for the stop condition definition         ##
## Output: files containing segmentation information for each chromosome       ##
##                                                                             ##
#################################################################################

runVEGAWES <- function(inputfile, chr.list = c(1:22), alpha = 0.001, beta = 0.7){
  
  params = initializeParams(inputfile)  
  normal.gatk.file = paste(params$normal.outfile, ".sample_interval_summary", sep="")
  tumor.gatk.file = paste(params$tumor.outfile, ".sample_interval_summary", sep="")
  
  
  ########## Run GATK to compute ARC values if they do not exist ###################
  if(!(file.exists(normal.gatk.file) && file.exists(tumor.gatk.file))){
    print("Running GATK...")
    runGATK(params)
  }
  
  ############################## Read in the normal-tumor RC files and normalize the data ###############
  print("Reading RC files from GATK and normalizing...")
  normal_RC = read.from.gatk(normal.gatk.file)
  tumor_RC = read.from.gatk(tumor.gatk.file)
  
  total.cov.normal = sum(as.numeric(normal_RC$coverage), na.rm=TRUE)
  total.cov.tumor = sum(as.numeric(tumor_RC$coverage), na.rm=TRUE)
  ratio = total.cov.tumor/total.cov.normal
  normal_RC$average.coverage = normal_RC$average.coverage*ratio
  ####################################### End of Normalization ##########################################
  
  
  ######################################## Processing each chromosome ###########################
  for (chr in chr.list){
    print(paste("Processing", chr ))
    
    ## Get the coverage values for this chromosome
    normal = normal_RC[normal_RC$chr==chr,]
    tumor = tumor_RC[tumor_RC$chr==chr,] 
        
    #####  GC Bias Removal #####
    load(paste(params$GCContent.folder,"//GCContent.",(chr),".RData", sep=''))
    
    ## Get the required exons
    GCData = GCData[GCData$probe_start %in% normal$probe_start,]
    
    # Compute the GC Content and update the Average Read Coverage (ARC) values
    GCContent.original <- GCData$GCContent       # Get the GCContent
    exon.lengths = normal$targeted.base             # Get the length of the exons
    GCContent.percent <- (GCContent.original/exon.lengths)*100
    step = 10    
    
    normal$average.coverage <- correct.GCContent(normal$average.coverage,GCContent.percent,step)
    tumor$average.coverage <- correct.GCContent(tumor$average.coverage,GCContent.percent,step)
    
    
    ##### Compute the LRR values ######
    logR =  log2(tumor$average.coverage/normal$average.coverage) 
    
    
    ########  Prepare the data for VEGAWES and run the segmentation ##########
    coverage.cutoff = 5
    covered.exons = (normal$average.coverage > coverage.cutoff) & (tumor$average.coverage > coverage.cutoff)
    vegaInput = cbind(as.numeric(as.character((normal$chr[covered.exons]))),normal$probe_start[covered.exons], normal$probe_end[covered.exons],logR[covered.exons])  
    segs = vegawes(vegaInput, chromosomes=unique((normal$chr[covered.exons])), alpha = as.double(alpha), beta = as.double(beta))
    
    cnv = data.frame(chr=segs[,1], probe_start=as.numeric(segs[,2]), probe_end=as.numeric(segs[,3]), num.mark=as.numeric(segs[,4]), logR=as.numeric(as.character(segs[,5])))
    cnv$copy.number = calculateCN(as.numeric(as.character(cnv$logR)))
    
    ######  Write it to file  ########
    write.table(cnv, file=paste(params$working.folder, "//", params$sample.name, "//", chr,".txt",sep=""), quote=FALSE, sep="\t", row.names=FALSE, )
  }  
}





