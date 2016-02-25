## Helpful functions used in VEGAWES pipeline


######################### Initialize all the parameters and file paths ##########################################
initializeParams <- function(inputfile){
  
  params.data = read.table(inputfile) 
  params.data = as.character(params.data$V1)
  ## working directory - Must have targets.interval_list, fasta file
  working.folder = params.data[1]
  reference.file = params.data[2]
  
  ## Parameters for GATK
  path.to.GATK = params.data[3]
  interval.list.file = params.data[4]
  
  sample.name = params.data[5]
  
  ## filenames of tumor and normal BAM files
  normal.inputfile = params.data[6]
  tumor.inputfile = params.data[7]
  
  ## Create directory for output
  dir.create(paste(working.folder, "/",sample.name, sep=""), showWarnings=FALSE)
  
  #filenames of GATK outputfiles
  normal.outfile = paste(working.folder,"/",sample.name,"/normal-", sample.name, sep="")
  tumor.outfile = paste(working.folder, "/",sample.name,"/tumor-", sample.name, sep="")
  
  GCContent.folder = params.data[8]
  return ("params" = data.frame(working.folder, reference.file, path.to.GATK, interval.list.file, sample.name, normal.inputfile, tumor.inputfile, normal.outfile, tumor.outfile, GCContent.folder))
}
########################################### End of Initialization ##########################################


################# Invoke the GATK command ################ 
runGATK <- function(params){
  
  # Run the GATK command
  commandString = paste("java -jar", params$path.to.GATK ,"-T DepthOfCoverage -omitBaseOutput -omitLocusTable -R", params$reference.file , "-L", params$interval.list.file, sep=" ")  
  
  system(paste(commandString, "-I", params$normal.inputfile, "-o", params$normal.outfile))  #Normal File
  system(paste(commandString, "-I", params$tumor.inputfile, "-o",params$tumor.outfile))     #Tumor File
  
}
################## End of GATK  ############################ 



########################### Read Data from the GATK format ################################
read.from.gatk <- 
  function(gatk.file) {
    gatk.data = read.table(gatk.file, header=TRUE)
    chrpos = matrix(unlist(strsplit(as.character(gatk.data$Target),":")), ncol=2, byrow=TRUE)
    ind = which(!(grepl("-",gatk.data$Target)))  ## get the targets with inconsistent format (due to size = 1)
    chrpos[ind,2]  = paste0(chrpos[ind,2], "-", chrpos[ind,2])
    chr = factor(chrpos[,1])
    pos = matrix(as.integer(unlist(strsplit(chrpos[,2],"-"))), ncol=2, byrow=TRUE)
    start = pos[,1]
    end = pos[,2]
    return(data.frame( probe=gatk.data$Target, 
                       chr=chr, 
                       probe_start=start, 
                       probe_end=end, 
                       targeted.base=end-start+1, 
                       coverage=as.numeric(gatk.data$total_coverage), 
                       average.coverage=as.numeric(gatk.data$average_coverage)
    ))
  }
#############################################################################


######################## Compute copy numbers for the given logR values ##############
calculateCN <- function(logR){
  CN = round(logR)
  CN[which(logR< -0.36)] = 1
  CN[which(logR<0.24&logR>=-0.36)] = 2
  CN[which(logR>=0.24)] = 3
  
  return(CN)
}


#####################################################################
## Correct the average read coverage values using the GC Content ###
## Adapted from EXCAVATOR tool 
####################################################################
correct.GCContent<-function(average.coverage,GCContent,step)
{
 
  stepseq<-seq(0,100,by=step)
  
  master.median<-median(average.coverage,na.rm=T)
  median.GC<-rep(0,length(stepseq)-1)
  average.coverage.normalized<-average.coverage
  for (i in 1:(length(stepseq)-1))
  {
    if (i==1) { ind<-which(GCContent>=stepseq[i] & GCContent<=stepseq[i+1]) }
    else { ind<-which(GCContent>stepseq[i] & GCContent<=stepseq[i+1])  }
   
    if (length(ind)>0)
    {
      m<-median(average.coverage[ind],na.rm=T)
      if (m>0)
      {
        median.GC[i]<-m
        average.coverage.normalized[ind]<-average.coverage[ind]*master.median/m
      }
    }
  } 
  return(average.coverage.normalized)
}



