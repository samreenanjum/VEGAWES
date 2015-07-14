vegawes <- function(CNVdata, chromosomes,  beta=0.5, alpha = 0.001, min_region_size=2){
  # Compute the segmentation for the data contained in 'CNVdata'
  # 'chromosomes' indicates the chromosomes that have to be analyzed
  # At the end of the computation the segmentation is saved in the tab delimited file 'out_file_name'
  
  
  # Create the structure of the output segmentation file
  segmentation <- matrix(0,0,6);
  colnames(segmentation) <- c("Chromosome", "bp Start", "bp End", "Num of Markers", "Mean", "Label");
  
  for( i in 1:length(chromosomes) ){
    
    curr_chr <- chromosomes[i];	
    chr_ids <- which(CNVdata[,1]==curr_chr);
    data_table <- CNVdata[chr_ids,];
    
    
    message("Processing Chromosome ", curr_chr);
    
    # Set the parameters 
    n <- nrow(data_table);
    al <- c(alpha);
    be <- c(beta);
    m <- c(min_region_size);
    start <- (0*c(1:(n+1)));
    end <- (0*c(1:(n+1)));
    size <- (0*c(1:(n+1)));
    mean <- (0*c(1:(n+1)));
    label <- (0*c(1:(n+1)));
    n_reg <- 0;
    std <- sd(data_table[,4]);
    del <- (0*c(1:n))-1;
    
    seg <- .C("vegawes", data=as.double(t(data_table[,4])),
              markers_start = as.integer(t(data_table[,2])),
              markers_end = as.integer(t(data_table[,3])),
              positions = as.double(((as.integer(t(data_table[,2])))+(as.integer(t(data_table[,3]))))/2),
              start=as.integer(start), 
              end=as.integer(end), 
              size=as.integer(size),
              mean=as.double(mean), 
              label=as.integer(label), 
              n = as.integer(n),
              al= as.double(alpha),
              be= as.double(beta),
              m = as.integer(min_region_size),
              std = as.double(std),
              n_reg = as.integer(n_reg)
    );
    # Save the segmentation for the current chromosome
    
    curr_seg = matrix(0, seg$n_reg, 6);
    curr_seg[,1] <- as.character(curr_chr);
    curr_seg[,2] <- seg$start[1:seg$n_reg];
    curr_seg[,3] <- seg$end[1:seg$n_reg];
    curr_seg[,4] <- seg$size[1:seg$n_reg];
    curr_seg[,5] <- seg$mean[1:seg$n_reg];
    curr_seg[,6] <- seg$label[1:seg$n_reg];
    segmentation <- rbind(segmentation, curr_seg);
    message("Done\n");
  }
  
  return(segmentation);	
  
}

