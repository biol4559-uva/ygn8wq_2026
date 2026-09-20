### install a new package; you only need to do this once.
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("Rsamtools")

library(ggplot2)
library(data.table)
library(foreach)
library(doMC)
registerDoMC(2)
library(Rsamtools)

### specify the bam file
fl2 <- system("ls -d /standard/BerglandTeach/mapping_output/*PRJNA304655*/*.original.bam", intern=T)
rd2 <- foreach(bamFile=fl2, .combine="rbind")%dopar%{
  ## bamFile=fl2[1]
  
  ### tell me what file we are working on
  message(bamFile)
  
  ### get the information about number of reads for each chromosome
  if(!file.exists(paste(bamFile, ".bai", sep=""))) indexBam(bamFile)
  stats <- as.data.table(idxstatsBam(bamFile))
  
  ### to subset to the main autosomal arms and X chromosome
  setkey(stats, seqnames)
  stats_small <- stats[J(c("2L", "2R", "3L", "3R", "X"))]
  
  ### separate autosomes and X chromosome
  autosomes <- stats_small[seqnames %in% c("2L", "2R", "3L", "3R")]
  xchrom <- stats_small[seqnames=="X"]
  
  
  ### coverage of the autosomes
  auto_cov <- sum(autosomes$mapped) / sum(autosomes$seqlength) #variables from stats data table
  ### coverage of the X chromosome
  x_cov <- sum(xchrom$mapped) / sum(xchrom$seqlength)
  
  
  ### X chromosome coverage to autosomal coverage ratio
  x_auto_ratio <- x_cov / auto_cov
  
  ### proportion female
  propFemale <- 2*x_auto_ratio - 1
  
  ### proportion male
  propMale <- 1 - propFemale
  
  ### format output
  out <- data.table(x_auto_ratio=x_auto_ratio, propFemale=propFemale, propMale=propMale)
  
  out[,samp:=last(tstrsplit(bamFile, "/"))]
  
  
  ### return output
  return(out)
}


### results
rd2


a