
### libraries
library(ggplot2)
library(data.table)
library(R.utils)
library(patchwork)
library(foreach)

### read the SYNC file
sync <- fread("/scratch/ygn8wq/ExpEvo_PRJNA304655_AK2_1_2009-03.tab")

setnames(sync, names(sync), c("chr", "pos", "ref", "A", "T", "C", "G", "N","del"))
sync <- sync[,-"del",]
sync[,id:=1:dim(sync)[1]]

### clean up a bit and only use the SNPs on the major chromosomes
setkey(sync, chr)
sync <- sync[J(c("2L", "2R", "3L", "3R", "X"))]

### Our first challenge is filtering down the SYNC file to only retain the polymorphic sites.

### One strategy is to use our aggregate skills to identify sites where there are less than three zeros
### in the A, C, T, G columns. This solution takes forever. How patient are you?
sync.ag <- sync[,list(nZero=sum(A==0 & T==0 & C==0 & G==0)), list(chr, pos)]

### another strategy is to ask whether the depth of the reference allele is equal to the total depth
sync[,depth:=A+T+C+G+N]

setkey(sync, ref)
sync.poly <- foreach(ref.allele=c("A", "C", "T", "G"))%do%{
  # ref.allele <- "A"
  message(ref.allele)
  tmp <- sync[J(ref.allele)] ### this line subsets the data to only sites where the ref is "ref.allele"
  setnames(tmp, ref.allele, "focal") ##
  polymorphic <- tmp[focal!=depth]
  polymorphic[,ref_count:=focal]
  setnames(polymorphic, "focal", ref.allele)
  polymorphic
}
sync.poly <- rbindlist(sync.poly)
sync.poly[,alt_freq:=1 - ref_count/depth]

### plot a histogram
ggplot2(data=sync.poly, aes(alt_freq, group=chr, color=chr)) + geom_density() + facet_grid(~chr)

### load in Gene definition
genes <- fread("/standard/BerglandTeach/data/genes.bed")
setnames(genes, names(genes), c("chr", "start", "stop", "gene"))
genes[,stop:=as.numeric(stop)]
genes <- na.omit(genes)
genes[,gene_region:=T]

### YOUR TURN:

sync.poly[,start:=pos]
sync.poly[,stop:=pos]

setkey(genes, chr, start, stop)
setkey(sync.poly, chr, start, stop)

sync.poly <- foverlaps(sync.poly, genes, nomatch=NA)
sync.poly[is.na(gene_region), gene_region:=F]

sync.poly[,list(meanRD=mean(depth), maxRD=max(depth)), list(chr, gene_region)]

coverage_density_poly_plot <- ggplot() + 
  geom_density(data=sync.poly, aes(alt_freq, group=gene_region, color=gene_region), alpha=0.5) + geom_density() + facet_grid(~chr)
  
coverage_density_poly_plot
### using the commands from our Coverage assignment, merge the genes file with the polymorphic SYNC file and make the output figure from the instructions