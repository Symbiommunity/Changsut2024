##code for differential gene expression analyses associated with Changsut et al. 2023##
##last edited by LEF 4/20/2024##

library(DESeq2)
library(ggplot2)
library(dplyr)
library(stringr)

##pre-processing
  ##read in and format our allcounts.txt file
  reads=read.delim(file="allcounts.txt", header=TRUE, sep = "\t")
  colnames(reads)=gsub(".trim.sam.counts", '', colnames(reads), fixed=TRUE)
  colnames(reads) = gsub(pattern = "\\.([[:alpha:]]).*", replacement= "", x=colnames(reads))
  reads2=reads[,-1]
  rownames(reads2)=reads[,1]
  
  final_reads =reads2 %>% 
    relocate(order(colnames(.)))
  final_reads=final_reads[,c(1:37)]
  
  ##read in and format our metadata file
  meta=read.csv("DESEQMeta.csv", header=TRUE, sep=",")
  meta=meta[order(meta$Sample),]
  ##center and scale symbiont density##
  meta$Symbiont_density=scale(meta$Symbiont_density, center=TRUE, scale=TRUE)
  meta$SQRT_Symbiont=scale(meta$SQRT_Symbiont, center=TRUE, scale=TRUE)
  colnames(meta)[5]= "Symbiont_density"
  colnames(meta)[6]= "SQRT_Symbiont"
  write.csv(meta,"final_meta.csv")

##filter based on mean minimum##
  means <- apply(final_reads,1,mean)
  table(means>3)
  
  means3 <- names(means[means>3])
  head(means3)
  length(means3)
  
  countFilt <- final_reads[row.names(final_reads) %in% means3,]
  head(countFilt)
  
  totalCountsFilt <- colSums(countFilt)
  totalCountsFilt
  
  ##figures out read count distributions##
  min(totalCountsFilt) #213393
  max(totalCountsFilt) #659619
  mean(totalCountsFilt) #424527.2

##NORMALIZING READS 
  dds <- DESeqDataSetFromMatrix(countData = countFilt, 
                              #colData = meta, 
                              #design = ~ Combo)

  dds <- estimateSizeFactors(dds)         
  dds <- estimateDispersions(dds)

vst <- getVarianceStabilizedData(dds)
write.csv(vst, file = "normalizedreads_filtered.csv")

##Now modeling for treament* symbiont density 
  dds2 <- model.matrix( ~ Treatment*Symbiont_density, meta)
  object2 = DESeqDataSetFromMatrix(countData= countFilt, colData= meta, design = ~ Treatment*Symbiont_density , ignoreRank= TRUE )
  
  object2 <- estimateSizeFactors(object2)
  
  object2 <- estimateDispersions(object2, modelMatrix = dds2)
  dds3 <- nbinomWaldTest(object2, maxit=500, modelMatrix = dds2)  
  
  resultsNames(dds3)
  ##Let's look at genes that are differentially expressed due to interaction##
  Int2=results(dds3,contrast=list(c("Treatmenttreatment.Symbiont_density")))
  Int2ordered = Int2[order(Int2$padj),]
  head(Int2ordered, 70) ##69
  write.csv(Int2ordered, file = "DESeq_Count_Int_filtered.csv")

##lets do an additive model to look at the main effects (treatment and symbiont density) independently, controlling for the others##
  dds4 <- model.matrix( ~ Treatment+Symbiont_density, meta)
  object4 = DESeqDataSetFromMatrix(countData= countFilt, colData= meta, design = ~ Treatment+Symbiont_density , ignoreRank= TRUE )
  
  object4 <- estimateSizeFactors(object4)
  
  object4 <- estimateDispersions(object4, modelMatrix = dds4)
  dds5 <- nbinomWaldTest(object4, maxit=500, modelMatrix = dds4)  
  
  resultsNames(dds5)
  
  ##Let's look at genes that are differentially expressed between browns and whites##
  Count=results(dds5,contrast=list(c("Symbiont_density")))
  Countordered = Count[order(Count$padj),]
  head(Countordered, 75) #71
  write.csv(Countordered, file = "DESeq_Count_Additive_Symbiont_filtered.csv")
  
  ##Let's look at genes that are differentially expressed due to interaction##
  Int2=results(dds5,contrast=list(c("Treatmenttreatment")))
  Int2ordered = Int2[order(Int2$padj),]
  head(Int2ordered, 275) ##271
  write.csv(Int2ordered, file = "DESeq_Count_Additive_Treatment_filtered.csv")
