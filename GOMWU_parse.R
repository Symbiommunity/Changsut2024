##code for parsing DESeq results to pass to GOMWU from Changsut et al. 2023##
##last edited by LEF 4/20/2024##

  library(tidyr)
  ##parse GOMWU results##
  infect=read.csv("DESeq_Count_Int_filtered.csv")
  infect=infect[,c(1,3,7)]
  ##adjust pvalues##
  infect$logpadj=-log(infect$padj,10)
  sign <- rep(1,nrow(infect))
  sign[infect$log2FoldChange<0]=-1
  table(sign)
  infect$logpadj=infect$logpadj*sign
  infect=infect[,c(1,4)]
  ##write it out##
  write.csv(infect, "Count_Int_filtered_GOMWU.csv", row.names = FALSE)
  
  ##parse GOMWU results##
  infect=read.csv("DESeq_Count_Additive_Symbiont_filtered.csv")
  infect=infect[,c(1,3,7)]
  ##adjust pvalues##
  infect$logpadj=-log(infect$padj,10)
  sign <- rep(1,nrow(infect))
  sign[infect$log2FoldChange<0]=-1
  table(sign)
  infect$logpadj=infect$logpadj*sign
  infect=infect[,c(1,4)]
  ##write it out##
  write.csv(infect, "Count_Symbiont_filtered_GOMWU.csv", row.names = FALSE)
  
  ##parse GOMWU results##
  infect=read.csv("DESeq_Count_Additive_Treatment_filtered.csv")
  infect=infect[,c(1,3,7)]
  ##adjust pvalues##
  infect$logpadj=-log(infect$padj,10)
  sign <- rep(1,nrow(infect))
  sign[infect$log2FoldChange<0]=-1
  table(sign)
  infect$logpadj=infect$logpadj*sign
  infect=infect[,c(1,4)]
  ##write it out##
  write.csv(infect, "Count_Treatment_filtered_GOMWU.csv", row.names = FALSE)
  
  ##make annos table##
  annos=read.delim("astrangia_annos_clean.txt", sep="\t")
  write.table(annos,"annos.tab", sep = "\t",row.names=FALSE)
