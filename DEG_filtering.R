##code for selecting candidate immune genes from DESeq results from Changsut et al. 2023##
##last edited by LEF 4/20/2024##

library(tidyr)
library(data.table)
library(dplyr)

##format our uniprot info##
annos=read.delim("astrangia_annos_clean.txt", sep="\t")
annos=annos[,c(1:3,5)]
colnames(annos)[4]="GO"
colnames(annos)[3]="name"

##Int DEGs
  Int=read.csv("DESeq_Count_Int_filtered.csv", header=TRUE)
  Int=Int[,c(1,3,7)]
  Int_DEGs <- Int %>% filter(padj < 0.1)
  colnames(Int_DEGs)[1]="transcript"
  Int_DEGs_annos=merge(Int_DEGs,annos,by="transcript") ##45 out of 69 are annotated##
  
  ##find those with GOs involved in immunity##
  immunity <- Int_DEGs_annos[Int_DEGs_annos$GO %like% "immun",]
  defense <- Int_DEGs_annos[Int_DEGs_annos$GO %like% "defense",]
  inflamm <- Int_DEGs_annos[Int_DEGs_annos$GO %like% "inflamm",]
  phagocy <- Int_DEGs_annos[Int_DEGs_annos$GO %like% "phagocy",]
  
  Int_immune=rbind(immunity,defense,inflamm,phagocy)
  Int_immune=distinct(Int_immune)
  write.csv(Int_immune, "Int_DEGs_immune.csv", row.names=FALSE)
  
##Symb DEGs
  Symb=read.csv("DESeq_Count_Additive_Symbiont_filtered.csv", header=TRUE)
  Symb=Symb[,c(1,3,7)]
  Symb_DEGs <- Symb %>% filter(padj < 0.1)
  colnames(Symb_DEGs)[1]="transcript"
  Symb_DEGs_annos=merge(Symb_DEGs,annos,by="transcript") ##49 out of 71 are annotated##
  
  ##find those with GOs involved in immunity##
  immunity <- Symb_DEGs_annos[Symb_DEGs_annos$GO %like% "immun",]
  defense <- Symb_DEGs_annos[Symb_DEGs_annos$GO %like% "defense",]
  inflamm <- Symb_DEGs_annos[Symb_DEGs_annos$GO %like% "inflamm",]
  phagocy <- Symb_DEGs_annos[Symb_DEGs_annos$GO %like% "phagocy",]
  
  Symb_immune=rbind(immunity,defense,inflamm,phagocy)
  Symb_immune=distinct(Symb_immune)
  write.csv(Symb_immune, "Symb_DEGs_immune.csv", row.names=FALSE)
  
##Treat DEGs
  Treat=read.csv("DESeq_Count_Additive_Treatment_filtered.csv", header=TRUE)
  Treat=Treat[,c(1,3,7)]
  Treat_DEGs <- Treat %>% filter(padj < 0.1)
  colnames(Treat_DEGs)[1]="transcript"
  Treat_DEGs_annos=merge(Treat_DEGs,annos,by="transcript") ##211 out of 271 are annotated##
  
  ##find those with GOs involved in immunity##
  immunity <- Treat_DEGs_annos[Treat_DEGs_annos$GO %like% "immun",]
  defense <- Treat_DEGs_annos[Treat_DEGs_annos$GO %like% "defense",]
  inflamm <- Treat_DEGs_annos[Treat_DEGs_annos$GO %like% "inflamm",]
  phagocy <- Treat_DEGs_annos[Treat_DEGs_annos$GO %like% "phagocy",]
  
  Treat_immune=rbind(immunity,defense,inflamm,phagocy)
  Treat_immune=distinct(Treat_immune)
  write.csv(Treat_immune, "Treat_DEGs_immune.csv", row.names=FALSE)
