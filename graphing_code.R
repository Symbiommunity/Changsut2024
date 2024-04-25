##code for making graphs with Changsut et al. 2023##
##last edited by LEF 4/20/2024##

library(VennDiagram)
library(ggplot2)
library(ggpubr)
library(tidyr)

##make the graphs that are in the publication##
  reads=read.csv("normalizedreads_filtered.csv", header=TRUE)  ##normalizd reads##
  samples_overlap=read.csv("final_meta_graph.csv",header=TRUE) ##metadata for graphs##
##then mak figure 1- main effects of treatment
  GOI_Treat=read.csv("GOI_Treat.csv", header=TRUE) ##list of immune genes, which were output based on GO terms then hand verified for immune function##
  GOI_Treat_Gene=GOI_Treat$transcript
  
  GOI_Treat_matrix <- reads[reads$X %in% GOI_Treat_Gene, ]
  GOI_Treat_matrix = merge(GOI_Treat_matrix,GOI_Treat, by.x="X",by.y="transcript")
  GOI_Treat_matrix=GOI_Treat_matrix[,c(39,2:38)]
  GOI_Treat_matrix_t=t(GOI_Treat_matrix)                          
  colnames(GOI_Treat_matrix_t)<-GOI_Treat_matrix_t[1,]
  GOI_Treat_matrix_t<-GOI_Treat_matrix_t[-1,]     
  GOI_Treat_matrix_t=as.data.frame(GOI_Treat_matrix_t)
  GOI_Treat_matrix_final <- tibble::rownames_to_column(GOI_Treat_matrix_t, "Sample")
  
  ##merge in metadata##
  final_matrix=merge(GOI_Treat_matrix_final, samples_overlap, by="Sample")
  
  final_matrix_long <- gather(final_matrix, key="gene", value="expression", 
                              c("eef-2","Selp"))
  sapply(final_matrix_long$expression, class)
  final_matrix_long$expression = as.numeric(final_matrix_long$expression)
  
  treat=ggplot(final_matrix_long, aes(x = Treatment, y = expression, fill=Treatment))+
    geom_boxplot()+
    facet_wrap(.~gene) + theme_classic() +
    scale_fill_manual(values=c("#7790a9", "#14314e")) +
    xlab("Treatment") + ylab("Normalized Gene Expression") +
    theme_classic() + 
    theme(strip.background =element_rect(fill="gray")) 
  treat

##make figure 2- main effects of symbiont density##
  GOI_Symb=read.csv("GOI_Symb.csv", header=TRUE) ##list of immune genes, which were output based on GO terms then hand verified for immune function##
  GOI_Symb_Gene=GOI_Symb$transcript
  
  GOI_Symb_matrix <- reads[reads$X %in% GOI_Symb_Gene, ]
  GOI_Symb_matrix = merge(GOI_Symb_matrix,GOI_Symb, by.x="X",by.y="transcript")
  GOI_Symb_matrix=GOI_Symb_matrix[,c(39,2:38)]
  GOI_Symb_matrix_t=t(GOI_Symb_matrix)                          
  colnames(GOI_Symb_matrix_t)<-GOI_Symb_matrix_t[1,]
  GOI_Symb_matrix_t<-GOI_Symb_matrix_t[-1,]     
  GOI_Symb_matrix_t=as.data.frame(GOI_Symb_matrix_t)
  GOI_Symb_matrix_final <- tibble::rownames_to_column(GOI_Symb_matrix_t, "Sample")
  
  ##merge in metadata##
  final_matrix=merge(GOI_Symb_matrix_final, samples_overlap, by="Sample")
  
  final_matrix_long <- gather(final_matrix, key="gene", value="expression", 
                              c("eef-2","PRDX5","Selp","Xbp1"))
  sapply(final_matrix_long$expression, class)
  final_matrix_long$expression = as.numeric(final_matrix_long$expression)
  
  symbplot<- ggplot(final_matrix_long, aes(Symbiont_density, expression)) +
    geom_smooth(method = "lm", color = "grey25", se=FALSE) +
    geom_point(color = "black", size = 2.5, alpha = 0.75, aes(fill = Type, shape = Type)) + facet_wrap(~gene) + 
    scale_fill_manual(values = c("brown" = "burlywood", "white" = "gray"), name = "Symbiotic State") +  ylab("Normalized Gene Expression") + xlab("Symbiont Density \n (cells/mL)") +
    theme_classic() + 
    theme(strip.background =element_rect(fill="gray")) + stat_cor() + scale_shape_manual(values=c(21,23))
  symbplot
 

##Make Figure 3- symb*immune interactions
  GOI_Int=read.csv("GOI_Int.csv", header=TRUE) ##list of immune genes, which were output based on GO terms then hand verified for immune function#
  GOI_Int_Gene=GOI_Int$transcript
  
  GOI_Int_matrix <- reads[reads$X %in% GOI_Int_Gene, ]
  GOI_Int_matrix = merge(GOI_Int_matrix,GOI_Int, by.x="X",by.y="transcript")
  GOI_Int_matrix=GOI_Int_matrix[,c(39,2:38)]
  GOI_Int_matrix_t=t(GOI_Int_matrix)                          
  colnames(GOI_Int_matrix_t)<-GOI_Int_matrix_t[1,]
  GOI_Int_matrix_t<-GOI_Int_matrix_t[-1,]     
  GOI_Int_matrix_t=as.data.frame(GOI_Int_matrix_t)
  GOI_Int_matrix_final <- tibble::rownames_to_column(GOI_Int_matrix_t, "Sample")
  
  ##merge in metadata##
  final_matrix=merge(GOI_Int_matrix_final, samples_overlap, by="Sample")
  
  final_matrix_long <- gather(final_matrix, key="gene", value="expression", 
                              c("eef-2","Ifi44","Selp"))
  sapply(final_matrix_long$expression, class)
  final_matrix_long$expression = as.numeric(final_matrix_long$expression)
  
  ggplot(final_matrix_long, aes(Symbiont_density, expression)) +
    geom_point(color = "black", size = 2.5, alpha = 0.75, aes(fill = Treatment, shape = Type)) + facet_wrap(~gene) + 
    geom_smooth(data=final_matrix_long,
                aes(Symbiont_density,expression,color=factor(Treatment)),method=lm,se=FALSE) +
    scale_fill_manual(values = c("treatment" = "#14314e", "Control" = "#7790a9"), name = "Treatment")+
    scale_color_manual(values = c("treatment" = "#14314e", "Control" = "#7790a9"), name = "Treatment")+
    ylab("Gene Expression") + xlab("Symbiont Density \n (cells/mL)") + theme_classic() + 
    theme(strip.background =element_rect(fill="gray")) + scale_shape_manual(values=c(21,23))
  
##Figure 4- gene ontology of symbiont*immune interactions
  goinclude<- read.csv("GOs_Included.csv")
  
  dataforinteractiongraph <- goinclude[c(1,6,7)]
  specificdataint <- dataforinteractiongraph[c(1:23),]
  specificdataint$name<- factor(specificdataint$name, levels= c("inner dynein arm assembly","sperm axoneme assembly","cilium-dependent cell motility","cilium movement involved in cell motility","cilium movement","microtubule bundle formation","axonemal dynein complex assembly","microtubule-based movement","microtubule cytoskeleton organization","microtubule-based process","cell projection assembly","microtubule-based transport","cell projection organization","carbohydrate derivative catabolic process","macromolecule biosynthetic process","organonitrogen compound biosynthetic process","amide metabolic process","peptide metabolic process","aminoglycan metabolic process","amino sugar metabolic process","humoral immune response","aminoglycan catabolic process","amino sugar catabolic process"))
  intbar<- ggplot(specificdataint, aes(x= delta.rank, y= name), fill = color) + scale_fill_manual(values= c('white'="#00A08A", 'brown' ="#FF0000"))+ geom_bar(stat = "identity") + theme_bw() + theme(panel.grid = element_blank()) + xlab("Delta Rank") +ylab("Biological Process")
  intbar
  ggsave("intbarsicb.pdf", intbar, height = 4, width = 9)

##Supplemental Figure 1##
  Immune_GO_graph=read.csv("TOI_graph_supp1.csv") ##list of sig GO terms for main effect of immunity##
  ggplot(Immune_GO_graph, aes(x = reorder(name, -delta.rank), y = delta.rank)) + geom_col()+coord_flip()+
    theme_classic() +  ylab("Delta Rank") + xlab("Biological Process")

##Supplemental Figure 2##
  Symb_GO_graph=read.csv("TOI_graph_supp2.csv") ##list of sig GO terms for main effect of symbiont##
  ggplot(Symb_GO_graph, aes(x = reorder(name, -delta.rank), y = delta.rank)) + geom_col()+coord_flip()+
    theme_classic() +  ylab("Delta Rank") + xlab("Biological Process")

