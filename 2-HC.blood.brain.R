setwd("~/gxue/PPMI/RNA-seq/PPMI-IR3/IR3_Analysis/")

library(tidyverse)
library(data.table)

library(ggplot2)
library(ggpointdensity) 
library(ggpubr)
library(Rmisc)

# load transcriptome data of PPMI
BL.TPM<-fread("~/gxue/PPMI/RNA-seq/PPMI-IR3/merge/TPM/IR3.BL.TPM.txt",sep ="\t",header = F,stringsAsFactors = F) 
  colnames(BL.TPM)<-c("FileName","gene_id","TPM")
  
#gene expression matrix transformed 
BL.TPM<- BL.TPM %>% tidyr::separate(FileName,into=c(NA,"PATNO"),sep="\\.",extra = "drop")%>%
    data.table::dcast(formula = gene_id~PATNO,value.var = "TPM") 

phen.Original<-read.csv("~/gxue/PPMI/Phenotype/PPMI_Original_Cohort_BL_to_Year_5_Dataset_Apr2020.csv",
                        na.strings = c("NA",""),stringsAsFactors = F)

#APPRDX
#1	PD Participant
#2	Healthy Control
#4	Prodromal 

BL.HC.phen<-phen.Original %>% filter(EVENT_ID=="BL") %>%
  filter(APPRDX==2) %>% dplyr::mutate(PATNO=as.character(PATNO))

# median value of expression of PD in blood
BL.HC.TPM<-BL.TPM %>% column_to_rownames("gene_id") %>% 
  dplyr::select(intersect(colnames(BL.TPM),BL.HC.phen$PATNO)) %>% 
  apply(1,median) %>% as.data.frame() %>% rownames_to_column("gene_id") %>% 
  tidyr::separate(gene_id,into=c("gene_id",NA),sep="\\.")
colnames(BL.HC.TPM)[2]<-"Blood"


# load median value of expression matrix of brain tissues from GTEx
GTEx <- fread("./Brain-blood/GTEx/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct",
              header = T,sep = "\t",stringsAsFactors = F)

Brain.GTEx<-GTEx %>% dplyr::select(starts_with("Brain"),Name) %>% 
  tidyr::separate(Name,into=c("gene_id",NA),sep="\\.") 

Blood_Brain<-inner_join(BL.HC.TPM,Brain.GTEx,by="gene_id") %>% dplyr::select(-gene_id)

#Comparison of gene expression levels between blood and each brain tissue
p<-list()

for(i in 2:ncol(Blood_Brain)){
  
  #keep the genes with TPM > 0.1
  df<-data.frame(Blood=Blood_Brain[,1],Brain=Blood_Brain[,i]) %>%
    filter(Blood>0.1 & Brain>0.1)
  
  df.log2<-apply(df,2,log2) %>% as.data.frame()
  
  p[[i-1]]<-ggplot(data = df.log2, mapping = aes(x = Blood, y = Brain)) +
    ylab(colnames(Blood_Brain[i]))+
    geom_pointdensity(adjust = 10,size = 0.1) +
    scale_colour_gradient(low = "#1874CD", high = "#D95F02") +
    stat_cor(method = "pearson", label.x = -3, label.y = 15,size = 5)+
    theme_bw()+
    theme(panel.grid = element_blank(),panel.border = element_rect(size=1))+
    theme(axis.ticks.length=unit(.25, "cm"))+
    theme(axis.ticks = element_line(size=0.8),
          axis.text = element_text(size=14,colour = "black",angle = 00),
          axis.title = element_text(size=14,color="black"),
          legend.title=element_text(size=12),legend.text=element_text(size=12))
}



pdf("./Brain-blood/GTEx/GTEx_blood_brain1-4.pdf",width = 10.5,height = 7)
multiplot(plotlist = p[1:4] , cols = 2)
dev.off()

pdf("./Brain-blood/GTEx/GTEx_blood_brain5-8.pdf",width = 10.5,height = 7)
multiplot(plotlist = p[5:8] , cols = 2)
dev.off()


pdf("./Brain-blood/GTEx/GTEx_blood_brain9.pdf",width = 5,height = 3.5)
p[9]
dev.off()

pdf("./Brain-blood/GTEx/GTEx_blood_brain10.pdf",width = 5,height = 4.5)
p[10]
dev.off()

pdf("./Brain-blood/GTEx/GTEx_blood_brain11-12.pdf",width = 10.5,height = 3.5)
multiplot(plotlist = p[11:12] , cols = 2)
dev.off()


pdf("./Brain-blood/GTEx/GTEx_blood_brain13.pdf",width = 5,height = 3.5)
p[13]
dev.off()

save.image("~/gxue/PPMI/RNA-seq/PPMI-IR3/IR3_Analysis/Brain-blood/GTEx/HC.blood.brain.RData")
