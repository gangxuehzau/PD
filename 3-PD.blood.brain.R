setwd("~/gxue/PPMI/RNA-seq/PPMI-IR3/IR3_Analysis/")
library(tidyverse)
library(data.table)
library(ggpointdensity) 
library(ggplot2)
library(ggpubr)

#load GSE68719 dataset: mRNA-Seq profiles normalized by DESeq2 of human post-mortem BA9 brain tissue of PD patients (N = 29)
expr.brain<-fread("./Brain-blood/GSE68719_mlpd_PCG_DESeq2_norm_counts.txt")
#load blood RNA-seq(PPMI) normalized by DESeq2
blood.normalized_counts<-readRDS("./Baseline_QA/blood.normalized.counts.rds")

PD.phen<-read.csv("~/gxue/PPMI/Phenotype/PPMI_Original_Cohort_BL_to_Year_5_Dataset_Apr2020.csv",
                        na.strings = c("NA",""),stringsAsFactors = F) %>% 
  filter(EVENT_ID=="BL") %>%
  filter(APPRDX==1) %>% dplyr::mutate(PATNO=as.character(PATNO))

# median value of expression
PD.blood<-blood.normalized_counts[,colnames(blood.normalized_counts) %in% PD.phen$PATNO] %>%
  apply(2,log2) %>% apply(1,median) %>% as.data.frame() %>% 
  rownames_to_column("gene_id") %>% tidyr::separate(gene_id,into=c("gene_id",NA),sep="\\.")
colnames(PD.blood)[2]<-"median_expr_blood"

PD.brain<- expr.brain %>% column_to_rownames("EnsemblID") %>% 
  dplyr::select(starts_with("P")) %>% apply(2,log2) %>% apply(1,median) %>% 
  as.data.frame() %>% rownames_to_column("EnsemblID") %>% 
  tidyr::separate(EnsemblID,into=c("gene_id",NA),sep="\\.")
colnames(PD.brain)[2]<-"median_expr_brian"

#keep the genes with counts >= 10 
PD.brain.blood<-inner_join(PD.brain,PD.blood,by="gene_id") %>% 
  filter(median_expr_blood>=log2(10) & median_expr_brian>=log2(10))

#Comparison of gene expression levels between blood and brain tissue
pdf("./Brain-blood/genome_wide_correlation_PD_blood_brain.pdf",width = 5.5,height = 4.5)
ggplot(data = PD.brain.blood, mapping = aes(x = median_expr_blood, y = median_expr_brian)) +
  geom_pointdensity(adjust = 10,size = 0.1) +
  scale_colour_gradient(low = "#1874CD", high = "#D95F02") +
  stat_cor(method = "pearson", label.x = 4, label.y = 19,size = 5.5)+
  ylab("Median log2(Normalized counts) of brain")+
  xlab("Median log2(Normalized counts) of blood")+
  theme_bw()+
  theme(panel.grid = element_blank(),panel.border = element_rect(size=1))+
  theme(axis.ticks.length=unit(.25, "cm"))+
  theme(axis.ticks = element_line(size=0.8),
        axis.text = element_text(size=14,colour = "black",angle = 00),
        axis.title = element_text(size=14,color="black"),
        legend.title=element_text(size=12),legend.text=element_text(size=12))
dev.off()

# save data
save.image("./Brain-blood/PD.blood.brain.RData")
