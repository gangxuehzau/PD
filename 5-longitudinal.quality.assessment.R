setwd("~/gxue/PPMI/RNA-seq/PPMI-IR3/IR3_Analysis/")

library(tidyverse)
library(data.table)
library(edgeR)
library(limma)
library(umap)
library(ggplot2)
library(ggpubr)
library(table1)
library(Rmisc)

# load transcriptome data 
get.rna<-function(data){
  
  RNA.count<-fread(data, sep ="\t",header = F,stringsAsFactors = F)
  colnames(RNA.count)<-c("FileName","gene_id","Count")
  
  #add bash informations
  Batch<-RNA.count[,1] %>% as.data.frame() %>% distinct() %>%
    tidyr::separate(FileName,
                    into=c("version","PATNO","visit","sampleID","sequencing_facility_ID","assay"),
                    sep="\\.",extra = "drop")
  
  #gene expression matrix transformed 
  RNA.count <- RNA.count %>% tidyr::separate(FileName,into=c(NA,"PATNO"),sep="\\.",extra = "drop")%>%
    data.table::dcast(formula = gene_id~PATNO,value.var = "Count") %>% 
    tibble::column_to_rownames("gene_id")
  
  return(list(RNA.count,Batch))
}

# BL: Baseline
# V04:  1 year 
# V06: 2 years
# V08: 3 years
BL<-get.rna("../merge/counts/IR3.BL.counts.txt")
V04<-get.rna("../merge/counts/IR3.V04.counts.txt")
V06<-get.rna("../merge/counts/IR3.V06.counts.txt")
V08<-get.rna("../merge/counts/IR3.V08.counts.txt")
save(BL,V04,V06,V08,file = "./longitudinal_QA/longitudinal.raw.RData")

#load demographics and characteristics of study cohort
phen.Original<-read.csv("~/gxue/PPMI/Phenotype/PPMI_Original_Cohort_BL_to_Year_5_Dataset_Apr2020.csv",
                        na.strings = c("NA",""),stringsAsFactors = F)

#select PD samples that have longitudinal (baseline, 1 year, 2 years and 3 years) RNA-seq and clinial assessment (phenotype)
phen.PD<-phen.Original %>% filter(APPRDX==1) %>%
  dplyr::select(PATNO,APPRDX,EVENT_ID,YEAR,age,age_cat,gen,educ,NHY,NHY_ON,PD_MED_USE,updrs_totscore,updrs_totscore_on)%>%
  dplyr::mutate(PATNO=as.character(PATNO)) %>%
  filter(EVENT_ID=="BL" | EVENT_ID=="V04" | EVENT_ID=="V06" | EVENT_ID=="V08")

phen.PD<-phen.PD %>% dplyr::count(PATNO) %>% filter(n==4) %>%
  inner_join(phen.PD,by="PATNO") %>%
  dplyr::select(-n)

# overlap phenotype and transcriptome
overlap.PD<-Reduce(intersect,  list(v1 = colnames(BL[[1]]),
                                    v2 = colnames(V04[[1]]),
                                    v3 = colnames(V06[[1]]),
                                    v4 = colnames(V08[[1]]),
                                    v5 = unique(phen.PD$PATNO)))


RNA.longitudinal<-cbind(BL[[1]] %>% dplyr::select(all_of(overlap.PD)),
                       V04[[1]] %>% dplyr::select(all_of(overlap.PD)),
                       V06[[1]] %>% dplyr::select(all_of(overlap.PD)),
                       V08[[1]] %>% dplyr::select(all_of(overlap.PD)))

phen.longitudinal<-phen.PD[phen.PD$PATNO %in% overlap.PD,] %>% arrange(YEAR,PATNO) %>%
  tidyr::unite("PATNO_visit",PATNO,EVENT_ID,remove = F)


# add batch
batch<-rbind(BL[[2]],V04[[2]],V06[[2]],V08[[2]]) %>%
  tidyr::separate(version,into=c(NA,"Phase"),sep="-",extra="drop")
batch<-batch[batch$PATNO%in%overlap.PD,]%>% arrange(visit,PATNO)

#identical(batch$PATNO,phen.longitudinal$PATNO)
phen.longitudinal<- mutate(phen.longitudinal,Phase=as.factor(batch$Phase))

phen.longitudinal<- phen.longitudinal %>%
  mutate_at(c("PATNO","YEAR","age_cat","gen","educ","NHY","NHY_ON","PD_MED_USE"),as.factor)

#identical(colnames(RNA.longitudinal),as.character(phen.longitudinal$PATNO))
colnames(RNA.longitudinal)<-phen.longitudinal$PATNO_visit

rm(phen.Original,phen.PD,batch,BL,V04,V06,V08)


# plot the demographic and characteristics table
phen.table<-phen.longitudinal

phen.table$PD_MED_USE<-factor(phen.table$PD_MED_USE,
                              levels=c("0","1","2","3","4","5","6","7"),
                              labels=c("Unmedicated","Levodopa","Dopamine Agonist",
                                       "Other","Levodopa + Other","Levodopa + Dopamine Agonist",
                                       "Dopamine Agonist + Other","Levodopa + Dopamine Agonist + Other"))

phen.table$YEAR<-factor(phen.table$YEAR,
                        levels=c("0","1","2","3"),
                        labels=c("Baseline","Year 1","Year 2","Year 3"))

phen.table$age_cat<-factor(phen.table$age_cat,
                           levels=c("1","2","3"),
                           labels=c("<56","56-65",">65"))

phen.table$gen<-factor(phen.table$gen,
                       levels=c("1","2"),
                       labels=c("Male","Female"))

phen.table$educ<-factor(phen.table$educ,
                        levels=c("1","2","3"),
                        labels=c("<13","13-23",">23"))

label(phen.table$YEAR) <- "Vist"
label(phen.table$gen) <- "Gender"
label(phen.table$age_cat) <- "Age Category"
label(phen.table$age) <- "Age"
label(phen.table$educ) <- "Education"
label(phen.table$PD_MED_USE)<-"Medications"
label(phen.table$NHY) <- "Hoehn & Yahr Stage (OFF)" 
label(phen.table$NHY_ON) <- "Hoehn & Yahr Stage (ON)"
label(phen.table$updrs_totscore) <- "MDS-UPDRS Total Score (OFF)"
label(phen.table$updrs_totscore_on) <- "MDS-UPDRS Total Score (ON)"
label(phen.table$Phase) <- "PPMI Project Phase"

units(phen.table$YEAR)  <- "years"
units(phen.table$age)  <- "years"
units(phen.table$age_cat)  <- "years"
units(phen.table$educ)  <- "years"

table1(~age + age_cat + gen + PD_MED_USE + NHY + NHY_ON + updrs_totscore + updrs_totscore_on + Phase| YEAR, phen.table,overall = F)


#estimate gene library sizes, remove the low expression genes 
#and estimate normalization factors
group<-phen.longitudinal$YEAR
y <- DGEList(counts=RNA.longitudinal, group=group)
keep <- filterByExpr(y,group = group)
y <- y[keep, , keep.lib.sizes=FALSE]
y <- calcNormFactors(y,method = "TMM")


design <- model.matrix(~0+group)
colnames(design)=levels(group)
rownames(design)=colnames(RNA.longitudinal)
#voom converts raw counts to log-CPM values by
#automatically extracting library sizes and normalisation factors
y <- voom(y, design,plot = T)
RNA.longitudinal.lcpm<-y$E

dim(RNA.longitudinal.lcpm)
# 26207  1020

#UMAP analysis
set.seed(1234)
RNA.umap = uwot::umap(t(RNA.longitudinal.lcpm),n_neighbors = 50,metric = "correlation")
colnames(RNA.umap)<-c("UMAP_1","UMAP_2")
phen.table.umap<-cbind(phen.table,RNA.umap)

# plot UMAP result
pdf("./longitudinal_QA/umap_longitudinal.pdf",width=9, height=5.5)
p<-list()
p[[1]]<-ggplot(phen.table.umap,aes(x=UMAP_1,y=UMAP_2,color=age_cat))+ geom_point()
p[[2]]<-ggplot(phen.table.umap,aes(x=UMAP_1,y=UMAP_2,color=gen))+ geom_point()
p[[3]]<-ggplot(phen.table.umap,aes(x=UMAP_1,y=UMAP_2,color=Phase))+ geom_point()
p[[4]]<-ggplot(phen.table.umap,aes(x=UMAP_1,y=UMAP_2,color=YEAR))+ geom_point()
Rmisc::multiplot(plotlist = p[1:4], cols=2)
dev.off()

# save data
save(overlap.PD,file = "./longitudinal_QA/overlap.PD.RData")
save(phen.longitudinal,RNA.longitudinal.lcpm,file = "./longitudinal_QA/simplified.longitudinal.QA.RData")
save.image("./longitudinal_QA/longitudinal.QA.RData")
