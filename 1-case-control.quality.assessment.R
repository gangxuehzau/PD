setwd("~/gxue/PPMI/RNA-seq/PPMI-IR3/IR3_Analysis/")
library(data.table)
library(tidyverse)
library(table1)
library(uwot)
library(Rmisc)
library(DESeq2)
library(BiocParallel)
library(ggpubr)

# load transcriptome data 
BL.RNA.count<-fread("~/gxue/PPMI/RNA-seq/PPMI-IR3/merge/counts/IR3.BL.counts.txt",
                    sep ="\t",header = F,stringsAsFactors = F)
colnames(BL.RNA.count)<-c("FileName","gene_id","Count")

#add bash informations
Batch<-BL.RNA.count[,1] %>% as.data.frame() %>% distinct() %>%
  tidyr::separate(FileName,
                  into=c("version","PATNO","visit","sampleID","sequencing_facility_ID","assay"),
                  sep="\\.",extra = "drop")

#gene expression matrix transformed  
BL.RNA.count <- BL.RNA.count %>% tidyr::separate(FileName,into=c(NA,"PATNO"),sep="\\.",extra = "drop")%>%
  data.table::dcast(formula = gene_id~PATNO,value.var = "Count") %>% 
  tibble::column_to_rownames("gene_id")

save(Batch,BL.RNA.count,file = "./Baseline_QA/BL.RNA.count.RData")

# load demographics and characteristics of study cohort
phen.Original<-read.csv("~/gxue/PPMI/Phenotype/PPMI_Original_Cohort_BL_to_Year_5_Dataset_Apr2020.csv",
                        na.strings = c("NA",""),stringsAsFactors = F)
Prodromal.phen<-read.csv("~/gxue/PPMI/Phenotype/PPMI_Prodromal_Cohort_BL_to_Year_1_Dataset_Apr2020.csv",
                         na.strings = c("NA",""),stringsAsFactors = F) 

sample.select<-function(phen,number,type){
  
  rna.phen<-phen %>% filter(EVENT_ID=="BL") %>%
    filter(APPRDX==number) %>%
    dplyr::mutate(PATNO=as.character(PATNO)) %>% 
    dplyr::select(PATNO,APPRDX,age_cat,age,gen,educ,NHY,PD_MED_USE,updrs_totscore) %>% 
    inner_join(Batch[,1:2],by="PATNO") %>%
    tidyr::separate(version,into=c(NA,"Phase"),sep="-",extra="drop") %>% 
    dplyr::mutate(condition=rep(type,nrow(.))) 
  
  rna<- BL.RNA.count[,match(rna.phen$PATNO,colnames(BL.RNA.count))]
  return(list(rna.phen,rna))
}

#APPRDX
#1	PD 
#2	Healthy Control
#4	Prodromal 

PD.rna<-sample.select(phen.Original,1,"PD")
HC.rna<-sample.select(phen.Original,2,"HC")
Prodromal.rna<-sample.select(Prodromal.phen,4,"Prodromal")

HC.Prodromal.PD.rna<-cbind(HC.rna[[2]],Prodromal.rna[[2]],PD.rna[[2]])
HC.Prodromal.PD.phen<- data.frame(row.names=colnames(HC.Prodromal.PD.rna),
                                  rbind(HC.rna[[1]],Prodromal.rna[[1]],PD.rna[[1]]) %>% 
                                    dplyr::select(-c("PATNO","APPRDX")) %>%
                                    mutate(Treatment=ifelse(PD_MED_USE==0,"Unmedicated","Medicated")) %>% 
                                    mutate_at(c("age_cat","gen","educ","Phase","condition",
                                                "Treatment","NHY"),as.factor))


# plot the demographic and characteristics table

HC.Prodromal.PD.phen$age_cat<-factor(HC.Prodromal.PD.phen$age_cat, 
                                     levels=c("1","2","3"),
                                     labels=c("<56","56-65",">65"))

HC.Prodromal.PD.phen$gen<-factor(HC.Prodromal.PD.phen$gen,
                                 levels=c("1","2"),
                                 labels=c("Male","Female"))

HC.Prodromal.PD.phen$educ<-factor(HC.Prodromal.PD.phen$educ, 
                                  levels=c("1","2","3"),
                                  labels=c("<13","13-23",">23"))

label(HC.Prodromal.PD.phen$gen) <- "Gender"
label(HC.Prodromal.PD.phen$age_cat) <- "Age Category"
label(HC.Prodromal.PD.phen$age) <- "Age"
label(HC.Prodromal.PD.phen$educ) <- "Education"
label(HC.Prodromal.PD.phen$NHY) <- "Hoehn & Yahr Stage (OFF)"
label(HC.Prodromal.PD.phen$updrs_totscore) <- "MDS-UPDRS Total Score (OFF)"
label(HC.Prodromal.PD.phen$Phase) <- "PPMI Project Phase"

units(HC.Prodromal.PD.phen$age)  <- "years"
units(HC.Prodromal.PD.phen$age_cat)  <- "years"
units(HC.Prodromal.PD.phen$educ)  <- "years"

table1(~ age + age_cat + gen + educ + Treatment + NHY +
         updrs_totscore + Phase| condition, HC.Prodromal.PD.phen,overall = F)


#Demographic and characteristics of sudy cohort
#were compared using student’s t-test for continuous variables,
#and fisher exact test for categorical variables.

rndr <- function(x, name, ...) {
  if (length(x) == 0) {
    y <- tmp[[name]]
    s <- rep("", length(render.default(x=y, name=name, ...)))
    if (is.numeric(y)) {
      p <- t.test(y ~ tmp$condition)$p.value
    } else {
      p <- fisher.test(table(y, droplevels(tmp$condition)))$p.value
    }
    s[2] <- sub("<", "&lt;", format.pval(p, digits=3, eps=0.001))
    s
  } else {
    render.default(x=x, name=name, ...)
  }
}

rndr.strat <- function(label, n, ...) {
  ifelse(n==0, label, render.strat.default(label, n, ...))
}

# PD vs Prodromal
tmp<-HC.Prodromal.PD.phen %>% filter(condition!="HC") %>% 
  mutate(condition=factor(condition, levels=c("HC", "PD", "Prodromal"),
                          labels=c("P-value", "PD", "Prodromal")))
table1(~ age + age_cat + gen  + educ + NHY + updrs_totscore + Phase|condition, 
       data=tmp, droplevels=F, render=rndr, render.strat=rndr.strat, overall=F)

# Prodromal vs HC
tmp<-HC.Prodromal.PD.phen %>% filter(condition!="PD") %>% 
  mutate(condition=factor(condition, levels=c("HC", "PD", "Prodromal"), 
                          labels=c("HC", "P-value", "Prodromal")))
table1(~ age + age_cat + gen  + educ + NHY + updrs_totscore + Phase|condition, 
       data=tmp, droplevels=F, render=rndr, render.strat=rndr.strat, overall=F)

# PD vs HC
tmp<-HC.Prodromal.PD.phen %>% filter(condition!="Prodromal") %>% 
  mutate(condition=factor(condition, levels=c("HC", "PD", "Prodromal"), 
                          labels=c("HC", "PD", "P-value")))
table1(~ age + age_cat + gen  + educ + NHY + updrs_totscore + Phase|condition, 
       data=tmp, droplevels=F, render=rndr, render.strat=rndr.strat, overall=F)



# normalize gene expression
dds <- DESeqDataSetFromMatrix(countData = HC.Prodromal.PD.rna,
                              colData = HC.Prodromal.PD.phen,
                              design = ~ condition)

# prepare for UMAP and heatmap
vsd <- vst(dds, blind=T)
vsd.expr<-assay(vsd)

# prepare for comparison with GSE68719 dataset
register(MulticoreParam(12))
dds2<- DESeq(dds,parallel=TRUE) 
blood.normalized.counts <- counts(dds2, normalized=TRUE)

# UMAP 
set.seed(1234)
RNA.umap = uwot::umap(t(vsd.expr),n_neighbors = 50,metric = "correlation")
colnames(RNA.umap)<-c("UMAP_1","UMAP_2")
RNA.umap<-cbind(RNA.umap,HC.Prodromal.PD.phen) 

# plot UMAP results
pdf("./Baseline_QA/umap_baseline.pdf",width = 9,height = 5)
p<-list()
p[[1]]<-ggplot(RNA.umap,aes(x=UMAP_1,y=UMAP_2,color=condition))+ geom_point() 
p[[2]]<-ggplot(RNA.umap,aes(x=UMAP_1,y=UMAP_2,color=age_cat))+ geom_point() 
p[[3]]<-ggplot(RNA.umap,aes(x=UMAP_1,y=UMAP_2,color=gen))+ geom_point() 
p[[4]]<-ggplot(RNA.umap,aes(x=UMAP_1,y=UMAP_2,color=Phase))+ geom_point() 
Rmisc::multiplot(plotlist = p[1:4], cols=2)
dev.off()

# save data
saveRDS(HC.rna,file = "./Baseline_QA/Baseline_HC.rna.phen.rds")
saveRDS(Prodromal.rna,file = "./Baseline_QA/Baseline_Prodromal.rna.phen.rds")
saveRDS(PD.rna,file = "./Baseline_QA//Baseline_PD.rna.phen.rds")
saveRDS(blood.normalized.counts,file = "./Baseline_QA/blood.normalized.counts.rds")
saveRDS(vsd.expr,file = "./Baseline_QA/Baseline_HC.Prodromal.PD.rna.rds")


save.image("./Baseline_QA/Baseline-RNA-seq.QA.RData")

