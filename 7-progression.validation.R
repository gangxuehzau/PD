setwd("~/gxue/PPMI/RNA-seq/PPMI-IR3/IR3_Analysis/")
library(tidyverse)
library(pheatmap)
library(RColorBrewer)
library(edgeR)
library(limma)

# load data
load("./longitudinal_QA/longitudinal.raw.RData")
load("./longitudinal_QA/overlap.PD.RData")
phen.Original<-read.csv("~/gxue/PPMI/Phenotype/PPMI_Original_Cohort_BL_to_Year_5_Dataset_Apr2020.csv",
                        na.strings = c("NA",""),stringsAsFactors = F)
progression.genes<-readRDS("./longitudinal_analysis/progression.genes.rds")

# select independent PD patients
# residue(res)

PD.Y<-function(Y,expr){
  
  phen<- phen.Original %>% filter(APPRDX==1) %>% 
    dplyr::select(PATNO,updrs_totscore_on,PD_MED_USE,YEAR,age,gen,EVENT_ID) %>% 
    tidyr::unite("PATNO_visit",PATNO,EVENT_ID,remove = F) %>% 
    filter(YEAR==Y)
  phen$updrs_totscore_on[is.na(phen$updrs_totscore_on)]<-median(phen$updrs_totscore_on, na.rm = T)
  
  phen.res<-phen[!(phen$PATNO %in% overlap.PD),]
  phen.res<-phen.res[phen.res$PATNO %in% colnames(expr),]
  phen.res$PATNO<-as.character(phen.res$PATNO)
  return(phen.res)
}

PD.BL.res<-PD.Y(0,BL[[1]])
PD.Y1.res<-PD.Y(1,V04[[1]]) 
PD.Y2.res<-PD.Y(2,V06[[1]])
PD.Y3.res<-PD.Y(3,V08[[1]])

RNA.res<-cbind(BL[[1]] %>% dplyr::select(all_of(PD.BL.res$PATNO)),
                        V04[[1]] %>% dplyr::select(all_of(PD.Y1.res$PATNO)),
                        V06[[1]] %>% dplyr::select(all_of(PD.Y2.res$PATNO)),
                        V08[[1]] %>% dplyr::select(all_of(PD.Y3.res$PATNO)))

phen.res<-rbind(PD.BL.res,PD.Y1.res,PD.Y2.res,PD.Y3.res)

identical(phen.res$PATNO,colnames(RNA.res))
colnames(RNA.res)<-phen.res$PATNO_visit

rm(BL,V04,V06,V08,PD.BL.res,PD.Y1.res,PD.Y2.res,PD.Y3.res)

#estimate gene library sizes, remove the low expression genes 
#and estimate normalization factors
group<-as.factor(phen.res$YEAR)
y <- DGEList(counts=RNA.res, group=group)
keep <- filterByExpr(y,group = group)
y <- y[keep, , keep.lib.sizes=FALSE]
y <- calcNormFactors(y,method = "TMM")

design <- model.matrix(~0+group)
colnames(design)=levels(group)
rownames(design)=colnames(RNA.res)
#voom converts raw counts to log-CPM values by
#automatically extracting library sizes and normalisation factors
y <- voom(y, design,plot = T)
RNA.res.lcpm<-y$E %>% as.data.frame()

#heatmap of progression genes in independent PD patients
set3.col<-brewer.pal(9,"Set3")
phen.res.sort<-arrange(phen.res,updrs_totscore_on) %>% 
  mutate_at(c("YEAR","gen","PD_MED_USE"),as.factor)

phen.res.sort$PD_MED_USE<-factor(phen.res.sort$PD_MED_USE,
                                          levels=c("0","1","2","3","4","5","6","7"),
                                          labels=c("Unmedicated","Levodopa","DopamineAgonist",
                                                   "Other","Levodopa_Other","Levodopa_DopamineAgonist",
                                                   "DopamineAgonist_Other","Levodopa_DopamineAgonist_Other"))

phen.res.sort$YEAR<-factor(phen.res.sort$YEAR,
                                    levels=c("0","1","2","3"),
                                    labels=c("Baseline","Year1","Year2","Year3"))

phen.res.sort$gen<-factor(phen.res.sort$gen,
                                   levels=c("1","2"),
                                   labels=c("Male","Female"))


RNA.res.progression<-RNA.res.lcpm[rownames(RNA.res.lcpm) %in% progression.genes$gene_id,]
RNA.res.progression<-RNA.res.progression[,match(phen.res.sort$PATNO_visit,colnames(RNA.res.progression))]

annotation_col = data.frame(
  MDS_UPDRS_total_score_on = phen.res.sort$updrs_totscore_on,
  Treatment = phen.res.sort$PD_MED_USE,
  Time = phen.res.sort$YEAR,
  Age = phen.res.sort$age,
  Gender=phen.res.sort$gen
)

rownames(annotation_col) = colnames(RNA.res.progression)

ann_colors = list(
  MDS_UPDRS_total_score_on= c("white","firebrick"),
  Age = c("#F0FFFF","#FA8072"),
  Time = c(Baseline="#5F9EA0",Year1="#FFDAB9",Year2="#F4A460",Year3="#FA8072"),
  Gender = c(Male= "#5F9EA0", Female = "#FF8C00"),
  Treatment = c(Unmedicated="#D9D9D9",Levodopa="#FCCDE5",DopamineAgonist="#FFFFB3",
                Other="#BEBADA",Levodopa_Other="#FB8072",Levodopa_DopamineAgonist="#80B1D3",
                DopamineAgonist_Other="#FDB462",Levodopa_DopamineAgonist_Other="#B3DE69")
)

pdf("./longitudinal_analysis/independent.heatmap.genes.pdf",width = 9,height = 5.7)
pheatmap(RNA.res.progression, 
         annotation_col = annotation_col,
         scale = "row",
         cluster_cols = F,show_colnames = F,show_rownames = F,
         treeheight_row = 25,
         fontsize = 9,
         cutree_rows = 2,
         border=T,
         annotation_colors = ann_colors,
         cluster_rows = T,
         color = colorRampPalette(c("#1874CD", "white", "#FA8072"))(20),
         breaks = seq(-1,1.3,length.out = 20)
)

dev.off()

#save data
save.image("./longitudinal_analysis/independent.RData")
