setwd("/home/kyzhou/gxue/PPMI/RNA-seq/PPMI-IR3/IR3_Analysis/")

library(tidyverse)
library(data.table)
library(ez)
library(nlme)
library(Rmisc)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggpubr)
library(enrichplot)
library(ggthemes)
library(pheatmap)

# load data
load("./longitudinal_QA/simplified.longitudinal.QA.RData")
gene.anno<-fread("./annotation/gencode.v29.gene.anno.txt",header = T)

#find which gene expression change along with the time
run.ezANOVA.f<-function(x){
  
  gene<-phen.longitudinal %>% mutate(expr=as.numeric(t(x)))
  model<-ezANOVA(data = gene,
                 dv = expr,
                 wid=PATNO,
                 within=YEAR,
                 type = 2)
  ANOVA.p<-ifelse(as.numeric(model$`Mauchly's Test for Sphericity`[3])>0.05,
                  as.numeric(model$ANOVA$p),as.numeric(model$`Sphericity Corrections`$`p[GG]`))
  
  return(ANOVA.p)
  
}

gene_ezANOVA<-apply(RNA.longitudinal.lcpm,1,run.ezANOVA.f)
gene_ezANOVA<-gene_ezANOVA %>% as.data.frame() %>% rownames_to_column("gene_id")
colnames(gene_ezANOVA)[2]<-c("ANOVA_P")

gene_ezANOVA.filter<-gene_ezANOVA %>%
  dplyr::mutate(ANOVA_fdr=p.adjust(ANOVA_P,method="BH")) %>%
  filter(ANOVA_fdr<0.05) %>% arrange(ANOVA_fdr)

# Missing values of updrs_totscore_on were imputed using the median value.
complement.f<-function(Y){
  
  df <- filter(phen.longitudinal,YEAR==Y) %>%
    dplyr::select("PATNO_visit","PATNO","age","educ","gen","PD_MED_USE","YEAR","Phase","updrs_totscore_on")
  df$updrs_totscore_on[is.na(df$updrs_totscore_on)]<-median(df$updrs_totscore_on, na.rm = T)
  return(df)
  
}

Y0<-complement.f(0)
Y1<-complement.f(1)
Y2<-complement.f(2)
Y3<-complement.f(3)

phen.longitudinal<-rbind(Y0,Y1,Y2,Y3)
rm(Y0,Y1,Y2,Y3)


# plot the clinical trajectory of PD patients

p1<-ggplot(data=phen.longitudinal, aes(x=YEAR, y=updrs_totscore_on, group=PATNO)) + 
  geom_line(size=0.5,alpha=0.2,colour="#F4A460") + 
  ylab("MDS-UPDRS total score (ON)")+xlab("Year") + theme_pubr()

p2<-ggline(data=phen.longitudinal, 
           x="YEAR", y="updrs_totscore_on",add = "mean_se",color = "#FF8C00")+
  ylab("MDS-UPDRS total score (ON)")+xlab("Year")

pdf("./longitudinal_analysis/clinical.progression.pdf",width = 14,height = 5)
Rmisc::multiplot(plotlist = list(p1,p2),cols=2)
dev.off()
rm(p1,p2)

#construct linear mixed-effects models for each gene
RNA.longitudinal.lcpm.ANOVA.filter<-RNA.longitudinal.lcpm[rownames(RNA.longitudinal.lcpm) %in% gene_ezANOVA.filter$gene_id,] %>%
  t() %>% as.data.frame()

gene.lme.res<-sapply(RNA.longitudinal.lcpm.ANOVA.filter,function(x){
  Data<-data.frame(gene=x,phen.longitudinal)
  lme.model<-nlme::lme(updrs_totscore_on ~ gene + age + gen + YEAR + PD_MED_USE,
                       random = ~1|PATNO,correlation = corCAR1(),method = "REML",data = Data)
  list("p.value"=anova(lme.model)$`p-value`[2],"coefficient"=fixef(lme.model)[[2]])
})

gene.lme.res <- gene.lme.res %>% t() %>% as.data.frame() %>% rownames_to_column(var="gene_id") %>%
  lapply(unlist) %>% as.data.frame() %>% left_join(gene.anno) %>%
  dplyr::mutate(FDR=p.adjust(p.value,method="fdr"))


# GSEA
gsea.f<-function(gene_df){
  
  gene_df$entrez <- mapIds(org.Hs.eg.db,
                           keys=gene_df$gene_name,
                           column="ENTREZID",
                           keytype="SYMBOL",
                           multiVals="first")
  
  gene_df <- gene_df %>%  arrange(p.value) %>% dplyr::select(coefficient,entrez) %>%
    filter(entrez!="NA") %>%
    dplyr::distinct(entrez,.keep_all = T)
  
  # prepare genelist for GSEA
  geneList <- gene_df$coefficient
  names(geneList) = gene_df$entrez
  set.seed(123)
  geneList = sort(geneList, decreasing = TRUE)
  
  gsea.kegg <- gseKEGG(geneList     = geneList,
                       organism     = 'hsa',
                       nPerm        = 2000,
                       minGSSize    = 20,
                       pvalueCutoff = 0.05,
                       pAdjustMethod = "BH",
                        verbose      = FALSE)
  return(gsea.kegg)
}

gsea<-gsea.f(gene.lme.res)
kegg.gsea.res<-gsea@result

# plot the results of GSEA
paths <- c("hsa05022", "hsa05010", "hsa05014","hsa04140",
"hsa05012","hsa05016","hsa04728")
pdf("./longitudinal_analysis/gsea.kegg.pdf",width = 8.4,height = 6.2)
gseaplot2(gsea,paths,
          color = c("#FFA500","#377EB8","#2E8B57", "#9ACD32", "#D95F02", "#7570B3","#FF69B4"),
          subplots = 1:2)
dev.off()

# identification of progression genes
progression.genes<-filter(gene.lme.res,FDR<=0.05) %>% arrange(p.value)
nrow(progression.genes)
#203

# heatmap of progression genes
phen.longitudinal.sort<-arrange(phen.longitudinal,updrs_totscore_on)
phen.longitudinal.sort$PD_MED_USE<-factor(phen.longitudinal.sort$PD_MED_USE,
                              levels=c("0","1","2","3","4","5","6","7"),
                              labels=c("Unmedicated","Levodopa","DopamineAgonist",
                                       "Other","Levodopa_Other","Levodopa_DopamineAgonist",
                                       "DopamineAgonist_Other","Levodopa_DopamineAgonist_Other"))

phen.longitudinal.sort$YEAR<-factor(phen.longitudinal.sort$YEAR,
                        levels=c("0","1","2","3"),
                        labels=c("Baseline","Year1","Year2","Year3"))

phen.longitudinal.sort$gen<-factor(phen.longitudinal.sort$gen,
                       levels=c("1","2"),
                       labels=c("Male","Female"))


RNA.progression<-RNA.longitudinal.lcpm.ANOVA.filter[,colnames(RNA.longitudinal.lcpm.ANOVA.filter) %in%
                                                      progression.genes$gene_id] 
RNA.progression<-RNA.progression[match(phen.longitudinal.sort$PATNO_visit,rownames(RNA.progression)),] %>% t()

annotation_col = data.frame(
  Updrs_totscore_on = phen.longitudinal.sort$updrs_totscore_on,
  Treatment = phen.longitudinal.sort$PD_MED_USE,
  Time = phen.longitudinal.sort$YEAR,
  Age = phen.longitudinal$age,
  Gender=phen.longitudinal.sort$gen
)

rownames(annotation_col) = colnames(RNA.progression)

ann_colors = list(
  Updrs_totscore_on= c("white","firebrick"),
  Age = c("#F0FFFF","#FA8072"),
  Time = c(Baseline="#5F9EA0",Year1="#FFDAB9",Year2="#F4A460",Year3="#FA8072"),
  Gender = c(Male= "#5F9EA0", Female = "#FF8C00"),
  Treatment = c(Unmedicated="#D9D9D9",Levodopa="#FCCDE5",DopamineAgonist="#FFFFB3",
                Other="#BEBADA",Levodopa_Other="#FB8072",Levodopa_DopamineAgonist="#80B1D3",
                DopamineAgonist_Other="#FDB462",Levodopa_DopamineAgonist_Other="#B3DE69")
)

pdf("./longitudinal_analysis/heatmap.progression.genes.pdf",width = 13.5,height = 6.7)
pheatmap(RNA.progression, annotation_col = annotation_col,
         scale = "row",
         cluster_cols = F,show_colnames = F,show_rownames = F,
         fontsize = 10,
         treeheight_row = 25,
         cutree_rows = 2,
         border=F,
         annotation_colors = ann_colors,
         cluster_rows = T,
         color = colorRampPalette(c("#1874CD", "white", "#FA8072"))(15),
         breaks = seq(-1,1.4,length.out = 15)
)

dev.off()

#gene type of progression genes
gene_class<-table(progression.genes$gene_type) %>% as.data.frame()  %>% 
  mutate(Frequency=round(Freq/sum(Freq) * 100,1)) %>% dplyr::arrange(desc(Freq)) %>% 
  mutate(Gene_type=as.character(Var1)) %>% dplyr::select(-c(Var1,Freq))
others<- c(sum(gene_class$Frequency[5:14]),"others")
gene_class<-rbind(others,gene_class)[1:5,] %>% dplyr::mutate(Frequency=as.numeric(Frequency))

pdf("./longitudinal_analysis/progression.gene.class.pdf")
ggdonutchart(gene_class, "Frequency",
             lab.pos = "in",
             fill = "Gene_type",
             color = "white")+
  theme(legend.position="right")
dev.off()
rm(others)


# funtional analysis of progression genes 
progression.genes$entrez  <- mapIds(org.Hs.eg.db,
                                           keys=progression.genes$gene_name,
                                           column="ENTREZID",
                                           keytype="SYMBOL",
                                           multiVals="first")

GO <- enrichGO(
  gene          = progression.genes$entrez,
  keyType = "ENTREZID",
  OrgDb         = org.Hs.eg.db,
  ont           = "ALL",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE)@result

kegg<-enrichKEGG(
  gene          = progression.genes$entrez,
  keyType     = "kegg",
  organism   = 'hsa',
  pvalueCutoff      = 0.05,
  pAdjustMethod     = "BH",
  qvalueCutoff  = 0.05
)@result



pdf("./longitudinal_analysis/GO.progression.genes.pdf",width = 7,height = 3)
ggplot(GO[c(2,4,14,27),] %>% mutate(logP=-log10(pvalue)), 
       aes(x=fct_reorder(Description,logP), y=logP)) +
  geom_bar(stat="identity",fill = '#6094FF') +
  labs(x=" ",y="-log10(Ajusted P)") +
  coord_flip() +
  theme_bw()+
  theme(panel.grid = element_blank(),panel.border = element_rect(size=1))+
  theme(axis.ticks = element_line(size=1),axis.text = element_text(size=14,colour = "black",angle = 00),
        axis.title = element_text(size=16,color="black"))+
  theme(legend.text = element_text(size=14,colour = "black"))
dev.off()

# cell marker enrichment of progression genes
cell_markers <- vroom::vroom('http://bio-bigdata.hrbmu.edu.cn/CellMarker/download/Human_cell_markers.txt') %>%
  tidyr::unite("cellMarker", tissueType, cellName, sep=", ") %>%
  dplyr::select(cellMarker, geneID) %>%
  dplyr::mutate(geneID=stringr::str_extract_all(geneID, "[:digit:]+"))

progression.genes.cellMarker <- enricher(progression.genes$entrez, TERM2GENE=cell_markers, minGSSize=1)@result

pdf("./longitudinal_analysis/cell.markers.progression.genes.pdf",width = 7,height = 4)
ggplot(progression.genes.cellMarker[c(2,3,4,6,8,11,12),] %>% mutate(logP=-log10(pvalue)), 
       aes(x=fct_reorder(Description,logP), y=logP)) +
  geom_bar(stat="identity",fill = '#6094FF') +
  labs(x=" ",y="-log10(Ajusted P)") +
  coord_flip() +
  theme_bw()+
  theme(panel.grid = element_blank(),panel.border = element_rect(size=1))+
  theme(axis.ticks = element_line(size=1),axis.text = element_text(size=14,colour = "black",angle = 00),
        axis.title = element_text(size=16,color="black"))+
  theme(legend.text = element_text(size=14,colour = "black"))
dev.off()


# investigate the relationship between progression genes and PD-associated genes in PPI
PD.DisGeNET<-fread("./annotation/curated_gene_disease_associations.tsv",
                   sep = "auto",stringsAsFactors = F,header = T)  %>%
  filter(diseaseName=="Parkinson Disease")
#intersect(progression.genes$gene_name,PD.DisGeNET$geneSymbol)

PPI<-fread("./annotation/9606.protein.links.v11.0.txt",header = T)
PPI.infor<-fread("./annotation/9606.protein.info.v11.0.txt",header = T)[,c(1,2)]
PPI.700<- PPI %>% filter(combined_score>=700) %>%
  inner_join(PPI.infor,by=c("protein1"="protein_external_id")) %>%
  inner_join(PPI.infor,by=c("protein2"="protein_external_id"))

rm(PPI,PPI.infor)

progression.ppi<-PPI.700[PPI.700$preferred_name.x %in% progression.genes$gene_name,]
progression.PD.gene.ppi<-progression.ppi[progression.ppi$preferred_name.y %in% PD.DisGeNET$geneSymbol,]

nodes.type<-data.frame(gene_name=c(unique(progression.PD.gene.ppi$preferred_name.x),
                      unique(progression.PD.gene.ppi$preferred_name.y)),
           type=c(rep("progression",length(unique(progression.PD.gene.ppi$preferred_name.x))),
           rep("PD.gene",length(unique(progression.PD.gene.ppi$preferred_name.y)))))

write.csv(progression.PD.gene.ppi,file="./longitudinal_analysis/progression.PD.gene.ppi.csv", row.names = F,quote = F)
write.csv(nodes.type,file = "./longitudinal_analysis/node.type.csv", row.names = F,quote = F)

# bootstrap test
permutation_test<-function(gene_number,edge_count){
  set.seed(123)
  N<-c()
  for(i in 1:1000){
    geneList<-sample(gene.anno$gene_name,gene_number)
    
    geneList.ppi<-PPI.700[PPI.700$preferred_name.x %in% geneList,]
    geneList.PD.genes.ppi<-geneList.ppi[geneList.ppi$preferred_name.y %in% PD.DisGeNET$geneSymbol,]
    N=c(N,nrow(geneList.PD.genes.ppi))
  }
  p<-sum(N>edge_count)/1000
  return(p)
}

p_value<-permutation_test(nrow(progression.genes),nrow(progression.PD.gene.ppi))

# save data
saveRDS(progression.genes,file = "./longitudinal_analysis/progression.genes.rds")
save.image("./longitudinal_analysis/longitudinal_analysis.RData")
