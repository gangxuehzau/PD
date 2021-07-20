setwd("~/gxue/PPMI/RNA-seq/PPMI-IR3/IR3_Analysis/")
library(tidyverse)
library(data.table)
library(DESeq2)
library(BiocParallel)
library(apeglm)
library(clusterProfiler)
library(org.Hs.eg.db)

library(ggpubr)
library(ggthemes)
library(pheatmap)
library(enrichplot)
library(VennDiagram)

# load data
HC<-readRDS("./Baseline_QA/Baseline_HC.rna.phen.rds")
PD<-readRDS("./Baseline_QA/Baseline_PD.rna.phen.rds")
Prodromal<-readRDS("./Baseline_QA/Baseline_Prodromal.rna.phen.rds")
#biotypes and gene names
gene.anno<-fread("./annotation/gencode.v29.gene.anno.txt")

# DESeq2 analysis
deseq2.f<-function(data,phen,Ref,Coef){
  
  #keep the genes with counts >= 10 in more than 10% of the samples
  keep<-rowSums(data>=10)>=0.1*ncol(data)
  data<- data[keep,]
  rm(keep)
  
  dds <- DESeqDataSetFromMatrix(countData = data,
                                colData = phen,
                                design = ~ age + gen + condition)
  dds$condition <- relevel(dds$condition, ref = Ref)
  
  register(MulticoreParam(12))
  dds<- DESeq(dds,parallel=TRUE) 
  res <- results(dds) 
  res.shr <- lfcShrink(dds, coef=Coef, parallel=TRUE,type="apeglm")
  
  res.shr<-as.data.frame(res.shr) %>% rownames_to_column("gene_id") %>% 
    left_join(gene.anno,by="gene_id")
  
  return(res.shr)
}


# Prodromal vs HC
HC.Prodromal.rna<-cbind(HC[[2]],Prodromal[[2]])
HC.Prodromal.phen<-rbind(HC[[1]],Prodromal[[1]]) %>% dplyr::select(PATNO,age,gen,condition) %>%
  mutate_at(c("condition","gen"),as.factor) %>% dplyr::mutate(age=scale(age)) %>% 
  tibble::column_to_rownames("PATNO") 
#identical(colnames(HC.Prodromal.rna),rownames(HC.Prodromal.phen))
Prodromal.HC.res<-deseq2.f(HC.Prodromal.rna,HC.Prodromal.phen,"HC","condition_Prodromal_vs_HC")


#PD vs Prodromal
Prodromal.PD.rna<-cbind(Prodromal[[2]],PD[[2]])
Prodromal.PD.phen <- rbind(Prodromal[[1]],PD[[1]]) %>% dplyr::select(PATNO,age,gen,condition) %>%
  mutate_at(c("condition","gen"),as.factor) %>% dplyr::mutate(age=scale(age)) %>% 
  tibble::column_to_rownames("PATNO") 
#identical(colnames(Prodromal.PD.rna),rownames(Prodromal.PD.phen))
PD.Prodromal.res<-deseq2.f(Prodromal.PD.rna,Prodromal.PD.phen,"Prodromal","condition_PD_vs_Prodromal")


# GSEA 
gsea.f<-function(gene_df){
  gene_df$entrez <- mapIds(org.Hs.eg.db,
                           keys=gene_df$gene_name,
                           column="ENTREZID",
                           keytype="SYMBOL",
                           multiVals="first")
  
  gene_df <- gene_df %>% filter(padj!="NA")%>%
    arrange(pvalue) %>% dplyr::select(log2FoldChange,entrez) %>% 
    filter(entrez!="NA") %>%
    dplyr::distinct(entrez,.keep_all = T)
  
  # prepare genelist for GSEA
  geneList <- gene_df$log2FoldChange
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

Prodromal.HC.gsea<-gsea.f(Prodromal.HC.res)
PD.Prodromal.gsea<-gsea.f(PD.Prodromal.res)
Prodromal.HC.gsea.res<-Prodromal.HC.gsea@result
PD.Prodromal.gsea.res<-PD.Prodromal.gsea@result

# plot the results of GSEA 
paths <- c("hsa04020", "hsa04360")
gseaplot2(Prodromal.HC.gsea,paths,color = c("#F16C6B","#6094FF"),subplots = 1:2)
ggsave("./Case-control/gsea.Prodromal.HC.pdf",width = 7,height = 4,units = "in")

gseaplot2(PD.Prodromal.gsea,paths,color = c("#F16C6B","#6094FF"),subplots = 1:2)
ggsave("./Case-control/gsea.PD.Prodromal.pdf",width = 7,height = 4,units = "in")


# cut-off of DEGs (Prodromal vs HC)
cut_off_pvalue = 0.05  
cut_off_logFC.Pro.HC = 0.25
Prodromal.HC.res$change = ifelse(Prodromal.HC.res$padj < cut_off_pvalue & abs(Prodromal.HC.res$log2FoldChange) >= cut_off_logFC.Pro.HC, 
                                 ifelse(Prodromal.HC.res$log2FoldChange >= cut_off_logFC.Pro.HC ,'Up','Down'),
                                 'Stable')

# Volcano Plot
pdf("./Case-control/Volcano_Prodromal.HC.pdf",width=7.5,height=5)  
ggscatter(Prodromal.HC.res%>% filter(padj!="NA") %>% dplyr::mutate(logP=-log10(padj)),
          x = "log2FoldChange", y = "logP",
          color = "change",
          palette = c("#2f5688","#BBBBBB","#CC0000"),
          size = 1,font.label = 8,repel = T) +
  theme_base()+ylab("-log10(Ajusted P-value)")+xlab("log2FoldChange (Prodromal vs HC)")+
  geom_hline(yintercept = 1.3,linetype="dashed")+
  geom_vline(xintercept = c(-cut_off_logFC.Pro.HC,cut_off_logFC.Pro.HC),linetype="dashed")+
  guides(color=guide_legend(title=NULL))+
  xlim(-1,1)+ylim(0,4)
dev.off()

#cut-off of DEGs (PD vs Prodromal)
cut_off_pvalue = 0.05 
cut_off_logFC.PD.pro = 0.35

PD.Prodromal.res$change = ifelse(PD.Prodromal.res$padj < cut_off_pvalue & abs(PD.Prodromal.res$log2FoldChange) >= cut_off_logFC.PD.pro, 
                                 ifelse(PD.Prodromal.res$log2FoldChange >= cut_off_logFC.PD.pro ,'Up','Down'),
                                 'Stable')

# Volcano Plot
pdf("./Case-control/Volcano_PD.prodromal.pdf",width=7.5,height=5)  
ggscatter(PD.Prodromal.res %>% filter(padj!="NA") %>% dplyr::mutate(logP=-log10(padj)),
          x = "log2FoldChange", y = "logP",
          color = "change",
          palette = c("#2f5688","#BBBBBB","#CC0000"),
          size = 1,font.label = 8,repel = T) +
  theme_base()+
  ylab("-log10(Ajusted P-value)")+xlab("log2FoldChange (PD vs Prodromal)")+
  geom_hline(yintercept = 1.3,linetype="dashed")+
  geom_vline(xintercept = c(-cut_off_logFC.PD.pro,cut_off_logFC.PD.pro),linetype="dashed")+
  guides(color=guide_legend(title=NULL))+
  xlim(-1.5,1)+ylim(0,7)
dev.off()


DEGs.Prodromal.HC<-filter(Prodromal.HC.res,padj<cut_off_pvalue & abs(log2FoldChange)>=cut_off_logFC.Pro.HC) 
table(DEGs.Prodromal.HC$change)
#Down   Up 
#354  530
DEGs.PD.Prodromal<-filter(PD.Prodromal.res,padj<cut_off_pvalue & abs(log2FoldChange)>=cut_off_logFC.PD.pro) 
table(DEGs.PD.Prodromal$change)
#Down   Up 
#6171   32 
overlap.DEGs<-inner_join(DEGs.Prodromal.HC,DEGs.PD.Prodromal,by="gene_id")
table(overlap.DEGs$change.x,overlap.DEGs$change.y)
#        Down  Up
# Down    0  16
# Up    398   0

venn.diagram(
  x = list(DEGs.Prodromal.HC$gene_id, DEGs.PD.Prodromal$gene_id),
  category.names = c("DEGs.Prodromal.HC" , "DEGs.PD.Prodromal"),
  output=TRUE,
  filename = "./Case-control/venn.tiff",
  fill = c("#85A6D8", "#FFCC99"),
  fontface = "bold"
)

write_csv(DEGs.Prodromal.HC %>% arrange(pvalue),file = "./Case-control/DEGs.Prodromal.HC.csv")
write_csv(DEGs.PD.Prodromal %>% arrange(pvalue),file = "./Case-control/DEGs.PD.Prodromal.csv")

# heatmap of overlapping DEGs (tipping-point genes)
HC.Prodromal.PD.counts<-readRDS("./Baseline_QA/Baseline_HC.Prodromal.PD.rna.rds")
overlap.DEGs.counts<-HC.Prodromal.PD.counts[rownames(HC.Prodromal.PD.counts) %in% overlap.DEGs$gene_id,] 

annotation_col<- bind_rows(HC[[1]],Prodromal[[1]],PD[[1]]) %>% 
  dplyr::select(PATNO,condition,age,gen) %>% mutate_at(c("condition","gen"),as.factor) %>%
  tibble::column_to_rownames("PATNO") %>% 
  plyr::rename(c("age"="Age")) %>% plyr::rename(c("gen"="Gender")) %>%
  plyr::rename(c("condition"="Condition"))
  
annotation_col$Gender<-factor(annotation_col$Gender,
                                   levels=c("1","2"),
                                   labels=c("Male","Female"))

#identical(rownames(annotation_col),colnames(HC.Prodromal.PD.counts))

ann_colors = list(
  Condition = c(HC = "#80B1D3", Prodromal = "#FDB462",PD = "#FB8072"),
  Age = c("#F0FFFF","#FA8072"),
  Gender = c(Male= "#377EB8", Female = "#4DAF4A")
)

pdf("./Case-control/heatmap.tipping-point.genes.pdf")
pheatmap(overlap.DEGs.counts,
         cluster_rows=T, cluster_cols=F,
         show_rownames=F, show_colnames = F,scale = "row",
         gaps_col = c(ncol(HC[[2]]),ncol(HC.Prodromal.rna)),
         cutree_rows = 2,
         annotation_col=annotation_col,
         annotation_colors=ann_colors,
         color = colorRampPalette(c("#1874CD", "white", "firebrick3"))(20),
         breaks = seq(-0.8,1.6,length.out = 20))
dev.off()

# gene type of overlapping DEGs (tipping-point genes) 
gene_class<-table(overlap.DEGs$gene_type.x) %>% as.data.frame()  %>% 
  mutate(Frequency=round(Freq/sum(Freq) * 100,1)) %>% dplyr::arrange(desc(Freq)) %>% 
  mutate(Gene_type=as.character(Var1)) %>% dplyr::select(-c(Var1,Freq))
others<- c(sum(gene_class$Frequency[5:14]),"others")
gene_class<-rbind(others,gene_class)[1:5,] %>% dplyr::mutate(Frequency=as.numeric(Frequency))

pdf("./Case-control/tipping.point.gene.class.pdf")
ggdonutchart(gene_class, "Frequency",
             lab.pos = "in",
             fill = "Gene_type",
             color = "white")+
  theme(legend.position="right")
dev.off()
rm(others)


# funtional analysis of overlapping DEGs (tipping-point genes)
overlap.DEGs$entrez  <- mapIds(org.Hs.eg.db,
                    keys=overlap.DEGs$gene_name.x,
                    column="ENTREZID",
                    keytype="SYMBOL",
                    multiVals="first")

GO <- enrichGO(
  gene          = overlap.DEGs$entrez,
  keyType = "ENTREZID",
  OrgDb         = org.Hs.eg.db,
  ont           = "ALL",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE)@result

kegg <- enrichKEGG(
  gene          = overlap.DEGs$entrez,
  keyType     = "kegg",
  organism   = 'hsa',
  pvalueCutoff      = 0.05,
  pAdjustMethod     = "BH",
  qvalueCutoff  = 0.05
)@result

pdf("./Case-control/GO.tipping.point.genes.pdf",width = 7,height = 3)
ggplot(GO[c(1,3,7,23),] %>% mutate(logP=-log10(pvalue)), 
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


# cell marker enrichment 
library(vroom)
cell_markers <- vroom::vroom('http://bio-bigdata.hrbmu.edu.cn/CellMarker/download/Human_cell_markers.txt') %>%
  tidyr::unite("cellMarker", tissueType, cellName, sep=", ") %>%
  dplyr::select(cellMarker, geneID) %>%
  dplyr::mutate(geneID=stringr::str_extract_all(geneID, "[:digit:]+"))

overlap.DEGs.cellMarker <- enricher(overlap.DEGs$entrez, TERM2GENE=cell_markers, minGSSize=1)@result

#investigate the relationship between progression genes and PD-associated genes in PPI
PD.DisGeNET<-fread("./annotation/curated_gene_disease_associations.tsv",
                   sep = "auto",stringsAsFactors = F,header = T)  %>%
  filter(diseaseName=="Parkinson Disease")
nrow(PD.DisGeNET)
#85
# overlap.DEGs[overlap.DEGs$gene_name.x%in%PD.DisGeNET$geneSymbol,]
# load PPI >= 700 
PPI<-fread("./annotation/9606.protein.links.v11.0.txt",header = T)
PPI.infor<-fread("./annotation/9606.protein.info.v11.0.txt",header = T)[,c(1,2)]
PPI.700<- PPI %>% filter(combined_score>=700) %>%
  inner_join(PPI.infor,by=c("protein1"="protein_external_id")) %>%
  inner_join(PPI.infor,by=c("protein2"="protein_external_id"))
rm(PPI,PPI.infor)

tipping.point.ppi<-PPI.700[PPI.700$preferred_name.x %in% overlap.DEGs$gene_name.x,]
tipping.point.PD.gene.ppi<-tipping.point.ppi[tipping.point.ppi$preferred_name.y %in% PD.DisGeNET$geneSymbol,]


nodes.type<-data.frame(gene_name=c(unique(tipping.point.PD.gene.ppi$preferred_name.x),
                                   unique(tipping.point.PD.gene.ppi$preferred_name.y)),
                       type=c(rep("tipping.point",length(unique(tipping.point.PD.gene.ppi$preferred_name.x))),
                              rep("PD.gene",length(unique(tipping.point.PD.gene.ppi$preferred_name.y)))))

write.csv(tipping.point.PD.gene.ppi,file="./Case-control/tipping.point.PD.gene.csv", row.names = F,quote = F)
write.csv(nodes.type,file = "./Case-control/node.type.csv", row.names = F,quote = F)

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

p_value<-permutation_test(nrow(overlap.DEGs),nrow(tipping.point.PD.gene.ppi))


# extract genes in the KEGG pathyway
pathyway.f<-function(id,gsea.res,gene_df){
  
  pathyway<-filter(gsea.res,ID==id) 
  pathyway.gene <-str_split(pathyway$core_enrichment,'/',simplify=T) %>% as.character() %>%
    mapIds(org.Hs.eg.db,
           keys=.,
           keytype="ENTREZID",
           column="SYMBOL",
           multiVals="first") %>% as.character()
  pathyway.gene<- gene_df[gene_df$gene_name %in% pathyway.gene,] 
  return(pathyway.gene)
}

#focus on pathyway: signaling pathways regulating pluripotency of stem cells (hsa04550)
gseaplot2(PD.Prodromal.gsea, geneSetID = "hsa04550",color = "red",subplots = 1:2)
ggsave("./Case-control/gsea.hsa04550.pdf",width = 5.5,height = 2.5,units = "in")
hsa04550.gene<-pathyway.f("hsa04550",PD.Prodromal.gsea.res,PD.Prodromal.res)
write_csv(hsa04550.gene,file = "./Case-control/hsa4550.gene.csv")

# save data
write_csv(overlap.DEGs %>% arrange(pvalue.x),file = "./Case-control/tipping.point.genes.csv")
saveRDS(overlap.DEGs,"./Case-control/tipping.point.genes.rds")
save.image("./Case-control/Case-control.RData")
