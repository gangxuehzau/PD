setwd("/home/kyzhou/gxue/PPMI/RNA-seq/PPMI-IR3/IR3_Analysis/")
library(data.table)
library(tidyverse)

# load data
tipping.point.genes<-readRDS("./Case-control/tipping.point.genes.rds")
progression.genes<-readRDS("./longitudinal_analysis/progression.genes.rds")
gene.anno<-fread("./annotation/gencode.v29.gene.anno.txt")

#overlap
library(VennDiagram)
venn.diagram(
  x = list(tipping.point.genes$gene_id,progression.genes$gene_id),
  category.names = c( "tipping.point.genes","progression.genes"),
  output=TRUE,
  filename = "./onset_progression/venn.tiff",
  fontface = "bold"
)


#investigate the relationship between tipping-point genes and progression genes in PPI
# load PPI >= 700 
PPI<-fread("./annotation/9606.protein.links.v11.0.txt",header = T)
PPI.infor<-fread("./annotation/9606.protein.info.v11.0.txt",header = T)[,c(1,2)]
PPI.700<- PPI %>% filter(combined_score>=700) %>%
  inner_join(PPI.infor,by=c("protein1"="protein_external_id")) %>%
  inner_join(PPI.infor,by=c("protein2"="protein_external_id"))

# ppi interactions
tipping.point.ppi<-PPI.700[PPI.700$preferred_name.x %in% tipping.point.genes$gene_name.x,]
tipping.point.progression.ppi<-tipping.point.ppi[tipping.point.ppi$preferred_name.y %in% progression.genes$gene_name,]

# nodes type
nodes.type<-data.frame(gene_name=c(unique(tipping.point.progression.ppi$preferred_name.x),
                                   unique(tipping.point.progression.ppi$preferred_name.y)),
                       type=c(rep("tipping.point",length(unique(tipping.point.progression.ppi$preferred_name.x))),
                              rep("Progression",length(unique(tipping.point.progression.ppi$preferred_name.y)))))

# bootstrap test
permutation_test<-function(gene_number,edge_count){
  set.seed(123)
  N<-c()
  for(i in 1:1000){
    geneList<-sample(gene.anno$gene_name,gene_number)
    
    geneList.ppi<-PPI.700[PPI.700$preferred_name.x %in% geneList,]
    geneList.PD.genes.ppi<-geneList.ppi[geneList.ppi$preferred_name.y %in% tipping.point.genes$gene_name.x,]
    N=c(N,nrow(geneList.PD.genes.ppi))
  }
  p<-sum(N>edge_count)/1000
  return(p)
}

p_value<-permutation_test(nrow(progression.genes),nrow(tipping.point.progression.ppi))

write_csv(tipping.point.progression.ppi,file="./onset_progression/ppi.tipping.point.progression.csv")
write_csv(nodes.type,file = "./onset_progression/ppi.tipping.point.progression.nodes.type.csv")


#TF regulatory_networks
human.regulatory<-fread("./annotation/regulatory_networks/human.source")[,c(1,3)]
colnames(human.regulatory)<-c("TF","target")

tipping.point.regulatory<-human.regulatory[human.regulatory$TF%in%tipping.point.genes$gene_name.x,]
tipping.point.progression.regulatory<-tipping.point.regulatory[tipping.point.regulatory$target%in%progression.genes$gene_name, ]


# nodes type
nodes.type<-data.frame(gene_name=c(unique(tipping.point.progression.regulatory$TF),
                                   unique(tipping.point.progression.regulatory$target)),
                       type=c(rep("tipping.point",length(unique(tipping.point.progression.regulatory$TF))),
                              rep("Progression",length(unique(tipping.point.progression.regulatory$target)))))


write.csv(tipping.point.progression.regulatory,
          file = "./onset_progression/TF.regulatory.csv",
          row.names = F,quote = F)

write.csv(nodes.type,
          file = "./onset_progression/regulatory.tipping.point.progression.nodes.type.csv",
          row.names = F,quote = F)

save.image("./onset_progression/onset.progression.RData")





