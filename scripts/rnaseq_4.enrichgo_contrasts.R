library(ohchibi)
library(clusterProfiler)
library(org.At.tair.db)


set.seed(seed = 130816)
setwd('/home/isai/Documents/results/github/variovoraxRGI/scripts/')
dir.create("../figures")


trinity <- readRDS(file = "../cleandata/dataset_rnaseq_contrast_dropout_inductedcl28_repressedcl28cl14.RDS") %$% 
  chosen_genes

dropout <- readRDS(file = "../cleandata/dataset_rnaseq_contrast_dropout_inducteddov_repressedfull.RDS") %$%
  chosen_genes

intersection <- intersect(trinity,dropout)
tonly <- which(!(trinity %in% intersection)) %>% trinity[.]
donly <-  which(!(dropout %in% intersection)) %>% dropout[.]


ego <- enrichGO(gene          = intersection,
                keyType = "TAIR",
                OrgDb         = org.At.tair.db,
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                qvalueCutoff  = 0.1,
                readable      = TRUE)
p <- emapplot(ego,showCategory = 42,color = "p.adjust") 
oh.save.pdf(p = p,outname = "rnaseq_enrichGO_18genesintersection.pdf",outdir = "../figures",
            width = 10,height = 10)


#Hypergeometric test
arabidopsis_genes_biggest_matrix = 22604
q = intersection %>% length
m = trinity %>% length
k = dropout %>% length
n = arabidopsis_genes_biggest_matrix - m
phyper(q = q-1,m = m,n = n,k = k,lower.tail = F)
