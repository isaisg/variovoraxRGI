library(ohchibi)
library(DESeq2)

set.seed(130816)
setwd('/home/isai/Documents/results/github/variovoraxRGI/scripts/')
dir.create("../cleandata")
dir.create("../figures")

#Create expression supplementary datasets S10 and S11
Dat_vdo <- readRDS(file = "../cleandata/dataset_rnaseq_contrast_dropout_inducteddov_repressedfull.RDS")
Dat_tri <- readRDS(file = "../cleandata/dataset_rnaseq_contrast_dropout_inductedcl28_repressedcl28cl14.RDS")

chosen_genes <- union(Dat_vdo$chosen_genes,Dat_tri$chosen_genes)


#Tripartite
exp <- Dat_tri$Expression
exp <- dcast(data = exp,formula = genes~variable,
             fun.aggregate = mean,value.var = "value")
dds <- Dat_tri$dds
res_rgi <- results(dds,contrast = c("Syncom","CL28","NB")) %>% as.data.frame
res_rgirev <- results(dds,contrast = c("Syncom","CL28","CL28_CL14")) %>% as.data.frame
res_rgi <- res_rgi[,c(2,5,6)]
res_rgirev <- res_rgirev[,c(2,5,6)]
res_rgi$genes <- rownames(res_rgi)
res_rgirev$genes <- rownames(res_rgirev)

exp <- match(chosen_genes,exp$genes) %>%
  exp[.,]


df <- merge(exp,res_rgi, by = "genes") %>%
  merge(res_rgirev, by = "genes")
colnames(df) <- c("Gene","Mean_zscore_VariovoraxCL14","Mean_zscore_ArthrobacterCL28","Mean_zscore_Duoassociation","Mean_zscore_NoBacteriaTripartite",
                  "GLM_ArthrobacterCL28vsNB_log2FoldChange","GLM_ArthrobacterCL28vsNB_pvalue","GLM_ArthrobacterCL28vsNB_padj",
                  "GLM_ArthrobacterCL28vsDuoassociation_log2FoldChange","GLM_ArthrobacterCL28vsDuoassociation_pvalue",
                  "GLM_ArthrobacterCL28vsDuoassociation_padj"
                  )
df_tri <- df[,c(1,5,2,3,4,6:11)]

#Dropout
exp <- Dat_vdo$Expression
exp <- dcast(data = exp,formula = genes~variable,
             fun.aggregate = mean,value.var = "value")
dds <- Dat_vdo$dds
res_rgi <- results(dds,contrast = c("Syncom","DOV","NB")) %>% as.data.frame
res_rgirev <- results(dds,contrast = c("Syncom","DOV","Full")) %>% as.data.frame
res_rgi <- res_rgi[,c(2,5,6)]
res_rgirev <- res_rgirev[,c(2,5,6)]
res_rgi$genes <- rownames(res_rgi)
res_rgirev$genes <- rownames(res_rgirev)

exp <- match(chosen_genes,exp$genes) %>%
  exp[.,]


df <- merge(exp,res_rgi, by = "genes") %>%
  merge(res_rgirev, by = "genes")
df <- df[,-c(2,4)]
colnames(df) <- c("Gene","Mean_zscore_VariovoraxDropout","Mean_zscore_Full","Mean_zscore_NoBacteriaDropout",
                  "GLM_VariovoraxDropoutvsNB_log2FoldChange","GLM_VariovoraxDropoutvsNB_pvalue","GLM_VariovoraxDropoutvsNB_padj",
                  "GLM_VariovoraxDropoutvsFull_log2FoldChange","GLM_VariovoraxDropoutvsFull_pvalue",
                  "GLM_VariovoraxDropoutvsFull_padj"
)
df_dov <- df[,c(1,4,3,2,5:10)]

df <- merge(df_tri,df_dov, by = "Gene")

write.table(x = df,file = "../rawdata/dataS11_56genes_RGI.tsv",
            append = F,quote = F,sep = "\t",row.names = F,col.names = T)


##################### 12 Genes ############################
df_auxin_regulon <- read.table(file = "../rawdata/rnaseq_regulon_iaa_natchembio.tsv",header = T,sep = "\t")
df_auxin_regulon <- with(df_auxin_regulon,order(-logFC)) %>%
  df_auxin_regulon[.,]
df_auxin_regulon %>% subset(logFC >=2) %>% dim
sig_genes <- df_auxin_regulon %>% head(n = 13) %$% AGI.code %>%
  as.character
sig_genes <- sig_genes[-which(!(sig_genes %in% Dat_tri$Expression$genes))]

chosen_genes <- sig_genes

#Tripartite
exp <- Dat_tri$Expression
exp <- dcast(data = exp,formula = genes~variable,
             fun.aggregate = mean,value.var = "value")
dds <- Dat_tri$dds
res_rgi <- results(dds,contrast = c("Syncom","CL28","NB")) %>% as.data.frame
res_rgirev <- results(dds,contrast = c("Syncom","CL28","CL28_CL14")) %>% as.data.frame
res_rgi <- res_rgi[,c(2,5,6)]
res_rgirev <- res_rgirev[,c(2,5,6)]
res_rgi$genes <- rownames(res_rgi)
res_rgirev$genes <- rownames(res_rgirev)

exp <- match(chosen_genes,exp$genes) %>%
  exp[.,]


df <- merge(exp,res_rgi, by = "genes") %>%
  merge(res_rgirev, by = "genes")
colnames(df) <- c("Gene","Mean_zscore_VariovoraxCL14","Mean_zscore_ArthrobacterCL28","Mean_zscore_Duoassociation","Mean_zscore_NoBacteriaTripartite",
                  "GLM_ArthrobacterCL28vsNB_log2FoldChange","GLM_ArthrobacterCL28vsNB_pvalue","GLM_ArthrobacterCL28vsNB_padj",
                  "GLM_ArthrobacterCL28vsDuoassociation_log2FoldChange","GLM_ArthrobacterCL28vsDuoassociation_pvalue",
                  "GLM_ArthrobacterCL28vsDuoassociation_padj"
)
df_tri <- df[,c(1,5,2,3,4,6:11)]

#Dropout
exp <- Dat_vdo$Expression
exp <- dcast(data = exp,formula = genes~variable,
             fun.aggregate = mean,value.var = "value")
dds <- Dat_vdo$dds
res_rgi <- results(dds,contrast = c("Syncom","DOV","NB")) %>% as.data.frame
res_rgirev <- results(dds,contrast = c("Syncom","DOV","Full")) %>% as.data.frame
res_rgi <- res_rgi[,c(2,5,6)]
res_rgirev <- res_rgirev[,c(2,5,6)]
res_rgi$genes <- rownames(res_rgi)
res_rgirev$genes <- rownames(res_rgirev)

exp <- match(chosen_genes,exp$genes) %>%
  exp[.,]


df <- merge(exp,res_rgi, by = "genes") %>%
  merge(res_rgirev, by = "genes")
df <- df[,-c(2,4)]
colnames(df) <- c("Gene","Mean_zscore_VariovoraxDropout","Mean_zscore_Full","Mean_zscore_NoBacteriaDropout",
                  "GLM_VariovoraxDropoutvsNB_log2FoldChange","GLM_VariovoraxDropoutvsNB_pvalue","GLM_VariovoraxDropoutvsNB_padj",
                  "GLM_VariovoraxDropoutvsFull_log2FoldChange","GLM_VariovoraxDropoutvsFull_pvalue",
                  "GLM_VariovoraxDropoutvsFull_padj"
)
df_dov <- df[,c(1,4,3,2,5:10)]

df <- merge(df_tri,df_dov, by = "Gene")

write.table(x = df,file = "../rawdata/dataS12_12genes_IAAmarkers.tsv",
            append = F,quote = F,sep = "\t",row.names = F,col.names = T)
