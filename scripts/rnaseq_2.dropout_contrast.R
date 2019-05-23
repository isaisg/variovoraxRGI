library(ohchibi)
library(pheatmap)
library(clusterProfiler)
library(org.At.tair.db)
library(DESeq2)
library(paletteer)
library(palettesPM)

setwd('/home/isai/Documents/results/github/variovoraxRGI/scripts/')
set.seed(seed = 130816)


Dat_rnaseq_comb_do <- readRDS(file = "../cleandata/dat_rnaseq_plants_tripartitecl28cl14_variovoraxdropout.RDS")

Dat_raw <- Dat_rnaseq_comb_do$Dat_rnaseq

Dat_sub <- Dat_raw %>% subset.Dataset(subset = Experiment == "VarioDO",drop = T,clean = T)


#Create design for the transformed object
dds<-DESeqDataSetFromMatrix(countData = Dat_sub$Tab,colData = Dat_sub$Map,
                            design =~ Syncom)


dds <- DESeq(dds)
vsd <- vst(dds, blind=FALSE)


#Direct contrasts
alpha_cutoff <- 0.1

res_dov <- results(object = dds,contrast = c("Syncom","DOV","NB")) %>%
  as.data.frame
inducted_genes_dov <- res_dov %>% subset(padj < alpha_cutoff & log2FoldChange > 0) %>% rownames %>%
  as.character

res_full <- results(object = dds,contrast = c("Syncom","DOV","Full")) %>%
  as.data.frame

res_projected <- match(inducted_genes_dov,rownames(res_full)) %>% res_full[.,] %>% droplevels

#Classification scheme
res_projected_diff <- res_projected %>% subset(padj < alpha_cutoff) %>% droplevels
res_projected_ndiff <- res_projected %>% subset(padj >= alpha_cutoff) %>% droplevels

#Create list
inducted_repressed <- res_projected_diff %>% rownames %>% as.character
inducted_inducted <- res_projected_ndiff %>% rownames %>% as.character
lista <- list(
  IR = inducted_repressed,
  II = inducted_inducted
)


####### Figures ####
Tab_vsd <- match(inducted_repressed,rownames(assay(vsd))) %>%
  assay(vsd)[.,]

Tab_vsd_all <- assay(vsd) %>% t %>% scale %>% t %>% as.matrix %>%
  melt
colnames(Tab_vsd_all) <- c("genes","Sample_Id","value")
Tab_vsd_all <- merge(Tab_vsd_all,Dat_sub$Map,by = "Sample_Id")
mdf_all <- dcast(Tab_vsd_all,genes ~Syncom,value.var = "value",fun.aggregate = mean) %>%
  melt
### Scale the data
Tab_vsd_Z <- Tab_vsd %>% t %>% scale %>% t %>% as.matrix
melted_Tab_vsd_Z <- Tab_vsd_Z %>% melt
colnames(melted_Tab_vsd_Z) <- c("genes","Sample_Id","value")
melted_Tab_vsd_Z <- merge(melted_Tab_vsd_Z,Dat_sub$Map,by = "Sample_Id")

mdf <- dcast(melted_Tab_vsd_Z,genes ~Syncom,value.var = "value",fun.aggregate = mean) %>%
  melt



mdf$variable <- mdf$variable %>% factor(levels = c("NB","DOB","DOV","DOVB","Full"))



#Determine the quantile for the Full community 
p <- ggplot(data = mdf, aes(x = variable, y = value)) +
  geom_boxplot(color="black",outlier.colour = NA,
               position = position_dodge(width = 0.9), 
               size=2,alpha = 0)
dat <- ggplot_build(p)$data[[1]]
p <- p + geom_jitter(position = position_jitter(0.2),
            size = 12, shape = 21, col = "black",stroke=2,alpha=1)
p <- p + geom_segment(data=dat, aes(x=xmin, xend=xmax,
                                    y=middle, yend=middle), colour="red", size=8)
p <- p + theme_ohchibi() + 
  theme(panel.background  =element_rect(fill = "#E8FFFA"),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.title.x = element_blank()
        ) + ylab(label = "z-score") 


oh.save.pdf(p = p,outname = "rnaseq_contrast_dropout_inducteddov_repressedfull.pdf",
            height = 20,width = 18,outdir = "../figures/")

#Save the plot
saveRDS(object = p,file = "../figures/rnaseq_contrast_dropout_inducteddov_repressedfull.RDS")

#Save structures of expression and glm results
Dat <- list(dds = dds, Expression = mdf_all,
            chosen_genes = inducted_repressed)

saveRDS(object = Dat,file = "../cleandata/dataset_rnaseq_contrast_dropout_inducteddov_repressedfull.RDS")


