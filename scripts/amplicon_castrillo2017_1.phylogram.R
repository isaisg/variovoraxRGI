library(ohchibi)
library(palettesPM)
library(extrafont)
loadfonts(device = "pdf")


setwd('/home/isai/Documents/results/github/variovoraxRGI/scripts/')

#Load dataset
Dat <- readRDS(file = "../cleandata/dat_amplicon_castrillo2017.RDS")
Dat_raw <- Dat$RawCounts
Map <- Dat_raw$Map
Map$Count <- rep(1,nrow(Map))

dcast(data = Map,formula = Rep~Genotype+Fraction,
      fun.aggregate = sum,value.var = "Count")

Dat_raw <- Dat_raw %>% 
  subset.Dataset(Genotype != "pho1_2",drop = T,clean = T)

Dat_raw$Map$group <- paste0(Dat_raw$Map$Fraction,"_",Dat_raw$Map$Genotype) %>%
  factor

original_reads <- Dat_raw$Tab %>% sum
#Identify measurable ESVs
Dat_raw <- measurable_taxa(Dat = Dat_raw,min_samples_otu =10,min_reads_otu = 12,method = "absolute",clean = T)
filtered_reads <- Dat_raw$Tab %>% sum

(filtered_reads/original_reads)*100



#Create the deseq object
dds <- DESeqDataSetFromMatrix(countData = Dat_raw$Tab,colData = Dat_raw$Map,
                              design = ~ Rep + group)
dds <- DESeq(dds)

res <- results(object = dds,contrast = c("group","EC_Col_0","Soil_Soil")) %>%
  as.data.frame
res$Id <- rownames(res)
rownames(res) <- NULL

df_taxonomy <- read.table("../cleandata/df_tax_castrillo2017.tsv",header = T,
                          sep = "\t",comment.char = "",quote = "")
res <- merge(res,df_taxonomy, by = "Id")

res_up <- res %>% subset(padj < 0.1 & log2FoldChange > 0) %>% droplevels
total <- nrow(res_up)
res_up$Count <- rep(1,nrow(res_up))
res_up$Perc <- res_up$Count/(total)
res_up$Type <- rep("Bar",nrow(res_up))
res_up$Phylum <- res_up$Phylum %>% as.character
tochange <- which(!(res_up$Phylum %in% palettesPM::pm.names.phyla())) %>% res_up$Phylum[.]
res_up$Phylum[which(res_up$Phylum %in% tochange)] <- "Other"
res_up$Phylum <- res_up$Phylum %>% factor(levels = pm.names.phyla())
p <- ggplot(data = res_up,aes(Type,Perc, fill = Phylum)) +
  geom_bar(stat = "identity") + scale_fill_phyla() +
  theme_ohchibi() + 
  theme(
    legend.position = "none",
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    aspect.ratio = 10,
    panel.background = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank()
  ) + coord_fixed(expand = F) 

oh.save.pdf(p = p, outdir = "../figures/",
            outname = "amplicon_castrillo2017_enrichedroot_deseq2.pdf")
