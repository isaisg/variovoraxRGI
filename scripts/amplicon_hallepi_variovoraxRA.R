library(ohchibi)
library(DESeq2)
library(palettesPM)
library(ggpubr)

setwd('/home/isai/Documents/results/github/variovoraxRGI/scripts/')
set.seed(seed = 130816)
#Load dataset
Dat <- readRDS(file = "../cleandata/dat_amplicon_hallepi.RDS")

##### ASV Level #####

Dat_raw <- Dat$RawCounts

#Remove BulkSoil Samples
Dat_raw <- subset.Dataset(x = Dat_raw,subset = Fraction!="BulkSoil",drop = T,clean = T) %>%
  subset.Dataset(x = .,Genotype == "Col-0",
                 drop = T,clean = T)  


#Read taxonomy
df_tax <- read.table(
  file = "../cleandata/df_tax_hallepi.tsv",
  header = T,sep = "\t",quote = "",comment.char = "")
df_tax <- df_tax[,-2]



original_reads <- Dat_raw$Tab %>% sum
#Identify measurable ESVs
Dat_raw <- measurable_taxa(Dat = Dat_raw,min_samples_otu =5,min_reads_otu = 12,method = "absolute",clean = T)
filtered_reads <- Dat_raw$Tab %>% sum

(filtered_reads/original_reads)*100


#Create the deseq object
dds <- DESeqDataSetFromMatrix(countData = Dat_raw$Tab,colData = Dat_raw$Map,
                              design = ~ Fraction + Phosphate)

gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans <- apply(counts(dds), 1, gm_mean)
dds <- estimateSizeFactors(dds, geoMeans = geoMeans)
dds <- DESeq(dds, fitType="local")

res <- results(object = dds,contrast = c("Fraction","Root","Soil")) %>% as.data.frame
res$ASV_Id <- rownames(res)
res <- merge(res,df_tax, by = "ASV_Id")
res %>% subset(Genus == "Variovorax")


res <- results(object = dds,contrast = c("Fraction","Shoot","Soil")) %>% as.data.frame
res$ASV_Id <- rownames(res)
res <- merge(res,df_tax, by = "ASV_Id")
res %>% subset(Genus == "Variovorax")


### Relative abundance
Dat_rab <- Dat$RelativeAbundance

#Remove BulkSoil Samples
Dat_rab <- subset.Dataset(x = Dat_rab,subset = Fraction!="BulkSoil",drop = T,clean = T) %>%
  subset.Dataset(x = .,Genotype == "Col-0",
                 drop = T,clean = T)  
melted <- Dat_rab$Tab %>% melt
colnames(melted) <- c("ASV_Id","DADA2_Header","RA")
#chosen <- df_tax %>% subset(Genus == "Variovorax") %$% ASV_Id %>% as.character
chosen <- res %>% subset(Genus == "Variovorax") %$% ASV_Id %>% as.character %>% unique
melted <- which(melted$ASV_Id %in% chosen) %>% 
  melted[.,] %>% droplevels

melted_mean <- aggregate(RA ~ DADA2_Header,melted, "mean")


melted_mean <- merge(melted_mean,Dat_rab$Map, by = "DADA2_Header")
melted_mean <- melted_mean %>% subset(Fraction != "Shoot") %>% droplevels


paleta <- palettesPM::pm.colors.fractions()
paleta <- paleta[c(10,7)]

p <- chibi.boxplot(Map = melted_mean,x_val = "Fraction",y_val = "RA",size_boxplot = 1,col_val = "Fraction",
                   mpalette = paleta,median_colored_as_points = T,
                   style = "mix",size_point = 8,size_axis_text.y = 30,alpha_point = 0.8,
                   size_axis_title.x = 40,size_axis_title.y = 40,size_median = 6,median_color = "black",
                   legend_proportion_size = 4,size_legend_text = 30,strip_text_size = 30) + 
  ylab(label = "Variovorax Relative Abundance")  + xlab(label = "Fraction") + 
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    legend.position = "none"
  ) +
  scale_y_continuous(limits = c(0,0.004))
p <- p+ stat_compare_means()


oh.save.pdf(p = p,outname = "amplicon_hallepi_variovoraxRA.pdf",outdir = "../figures/",height = 15,width = 10)
