library(ohchibi)
library(DESeq2)
library(palettesPM)
library(ggpubr)


setwd('/home/isai/Documents/results/github/variovoraxRGI/scripts/')
set.seed(seed = 130816)
#Load dataset
Dat <- readRDS(file = "../cleandata/dat_amplicon_castrillo2017.RDS")

##### ASV Level #####

Dat_rab <- Dat$RelativeAbundance
Dat_rab <- subset.Dataset(x = Dat_rab,Genotype == "Col_0" | 
                            Genotype == "Soil",drop = T,clean = T)

#Read taxonomy
df_tax <- read.table(
  file = "../cleandata/df_tax_castrillo2017.tsv",
  header = T,sep = "\t",quote = "",comment.char = "")

melted <- Dat_rab$Tab %>% melt
colnames(melted) <- c("ASV_Id","DADA2_Header","RA")
chosen <- df_tax %>% subset(Genus == "Variovorax") %$% Id %>% as.character
melted <- which(melted$ASV_Id %in% chosen) %>% melted[.,] %>% droplevels

melted_mean <- aggregate(RA ~ DADA2_Header,melted, "mean")
melted_mean <- merge(melted_mean,Dat_rab$Map, by = "DADA2_Header")
melted_mean$Fraction <- melted_mean$Fraction %>% gsub(pattern = "EC",replacement = "Root") %>%
  factor(levels = c("Soil","Root"))
paleta <- palettesPM::pm.colors.fractions()
paleta <- paleta[c(10,7)]


p <- chibi.boxplot(Map = melted_mean,x_val = "Fraction",y_val = "RA",size_boxplot = 1,col_val = "Fraction",
                   mpalette = paleta,median_colored_as_points = T,
                   style = "open",size_point = 8,size_axis_text.y = 30,alpha_point = 0.8,
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
oh.save.pdf(p = p,outname = "amplicon_castrillo2017_variovoraxRA.pdf",outdir = "../figures/",height = 15,width = 10)

