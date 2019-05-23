library(ohchibi)
library(palettesPM)
library(extrafont)
loadfonts(device = "pdf") 


setwd('/home/isai/Documents/results/github/variovoraxRGI/scripts/')
set.seed(130816)
## Abundance
Dat <- readRDS(file = "../cleandata/dat_amplicon_4stresses_useq97.RDS")
Dat_rab <- Dat$RelativeAbundance
taxons_remove <- Dat_rab$Tab %>% rownames %>% grep(pattern = "OTU",value = T)
Dat_rab <- remove_taxons(Dat = Dat_rab,taxons = taxons_remove)

#Collapse taxonomy
Dat_phyla <- Dat_rab %>% collapse_by_taxonomy.Dataset(Dat = .,level = 3)

#Split Tab and Map to create another structure
Tab  <- Dat_phyla$Tab
Map <- Dat_phyla$Map

rownames(Tab) <- Tab %>% rownames %>%
  strsplit(split = "\\;") %>% unlist %>%
  grep(pattern = "p__",value = T) %>% 
  gsub(pattern = "p__",replacement = "") %>%
  gsub(pattern = " ",replacement = "")

#rownames(Tab)[4] <- "Alphaproteobacteria"
#rownames(Tab)[5] <- "Gammaproteobacteria"
#rownames(Tab)[6] <- "Gammaproteobacteria"

#Gammaproteobacteria <- Tab[5:6,] %>% colSums()
#Tab_other <- Tab[1:4,]
#Tab <- rbind(Tab_other,Gammaproteobacteria)

mphyla <- palettesPM::pm.names.phyla()
paleta <- palettesPM::pm.colors.phyla()
extra <- c("#CC873A","#FF6100")
names(extra) <- c("Alphaproteobacteria","Gammaproteobacteria")
paleta <- c(paleta,extra)
#Create new phyla dataset
Dat_phyla <- create_dataset(Tab = Tab,Map = Map)

Dat_phyla$Map$typebyTissue <- Dat_phyla$Map$typebyTissue %>%
  gsub(pattern = "Inoculum_Inoculum",replacement = "Inoculum") %>%
  factor(levels = c("Inoculum","AgarNoPlant","AgarPlant","Root","Shoot"))


res <- chibi.phylogram(Tab = Dat_phyla$Tab,Map = Dat_phyla$Map,facet_formula = "typebyTissue",
                       size_ticks_x = 0,size_strip_text = 35,size_axis_text = 25,
                       legend_proportion_size = 4,size_legend_text = 30,
                       size_axis_title = 0,font_family = "Arial",size_title_text = 35)
p <- res$p_mean + scale_fill_manual(values = paleta) +
  theme(aspect.ratio = 10) + theme(legend.position = "none")

oh.save.pdf(p = p, outdir = "../figures/",outname = "amplicon_4stresses_figure_syncom_phylogram_fraction.pdf")

