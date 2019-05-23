library(ohchibi)
library(emmeans)
library(ape)
library(ggdendro)
library(paletteer)
library(palettesPM)
library(egg)
library(scales)

setwd('/home/isai/Documents/results/github/variovoraxRGI/scripts/')
set.seed(130816)

dir.create("../cleandata/")
dir.create("../figures//")

Map_sub <- read.table(file = "../rawdata/dataS6_primaryroot_elongation_variovoraxburkholderiavsRGIs.tsv",
                      header = T,sep = "\t")


#Test with ANOVA 
Res <- NULL
for(st in levels(Map_sub$InhibitedStrain)){
  Map_in <- subset(Map_sub,InhibitedStrain == st)
  rdf <- aov(formula = MainRoot ~ InhibitorStrain + Rep,data =Map_in ) %>% emmeans(specs = "InhibitorStrain") %>%
    pairs(reverse = TRUE,adjust = "none") %>% data.frame 
  rdf$contrast <- rdf %$% contrast %>% as.character %>% 
    gsub(pattern = " ",replacement = "") %>% factor
  rdf$InhibitedStrain <- rep(st,nrow(rdf))
  rdf$Rep <- rep("Yes",nrow(rdf))
  Res <- rbind(Res, rdf)
}


Res <- Res %$% contrast %>% grep(pattern = "NB") %>% Res[.,] %>% droplevels
Res$p.adj <- Res %$% p.value %>% p.adjust(method = "fdr")
Res$InhibitorStrain <- Res %$% contrast %>% gsub(pattern = ".*-",replacement = "") %>% factor


#Map the data to the tree
meta <- read.table(file = "../rawdata/metadata_97useq_to_ids.tsv",
                   header = T,sep = "\t")
meta$Freezer_Id <- meta %$% Freezer_Id %>% as.character %>% 
  gsub(pattern = "\\_.*",replacement = "")
meta <- meta[,c(2,4,5,6)] %>% unique

#Identify taxon_oid of inhibitors
meta_chosen <- Map_sub$InhibitedStrain %>% levels %>%
  grep(pattern = "NB",invert = T,value = T) %>%
  match(meta$Freezer_Id) %>% meta[.,] %>% droplevels

#Load phylogenetic tree
tree <- read.tree(file = "../rawdata/tree_201_uncisolates_47markers.rooted.newick")
to_remove <- tree$tip.label[which(!(tree$tip.label %in% meta_chosen$taxon_oid))]
tree <- ape::drop.tip(phy = tree,tip = to_remove)

tree <- tree %>% ladderize
meta <- match(tree$tip.label,meta$taxon_oid) %>% meta[.,] %>% droplevels
mphyla <- palettesPM::pm.colors.phyla()

taxa_name <- meta$Genome_Classification %>% as.character
df_tree_tax <- taxa_name %>% strsplit(split = "\\;") %>% unlist %>%
  gsub(pattern = "[k|p|c|o|f|g]__",replacement = "") %>%
  matrix(data = .,ncol = 6,byrow = T) %>% as.data.frame
colnames(df_tree_tax) <- c("Kingdom","Phylum","Class","Order","Family","Genus")
df_tree_tax$taxon_oid <- meta$taxon_oid
meta <- merge(meta,df_tree_tax, by = "taxon_oid")

meta$PhylaColor <- match(meta$Phylum,names(mphyla)) %>% mphyla[.] %>% as.character
meta$ClassColor <- meta$PhylaColor
meta$ClassColor[which(meta$ClassColor == "#FF8000")] <- "#FF6100"
meta$ClassColor[which(meta$Class == "Alphaproteobacteria")] <- "#CC873A"
meta <- match(tree$tip.label,meta$taxon_oid) %>% meta[.,]

tree$tip.label <- meta$Freezer_Id %>% as.character


cairo_pdf(filename = "../figures/primaryroot_elongation_tripartite_variovoraxburkholderiavsRGIs_phylogeny.pdf",onefile = FALSE, fallback_resolution = 1200,
          width = 15, height = 2, family = "Arial", antialias = "default",
          pointsize = 12)
plot(tree,align.tip.label = T,show.tip.label = F,adj = 1,
     no.margin = TRUE,cex = 1,direction = "downwards")
tiplabels(pch = 21,col = "#414141",
          bg = meta$ClassColor,cex = 3)
dev.off()

is_tip <- tree$edge[,2] <= length(tree$tip.label)
ordered_tips <- tree$edge[is_tip, 2]  %>% tree$tip.label[.] 

meta <- match(ordered_tips,meta$Freezer_Id) %>% meta[.,]



#
Map_sub$InhibitedStrain <- Map_sub %$% InhibitedStrain %>% factor(levels = meta$Freezer_Id)
Map_sub$InhibitorStrain <- factor(Map_sub$InhibitorStrain,levels = c("NB","CL14","MF160","CL11"))

#Add the significance
Map_sub$Comparison <- paste(Map_sub$InhibitedStrain,Map_sub$InhibitorStrain,sep="-")
Res$Comparison <- paste(Res$InhibitedStrain,Res$InhibitorStrain,sep = "-")
Res_tomerge <- Res[match(Map_sub$Comparison,Res$Comparison),]
Map_sub$p.adj <- Res_tomerge$p.adj
Map_sub$p.adj[is.na(Map_sub$p.adj)] <- 1
Map_sub$estimate <- Res_tomerge$estimate
Map_sub$estimate[is.na(Map_sub$estimate)] <- 0
Map_sub$Significance <- rep("NoSignificant",nrow(Map_sub))
num_sig <- subset(Map_sub,p.adj < 0.05 & abs(estimate)>=2) %>% rownames
Map_sub$Significance[Map_sub %>% rownames %in% num_sig %>% which] <- "Significant"
Map_sub$Significance <- Map_sub$Significance %>% factor

#https://stackoverflow.com/questions/42641106/add-a-second-geom-tile-layer-in-ggplot
#parula
#isol
#linearl

Map_sub$InhibitedStrain <- Map_sub$InhibitedStrain  %>% 
  factor(levels = Map_sub$InhibitedStrain %>% levels %>% rev)
Map_sub$InhibitorStrain <- Map_sub$InhibitorStrain %>% 
  factor(levels = c("NB","CL11","CL14","MF160"))
p_tile <- ggplot(data = Map_sub,aes(x=InhibitedStrain,y=InhibitorStrain)) +
  geom_raster(aes(fill = MainRoot)) +
  geom_tile(aes(color = Significance),fill = '#00000000', size = 3,width = 0.85,height = 0.85) +
  scale_fill_paletteer_c(package = "pals",palette = "parula",limits = c(0,8.2),oob = squish) +
  scale_color_manual(values = c("#00000000","#414141"))+
  xlab(label = "Short Root Inhibitors") + ylab(label = "Short Root Inducers")+
  theme(
    aspect.ratio = 0.255,
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    #axis.ticks.length = element_blank()
    axis.text.y = element_text(family = "Arial", 
                               face = "bold", size = 12, colour = "#414141"),
    axis.text.x = element_text(family = "Arial", 
                               face = "bold", size = 12, colour = "#414141"),
    axis.title.x = element_text(family = "Arial", face = "plain", 
                                size = 15, colour = "#414141"),
    axis.title.y = element_text(family = "Arial", face = "plain", 
                                size = 15, colour = "#414141"),
    axis.ticks = element_blank()
    
  )

oh.save.pdf(p = p_tile,outname = "primaryroot_elongation_tripartite_variovoraxburkholderiavsRGIs_heatmap_legends.pdf",
            outdir = "../figures/")
p_tile <- p_tile + theme(legend.position = "none")

oh.save.pdf(p = p_tile,outname = "primaryroot_elongation_tripartite_variovoraxburkholderiavsRGIs_heatmap_nolegends.pdf",
            outdir = "../figures/",width = 25)


saveRDS(object = Map_sub,file = "../cleandata/dataset_tripartite_variovoraxburkholderiavsRGIs.RDS")
