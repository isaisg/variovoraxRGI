library(ohchibi)
library(paletteer)
library(egg)
library(RColorBrewer)
library(extrafont)
loadfonts(device = "pdf") 
library(scales)
library(dendextend)
library(adephylo)


setwd('/home/isai/Documents/results/github/variovoraxRGI/scripts/')

set.seed(seed = 130816)

#Read the dataset 
Dat_amplicon <- readRDS(file = "../cleandata/dat_amplicon_4stresses_useq97.RDS")
Dat <- Dat_amplicon$RelativeAbundance

########### Definition of modules ################
Tab <- Dat$Tab
Tab <- Tab[grep(pattern = "OTU",x = rownames(Tab),invert = T),]

Tab_ra <- Tab
Tab_to_plot<-Tab_ra*1000
Tab_to_plot<-Tab_to_plot+1
Tab_to_plot<-log10(Tab_to_plot)

Tab_rowmeans<-Tab_ra/rowMeans(Tab_ra)
#This is a constant
Tab_rowmeans<-log10((Tab_rowmeans*1000)+1)


distfun <- function(x,method) as.dist(1-cor(t(x)))
hclustfun <- function(x,method) hclust(d = x,method = "ward.D2")

d <- hclustfun(distfun(Tab_rowmeans))

#Rotate C and D branches for legibility
label_order <- d$labels[d$order]
new_order <- c(label_order[14:44],label_order[1:13],label_order[45:94])
new_clust_samples <- dendextend::rotate(x = d,new_order) 
d_new <- new_clust_samples %>% as.dendrogram
d_new <- d_new %>% as.hclust
d_order <- rev(d_new$labels[d_new$order])



write.table(x = d_order,
            file = "../rawdata/amplicon_4stresses_order_strains_ra_heatmap.tsv",
            append = F,quote = F,sep = "\t",row.names = F,col.names = F)


# dendogram representation in ggplot2
dend <- as.dendrogram(d_new)
tree_dend <- as.phylo.dendrogram(d_new)


#Read phylogenetic tree and metadata
tree <- read.tree(file = "../rawdata/tree_201_uncisolates_47markers.rooted.newick")
meta <- read.table(file = "../rawdata/metadata_97useq_to_ids.tsv",
                   header = T,sep = "\t",quote = "",comment.char = "")

meta <- which(meta$Useq %in% label_order) %>% meta[.,] %>%
  droplevels
df_useq2taxon <- NULL
for(useq in levels(meta$Useq)){
  tids <- meta %>% subset(Useq == useq) %$% taxon_oid %>%
    na.omit %>% as.character %>% unique
  if(length(tids) ==0){
  }else if(length(tids) ==1){
    df_useq2taxon <- data.frame(USeq = useq,Representative = tids) %>%
      rbind(df_useq2taxon,.)
  }else{
    dfdist <- distRoot(tree,tids) %>% 
      data.frame(taxon_oid = names(.),Distance = .,row.names = NULL)
    dfdist <- with(dfdist,order(Distance)) %>% dfdist[.,]
    rep <- dfdist[1,1] %>% as.character
    df_useq2taxon <- data.frame(USeq = useq,Representative = rep) %>%
      rbind(df_useq2taxon,.)
  }
}
df_useq2taxon$Representative <- df_useq2taxon$Representative %>% as.character
df_useq2taxon$Representative[which(df_useq2taxon$USeq == "Sequence_86")] <- "2639762526"


#Change some labels 
taxa_name <- match(df_useq2taxon$Representative,meta$taxon_oid) %>% 
  meta$Genome_Classification[.] %>% as.character
df_tree_tax <- taxa_name %>% strsplit(split = "\\;") %>% unlist %>%
  gsub(pattern = "[k|p|c|o|f|g]__",replacement = "") %>%
  matrix(data = .,ncol = 6,byrow = T) %>% as.data.frame
colnames(df_tree_tax) <- c("Kingdom","Phylum","Class","Order","Family","Genus")
df_tree_tax$Representative  <- df_useq2taxon$Representative
df_useq2taxon <- merge(df_useq2taxon,df_tree_tax, by = "Representative")


mphyla <- palettesPM::pm.colors.phyla()
df_useq2taxon$PhylaColor <- match(df_useq2taxon$Phylum,names(mphyla)) %>%
  mphyla[.] %>% as.character


df_useq2taxon <- match(tree_dend$tip.label,df_useq2taxon$USeq) %>% df_useq2taxon[.,]

#Add a color for NA
df_useq2taxon$PhylaColor[which(is.na(df_useq2taxon$PhylaColor))] <- "#D9D9D9"

#Create another latyer where we differentiate the alpha proteobacteria from the gamma proteobacteria
df_useq2taxon$PhylaColorClass <- df_useq2taxon$PhylaColor

df_useq2taxon$PhylaColorClass[which(df_useq2taxon$Class == "Alphaproteobacteria")] <- "#CC873A"
df_useq2taxon$PhylaColorClass[which(df_useq2taxon$PhylaColorClass == "#FF8000")] <- "#FF6100"


cairo_pdf(filename = "../figures/amplicon_4stresses_figure_tree_heatmap.pdf",onefile = FALSE, fallback_resolution = 1200,
          width = 2, height = 15, family = "Arial", antialias = "default",
          pointsize = 12)
plot(tree_dend,align.tip.label = T,show.tip.label = F,adj = 1,
     no.margin = TRUE,cex = 1)
tiplabels(pch = 21,col = "#414141",
          bg = df_useq2taxon$PhylaColorClass,cex = 1.2)
dev.off()

#Compute family enrichments per module
map_useq2cluster <- read.table(file = "../rawdata/amplicon_4stresses_map_useq2cluster_abcd.tsv",header = T,sep = "\t")
colnames(map_useq2cluster)[1] <- "USeq"
df_enrichments <- merge(df_useq2taxon, map_useq2cluster, by = "USeq",all.x = TRUE)

#Test using hypergeometric test 
df_enrichments$Count <- rep(1,nrow(df_enrichments))
Tab_enrichments <- acast(data = df_enrichments %>% na.omit,formula = Cluster ~ Family,fun.aggregate = sum,
      value.var  = "Count") 
total_universe <- sum(Tab_enrichments)
Total <- colSums(Tab_enrichments)

Tab_enrichments <- rbind(Tab_enrichments,Total)
Total <- rowSums(Tab_enrichments)
Tab_enrichments <- cbind(Tab_enrichments,Total)

#COmpute enrichments using hypergeometric
df_overrep <- NULL
df_singletons <- NULL
for(i in 1:(ncol(Tab_enrichments)-1)){
  for(j in 1:(nrow(Tab_enrichments)-1)){
    val <-Tab_enrichments[j,i]
    if(val > 0){
          if(val > 1){
            module <- rownames(Tab_enrichments)[j]  
            family <- colnames(Tab_enrichments)[i]
            total_module <- Tab_enrichments[j,ncol(Tab_enrichments)]
            total_family <- Tab_enrichments[nrow(Tab_enrichments),i]
            failInPop <- total_universe - total_family
            pval <- phyper(val-1, total_family, failInPop, 
                           total_module, lower.tail= FALSE);
            temp <- data.frame(Module = module, Family = family, pvalue = pval)  
            df_overrep <- rbind(df_overrep,temp)
          }else{
            module <- rownames(Tab_enrichments)[j]  
            family <- colnames(Tab_enrichments)[i]
            temp <- data.frame(Module = module, Family = family)  
            df_singletons <- rbind(df_singletons,temp)
          }

    }
  }
}

df_overrep$padj <- df_overrep$pvalue %>% p.adjust(method = "fdr")
df_overrep %>% subset( padj <0.1)
######### Heatmap #################
#Read enrichment results
pthres <- 0.05
res <- read.table("../rawdata/dataS2_fraction_enrichments.tsv",
                  header = T,sep = "\t",quote = "",comment.char = "")
res <- res$Id %>% grep(pattern = "OTU",invert = T) %>% res[.,] %>%
  droplevels
#Do a general correction pf pvalues
res$padj <- res$pvalue %>% p.adjust(p = .,method = "fdr")
res$Id <- res$Id %>% factor(levels = rev(d_order))
res$Significance <- rep("NoSignificant",nrow(res))
res$Significance[which(res$padj < pthres)] <- "Significant"
res$Significance <- res$Significance %>% factor
res$Facet <- res$Contrast %>% as.character
res$Facet[res$Facet %>% grep(pattern = "Pi")] <- "Phosphate"
res$Facet[res$Facet %>% grep(pattern = "NaCl")] <- "Salinity"
res$Facet[res$Facet %>% grep(pattern = "pH")] <- "pH"
res$Facet[res$Facet %>% grep(pattern = "C$")] <- "Temperature"
res$Facet <- res$Facet %>% factor(levels = c("Phosphate","Salinity","pH","Temperature"))
res$log2FoldChange[which(res$log2FoldChange > 5)] <- 5
res$log2FoldChange[which(res$log2FoldChange < -5)] <- -5


##### Root vs Agar ####
res_sub <- res %>% subset(OverallContrast == "RootvsAgar") %>%
  droplevels
res_sub$Contrast <- res_sub$Contrast %>% factor(levels = c(
  "Root_0Pi_vs_AgarPlant_0Pi","Root_10Pi_vs_AgarPlant_10Pi",
  "Root_30Pi_vs_AgarPlant_30Pi","Root_50Pi_vs_AgarPlant_50Pi",
  "Root_100Pi_vs_AgarPlant_100Pi","Root_1000Pi_vs_AgarPlant_1000Pi",
  "Root_50NaCl_vs_AgarPlant_50NaCl","Root_100NaCl_vs_AgarPlant_100NaCl",
  "Root_150NaCl_vs_AgarPlant_150NaCl","Root_200NaCl_vs_AgarPlant_200NaCl",
  "Root_5.5pH_vs_AgarPlant_5.5pH","Root_7pH_vs_AgarPlant_7pH",
  "Root_8.2pH_vs_AgarPlant_8.2pH","Root_10C_vs_AgarPlant_10C",
  "Root_21C_vs_AgarPlant_21C","Root_31C_vs_AgarPlant_31C"
))

p_heatmap_root <- ggplot(data = res_sub,mapping = aes(x = Contrast, y = Id)) + 
  geom_raster(aes(fill = log2FoldChange)) +  
  geom_tile(aes(color = Significance),fill = '#00000000', size = 0.3,width = 0.85,height = 0.85) + 
  facet_grid(~Facet,space = "free",scales = "free") +
  scale_fill_paletteer_c(package = "pals",palette = "kovesi.diverging_bwr_55_98_c37",limits = c(-5,5),oob = squish) +
  scale_color_manual(values = c("#00000000","#414141"))+
  theme_ohchibi(size_axis_text.x = 0,size_axis_text.y = 0,size_axis_title.x = 0,
                size_axis_title.y = 0,size_strip_text.x = 0,
                legend_proportion_size = 0,size_title_text = 0,size_legend_text = 0,
                size_panel_border = 0.1,size_lines_panel = 0) +
  theme(aspect.ratio = 0.125,axis.ticks = element_blank(),
        panel.spacing = unit(0.2, "lines"),
        legend.position = "none")

##### Shoot vs Agar ####
res_sub <- res %>% subset(OverallContrast == "ShootvsAgar") %>%
  droplevels
res_sub$Contrast <- res_sub$Contrast %>% factor(levels = c(
  "Shoot_0Pi_vs_AgarPlant_0Pi","Shoot_10Pi_vs_AgarPlant_10Pi",
  "Shoot_30Pi_vs_AgarPlant_30Pi","Shoot_50Pi_vs_AgarPlant_50Pi",
  "Shoot_100Pi_vs_AgarPlant_100Pi","Shoot_1000Pi_vs_AgarPlant_1000Pi",
  "Shoot_50NaCl_vs_AgarPlant_50NaCl","Shoot_100NaCl_vs_AgarPlant_100NaCl",
  "Shoot_150NaCl_vs_AgarPlant_150NaCl","Shoot_200NaCl_vs_AgarPlant_200NaCl",
  "Shoot_5.5pH_vs_AgarPlant_5.5pH","Shoot_7pH_vs_AgarPlant_7pH",
  "Shoot_8.2pH_vs_AgarPlant_8.2pH","Shoot_10C_vs_AgarPlant_10C",
  "Shoot_21C_vs_AgarPlant_21C","Shoot_31C_vs_AgarPlant_31C"
))

p_heatmap_shoot <- ggplot(data = res_sub,mapping = aes(x = Contrast, y = Id)) + 
  geom_raster(aes(fill = log2FoldChange)) +  
  geom_tile(aes(color = Significance),fill = '#00000000', size = 0.3,width = 0.85,height = 0.85) + 
  facet_grid(~Facet,space = "free",scales = "free") +
  scale_fill_paletteer_c(package = "pals",palette = "kovesi.diverging_bwr_55_98_c37",limits = c(-5,5),oob = squish) +
  scale_color_manual(values = c("#00000000","#414141"))+
  theme_ohchibi(size_axis_text.x = 0,size_axis_text.y = 0,size_axis_title.x = 0,
                size_axis_title.y = 0,size_strip_text.x = 0,
                legend_proportion_size = 0,size_title_text = 0,size_legend_text = 0,
                size_panel_border = 0.1,size_lines_panel = 0) +
  theme(aspect.ratio = 0.125,axis.ticks = element_blank(),
        panel.spacing = unit(0.2, "lines"),
        legend.position = "none")

######## Root  vs Shoot ####
res_sub <- res %>% subset(OverallContrast == "RootvsShoot") %>%
  droplevels
res_sub$Contrast <- res_sub$Contrast %>% factor(levels = c(
  "Root_0Pi_vs_Shoot_0Pi","Root_10Pi_vs_Shoot_10Pi",
  "Root_30Pi_vs_Shoot_30Pi","Root_50Pi_vs_Shoot_50Pi",
  "Root_100Pi_vs_Shoot_100Pi","Root_1000Pi_vs_Shoot_1000Pi",
  "Root_50NaCl_vs_Shoot_50NaCl","Root_100NaCl_vs_Shoot_100NaCl",
  "Root_150NaCl_vs_Shoot_150NaCl","Root_200NaCl_vs_Shoot_200NaCl",
  "Root_5.5pH_vs_Shoot_5.5pH","Root_7pH_vs_Shoot_7pH",
  "Root_8.2pH_vs_Shoot_8.2pH","Root_10C_vs_Shoot_10C",
  "Root_21C_vs_Shoot_21C","Root_31C_vs_Shoot_31C"
  
))

p_heatmap_rs <- ggplot(data = res_sub,mapping = aes(x = Contrast, y = Id)) + 
  geom_raster(aes(fill = log2FoldChange)) +  
  geom_tile(aes(color = Significance),fill = '#00000000', size = 0.3,width = 0.85,height = 0.85) + 
  facet_grid(~Facet,space = "free",scales = "free") +
  scale_fill_paletteer_c(package = "pals",palette = "kovesi.diverging_bwr_55_98_c37",limits = c(-5,5),oob = squish) +
  scale_color_manual(values = c("#00000000","#414141"))+
  theme_ohchibi(size_axis_text.x = 0,size_axis_text.y = 0,size_axis_title.x = 0,
                size_axis_title.y = 0,size_strip_text.x = 0,
                legend_proportion_size = 0,size_title_text = 0,size_legend_text = 0,
                size_panel_border = 0.1,size_lines_panel = 0) +
  theme(aspect.ratio = 0.125,axis.ticks = element_blank(),
        panel.spacing = unit(0.2, "lines"),
        legend.position = "none")


########## Panels of abundance in the three fractions #######

Dat_rab <- Dat_amplicon$RelativeAbundance
melted_rab <- Dat_rab$Tab %>% melt
colnames(melted_rab) <- c("Id","ID_Matrix","RA")
melted_rab <- merge(melted_rab,Dat_rab$Map, by = "ID_Matrix")

#Create data frame with mean RA 
df_rab <- dcast(data = melted_rab,formula = Id~typebyTissue,
                fun.aggregate = mean,value.var = "RA") %>% melt()
colnames(df_rab) <- c("Id","typebyTissue","MeanRA")

df_rab <-  which(df_rab$Id %in% levels(res$Id)) %>%
  df_rab[.,] %>% droplevels()
df_rab$Id <- df_rab$Id %>% factor(levels = res$Id %>% levels)
df_rab <- df_rab %>% subset(typebyTissue == "AgarPlant" | 
                              typebyTissue == "Root" | typebyTissue == "Shoot") %>%
  droplevels

#Create abundance bargraph
p_abundance <- ggplot(data = df_rab,mapping = aes(x= Id, y = (MeanRA)  )) + 
  geom_bar(stat = "identity") +facet_grid(~typebyTissue) + 
  coord_flip() + theme_ohchibi(size_axis_text.y = 0,size_axis_title.x = 0,
                               size_axis_title.y = 0,size_axis_text.x = 2.5,size_panel_border = 0.1,size_strip_text.x = 0) +
  theme(axis.ticks = element_blank(), 
        panel.spacing = unit(0.1, "lines"),
        legend.position = "none",aspect.ratio = 0.5,
        panel.grid.major.x = element_line(size = 0.3, color = "#D9D9D9"),
        panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank()) +
  scale_y_sqrt(breaks = c(0.01,0.05,0.1,0.2,0.3))  +
  geom_hline(yintercept = 0.01,color = "firebrick",size = 0.3) 



########## Phylogenetic signal ###################

#Read phylogenetic tree
tree <- read.tree(file = "../rawdata/tree_201_uncisolates_47markers.rooted.newick")
meta <- read.table(file = "../rawdata/metadata_97useq_to_ids.tsv",
                   header = T,sep = "\t",quote = "",comment.char = "")

ids <- res$Id %>% as.character %>% unique
meta <- which(meta$Useq %in% ids) %>% meta[.,] %>%
  droplevels
df_useq2taxon <- NULL
for(useq in levels(meta$Useq)){
  tids <- meta %>% subset(Useq == useq) %$% taxon_oid %>%
    na.omit %>% as.character %>% unique
  if(length(tids) ==0){
  }else if(length(tids) ==1){
    df_useq2taxon <- data.frame(USeq = useq,Representative = tids) %>%
      rbind(df_useq2taxon,.)
  }else{
    dfdist <- distRoot(tree,tids) %>% 
      data.frame(taxon_oid = names(.),Distance = .,row.names = NULL)
    dfdist <- with(dfdist,order(Distance)) %>% dfdist[.,]
    rep <- dfdist[1,1] %>% as.character
    df_useq2taxon <- data.frame(USeq = useq,Representative = rep) %>%
      rbind(df_useq2taxon,.)
  }
}

#Subset the tree
todrop <- which(!(tree$tip.label %in% df_useq2taxon$Representative)) %>%  
  tree$tip.label[.]
tree <- ape::drop.tip(phy = tree,tip = todrop)

tree$tip.label <- match(tree$tip.label,df_useq2taxon$Representative) %>% 
  df_useq2taxon$USeq[.] %>% as.character

#Create a matrix with log2FoldChanges
Tab_ps <- acast(data = res,formula =Id ~ Contrast ,value.var = "log2FoldChange") 
Tab_ps <- which(rownames(Tab_ps) %in% tree$tip.label) %>% Tab_ps[.,]
df_ps <- apply(X = Tab_ps,MARGIN = 2,
               FUN = function(x)phylosig(tree = tree,x,method = "lambda",test = TRUE) %>% unlist) %>%
  melt
colnames(df_ps) <- c("Type","Contrast","value")

#Subset pvalues
df_ps_p <- df_ps %>% subset(Type == "P") 
#Correct for multiple testing
df_ps_p$value <- df_ps_p$value %>% p.adjust(method = "fdr")
df_ps_p$Type <- df_ps_p$Type %>% as.character %>% 
  gsub(pattern = "P",replacement = "padj")
df_ps <- rbind(df_ps,df_ps_p )
df_ps <- colnames(res) %>% grep(pattern = "Contrast") %>% res[,.] %>%
  unique %>% merge(df_ps,., by = "Contrast")

#Subset padj
df_ps_p <- df_ps %>% subset(Type == "padj" & OverallContrast == "RootvsAgar") %>% droplevels
df_ps_p$lvalue <- -log10(df_ps_p$value)
df_ps_p$Facet <- rep("Temperature",nrow(df_ps_p))
df_ps_p$Facet[(df_ps_p$Contrast %>% grep(pattern = "Pi"))] <- "Phosphate"
df_ps_p$Facet[(df_ps_p$Contrast %>% grep(pattern = "NaCl"))] <- "Salinity"
df_ps_p$Facet[(df_ps_p$Contrast %>% grep(pattern = "pH"))] <- "pH"
df_ps_p$Facet <- df_ps_p$Facet %>% factor(levels = c("Phosphate","Salinity","pH","Temperature"))

p_signal_root <- df_ps_p %>% ggplot(data = .,mapping = aes(x = Contrast, y = lvalue)) + 
  geom_bar(stat = "identity") + facet_grid(~Facet,scales = "free",space = "free") +
  ylab(label = "-log10(padj)") +
  theme_ohchibi(size_axis_text.x = 0,size_axis_title.y =10,size_axis_title.x = 0,
                size_strip_text.x = 0,size_axis_text.y = 6,size_panel_border = 0.1) +
  theme(aspect.ratio = 0.5,axis.ticks.x = element_blank(),
        panel.spacing = unit(0.1, "lines"),
        legend.position = "none") + 
  geom_hline(yintercept = 1.30103,color = "firebrick", size = 0.5,linetype = "dashed") + 
 scale_y_reverse(limits = c(12,0)) 

df_ps_p <- df_ps %>% subset(Type == "padj" & OverallContrast == "ShootvsAgar") %>% droplevels
df_ps_p$lvalue <- -log10(df_ps_p$value)
df_ps_p$Facet <- rep("Temperature",nrow(df_ps_p))
df_ps_p$Facet[(df_ps_p$Contrast %>% grep(pattern = "Pi"))] <- "Phosphate"
df_ps_p$Facet[(df_ps_p$Contrast %>% grep(pattern = "NaCl"))] <- "Salinity"
df_ps_p$Facet[(df_ps_p$Contrast %>% grep(pattern = "pH"))] <- "pH"
df_ps_p$Facet <- df_ps_p$Facet %>% factor(levels = c("Phosphate","Salinity","pH","Temperature"))


p_signal_shoot <- df_ps_p %>% ggplot(data = .,mapping = aes(x = Contrast, y = lvalue)) + 
  geom_bar(stat = "identity") + facet_grid(~Facet,scales = "free",space = "free") +
  ylab(label = "-log10(padj)") +
  theme_ohchibi(size_axis_text.x = 0,size_axis_title.y =10,size_axis_title.x = 0,
                size_strip_text.x = 0,size_axis_text.y = 6,size_panel_border = 0.1) +
  theme(aspect.ratio = 0.5,axis.ticks.x = element_blank(),
        panel.spacing = unit(0.1, "lines"),
        legend.position = "none")  +
  geom_hline(yintercept = 1.30103,color = "red", size = 0.5,linetype = "dashed") +
  scale_y_reverse(limits = c(12,0)) 


#Create a heder with the descriptions of each column
res_header <- res$Contrast %>% grep(pattern = "Root.*AgarPlant") %>% res[.,] %>% droplevels
colores_phosphate <- c("#ffffe5","#f7fcb9","#addd8e","#41ab5d","#238443","#005a32")
colores_salinity <- c("#fde0dd","#fa9fb5","#dd3497", "#7a0177")
colores_temp <- c("#FFC88C","#FFA33F","#7F6446")
colores_ph <- c("#C3C1EA","#7E78E5","#555466")
colores_header <- c(colores_phosphate,colores_salinity,colores_ph,colores_temp)


df_header <- data.frame(Contrast = res_header$Contrast %>% levels,
                        Type = rep("Header",16),
                        Condition = res_header$Contrast %>% levels %>% gsub(pattern = "Root|Shoot|_vs_.*",replacement = "") %>%
                          gsub(pattern = "_",replacement = ""))
df_header$Contrast <- df_header$Contrast %>% factor(levels = c("Root_0Pi_vs_AgarPlant_0Pi","Root_10Pi_vs_AgarPlant_10Pi",
                                                               "Root_30Pi_vs_AgarPlant_30Pi","Root_50Pi_vs_AgarPlant_50Pi",
                                                               "Root_100Pi_vs_AgarPlant_100Pi","Root_1000Pi_vs_AgarPlant_1000Pi",
                                                               "Root_50NaCl_vs_AgarPlant_50NaCl","Root_100NaCl_vs_AgarPlant_100NaCl",
                                                               "Root_150NaCl_vs_AgarPlant_150NaCl","Root_200NaCl_vs_AgarPlant_200NaCl",
                                                               "Root_5.5pH_vs_AgarPlant_5.5pH","Root_7pH_vs_AgarPlant_7pH",
                                                               "Root_8.2pH_vs_AgarPlant_8.2pH","Root_10C_vs_AgarPlant_10C",
                                                               "Root_21C_vs_AgarPlant_21C","Root_31C_vs_AgarPlant_31C"))
df_header$Condition <- df_header$Condition %>% factor(levels = c("0Pi","10Pi","30Pi","50Pi","100Pi","1000Pi",
                                                                 "50NaCl","100NaCl","150NaCl","200NaCl",
                                                                 "5.5pH","7pH","8.2pH","10C","21C","31C"))
df_header <- match(levels(df_header$Condition),df_header$Condition) %>% df_header[.,]

df_header$Facet <- c(rep("Phosphate",6),rep("Salinity",4),rep("pH",3),rep("Temperature",3)) %>%
  factor(levels = c("Phosphate","Salinity","pH","Temperature"))
names(colores_header) <- df_header$Condition %>% levels

p_header <- ggplot(data = df_header,mapping = aes(Contrast, y = Type, fill = Condition)) +
  geom_tile() + facet_grid(~Facet,scales = "free",space = "free") +
  scale_fill_manual(values = colores_header) +
  theme_ohchibi(size_axis_text.x = 0,size_axis_text.y = 0,
                size_axis_title.x = 0,size_axis_title.y = 0,legend_proportion_size = 0,
                size_title_text = 0,size_legend_text = 0,
                size_strip_text.x = 0,size_panel_border = 0.1,size_lines_panel = 0)+
  theme(axis.ticks = element_blank(),
        panel.spacing = unit(0.5, "lines"),
        aspect.ratio = 2.5)

#oh.save.pdf(p = p_header,outname =  "legends_header_heatmap.pdf")
p_header <- p_header + theme(legend.position = "none")

##Create bar to name the genus
map_useq2cluster <- read.table(file = "../rawdata/amplicon_4stresses_map_useq2cluster_abcd.tsv",header = T,sep = "\t")
df_useq2cluster <- data.frame(Useq = levels(df_rab$Id))
df_useq2cluster <- merge(df_useq2cluster, map_useq2cluster, by = "Useq")
df_useq2cluster$Useq <- df_useq2cluster$Useq %>% factor(levels = levels(df_rab$Id))
df_useq2cluster <- match(levels(df_useq2cluster$Useq),df_useq2cluster$Useq) %>% 
  df_useq2cluster[.,]
df_useq2cluster$Cluster %>% table
df_useq2cluster$Color <- c(rep("black",18),rep("white",32),rep("black",13),rep("white",31)) %>% rev

paleta_blanco_negro <- c("white","black")
names(paleta_blanco_negro) <- c("white","black")
df_useq2cluster$Type <- rep("Bar",nrow(df_useq2cluster))

p_colorbar <- ggplot(df_useq2cluster,aes(Type,Useq, fill = Color,color = Color)) +
  geom_tile(size = 0)  + theme_ohchibi(size_axis_text.y = 0,size_axis_title.x = 0,
                                       size_axis_title.y = 0,size_axis_text.x = 0,
                                       size_panel_border = 0.1,size_strip_text.x = 0) +
  theme(axis.ticks = element_blank(), axis.text.x = element_blank(),
        panel.spacing = unit(0.25, "lines"),axis.text.y = element_blank(),
        legend.position = "none",
        panel.grid.major.x = element_line(size = 0.4, color = "#D9D9D9"),
        panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank())+
  coord_cartesian(expand = FALSE) + scale_fill_manual(values = paleta_blanco_negro) +
  scale_color_manual(values = paleta_blanco_negro)



######### Construct composite figure ##########
p_blank <- ggplot() + theme_void() + 
  theme(legend.position = "none")

composition <- egg::ggarrange(p_header,p_header,p_blank,
                              p_heatmap_root,p_heatmap_shoot,p_abundance,
                              p_signal_root,p_signal_shoot,p_blank,nrow = 3,ncol = 3,
                              widths = c(1.5,1.5,1),heights = c(0.05,1,0.3),debug = F,byrow = TRUE)


oh.save.pdf(p = composition,outname = "amplicon_4stresses_composition_heatmap_4stresses.pdf",outdir = "../figures/",
            height = 20,width = 15)



