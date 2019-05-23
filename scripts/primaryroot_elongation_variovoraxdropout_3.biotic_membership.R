library(ohchibi)
library(paletteer)
library(extrafont)
loadfonts(device = "pdf") 
library(palettesPM)
library(ggtree)


setwd('/home/isai/Documents/results/github/variovoraxRGI/scripts/')

#Read metadata information for 201 unc isolates
Map <- read.table(file = "../rawdata/metadata_201_isolates_unc.tsv",
                header = F,sep = "\t",quote = "",comment.char = "")
Map <- Map[,c(1,4:8,3)] 
colnames(Map) <- c("taxon_oid","Phylum","Class","Order","Family","Genus","Name")

tree <- read.tree(file = "../rawdata/tree_201_uncisolates_47markers.rooted.newick")
tree <- ladderize(tree)

Map <- tree$tip.label %>% match(.,Map$taxon_oid) %>% Map[.,]



#Read the 185 strains in the inoculum
strains_185 <- read.table(file = "../rawdata/185_strains_innoculum.txt",header = F,sep = "\t") %$% 
V1  %>% as.character
meta <- read.table(file = "../rawdata/metadata_97useq_to_ids.tsv",
                 header = T,sep = "\t",quote = "",comment.char = "")
df_found <- which((strains_185 %in% meta$Freezer_Id)) %>% strains_185[.] %>%
match(.,meta$Freezer_Id) %>% meta[.,] %>% droplevels
df_found <- df_found[,c(4,2)]

df_cl144 <- data.frame(taxon_oid = Map$Name %>% grep(pattern = "CL144") %>%
                       Map$taxon_oid[.],Freezer_Id = "CL144")
df_cl72 <- data.frame(taxon_oid = Map$Name %>% grep(pattern = "CL72") %>%
                      Map$taxon_oid[.],Freezer_Id = "CL72")
df_cl9 <- data.frame(taxon_oid = Map$Name %>% grep(pattern = "UNCCL9") %>%
                     Map$taxon_oid[.],Freezer_Id = "CL9")
df_mf107 <- data.frame(taxon_oid = Map$Name %>% grep(pattern = "107") %>%
                       Map$taxon_oid[.],Freezer_Id = "MF107")

df_mf140 <- data.frame(taxon_oid = Map$Name %>% grep(pattern = "140Col") %>%
                       Map$taxon_oid[.],Freezer_Id = "MF140")
df_mf142 <- data.frame(taxon_oid = Map$Name %>% grep(pattern = "142") %>%
                       Map$taxon_oid[.],Freezer_Id = "MF142")

df_mf333 <- data.frame(taxon_oid = Map$Name %>% grep(pattern = "333") %>%
                       Map$taxon_oid[.],Freezer_Id = "MF333")

df_mf451 <- data.frame(taxon_oid = Map$Name %>% grep(pattern = "451") %>%
                       Map$taxon_oid[.],Freezer_Id = "MF451")

df_mf160 <- data.frame(taxon_oid = Map$Name %>% grep(pattern = "160") %>%
                       Map$taxon_oid[.],Freezer_Id = "MF160")

df_mf501a <- data.frame(taxon_oid = Map$Name %>% grep(pattern = "501") %>%
                        Map$taxon_oid[.],Freezer_Id = "MF501A")

df_mf72 <- data.frame (taxon_oid = "2510065092",Freezer_Id = "MF72")

df_all <- rbind(df_found,df_cl144,df_cl72,df_cl9,df_mf107,df_mf140,df_mf142,df_mf333,df_mf451,
              df_mf160,df_mf501a,df_mf72)

df_all <- rbind(df_all,data.frame(taxon_oid = "2636416056", Freezer_Id = "MF47"))
#which(!(strains_185 %in% df_all$Freezer_Id)) %>% strains_185[.]
chosen_taxons <- df_all$taxon_oid %>% na.omit %>% as.character 
todrop <- which(!(tree$tip.label %in% chosen_taxons)) %>% tree$tip.label[.]
tree <- ape::drop.tip(phy = tree,tip = todrop)

Map <- match(tree$tip.label,Map$taxon_oid) %>%
Map[.,] %>% droplevels

#There are two genomes that had inconsistent label with the phylogeny
#Mislabelled genomes
Map$Phylum <- Map$Phylum %>% as.character
Map$Order <- Map$Order %>% as.character
Map$Phylum[which(Map$taxon_oid == "2643221531")] <- "Actinobacteria"
Map$Order[which(Map$taxon_oid == "2643221531")] <- "Corynebacteriales"

Map$Phylum[which(Map$taxon_oid == "2639762635")] <- "Actinobacteria"
Map$Order[which(Map$taxon_oid == "2639762635")]  <- "Micrococcales"

Map$Order[which(Map$taxon_oid == "2548877166")]  <- "Micrococcales"


Map$Phylum <- Map$Phylum %>% factor
Map$Order <- Map$Order %>% factor

#Load palette of colors for phyla
mphyla <- palettesPM::pm.colors.phyla()
Map$PhylaColor <- match(Map$Phylum,names(mphyla)) %>%
mphyla[.] %>% as.character

Map$NLabel <- paste(Map$Genus,Map$taxon_oid,sep = " ")


#Read the useq to cluster df
df_useq <- read.table(file = "../rawdata/amplicon_4stresses_map_useq2cluster_abcd.tsv",header = T,sep = "\t")
df_end <- match(Map$taxon_oid,meta$taxon_oid) %>% meta$Useq[.] %>%
match(df_useq$Useq) %>% df_useq[.,] 
df_end <- cbind(Map,df_end)
df_end$Cluster <- df_end$Cluster %>% as.character
df_end$Cluster[which(df_end$taxon_oid == "2636416056")] <- "35"
df_end$Useq <- df_end$Useq %>% as.character
df_end$Useq[which(df_end$taxon_oid == "2636416056")] <- "Un"
df_end <- na.omit(df_end) 
todrop <- which(!(tree$tip.label %in% (df_end$taxon_oid %>% as.character))) %>%
  tree$tip.label[.]
tree <- ape::drop.tip(phy = tree,tip = todrop)



plot(tree,type = "fan",align.tip.label = T,show.tip.label = T,
    no.margin = TRUE,cex = 1,rotate.tree = 160,adj = 0.5,
    edge.width = 2,font = 4,label.offset = 0.02)
tiplabels(pch = 21,col = "#414141",
         bg = df_end$PhylaColor,cex = 2)



is_tip <- tree$edge[,2] <= length(tree$tip.label)
ordered_tips <- tree$edge[is_tip, 2]
mtips <- tree$tip.label[ordered_tips]

#Create figure
df_end$Freezer_Id <- match(df_end$taxon_oid,meta$taxon_oid) %>%
  meta$Freezer_Id[.] %>% gsub(pattern = "_.*",replacement = "") %>%
  as.character
df_end$Freezer_Id[72] <- "MF47"
df_end$NLabel <- paste(df_end$Genus,df_end$Freezer_Id,sep = " ")




tree$tip.label <- df_end$NLabel %>% as.character
cairo_pdf(filename = "../figures/primaryroot_elongation_variovoraxdropout_biotic_phylogeny.pdf",onefile = FALSE, fallback_resolution = 1200,
          width = 5, height = 15, family = "Arial", antialias = "default",
          pointsize = 12)
plot(tree,align.tip.label = T,show.tip.label = T,
     no.margin = TRUE,cex = 0.6,
     edge.width = 2,font = 4,label.offset = 0.01)
tiplabels(pch = 21,col = "black",
          bg = df_end$PhylaColor,cex = 1.2)
dev.off()
#Define a vector with members of the 35 syncom
memb35 <- c("CL69","MF33","MF374","MF376","MF161","CL21","MF138","MF79",
            "MF8","MF125","MF40","MF360","MF339","MF2","MF29","MF105",
            "MF50","MF370","MF109","MF181","CL14","MF47","CL18",
            "MF314","MF345","MF41","MF302","MF303","MF362","MF27",
            "MF299","MF327","MF273","MF136")
meta_35 <- match(memb35,meta$Freezer_Id) %>% meta[.,]
meta_35 <- meta_35[,c(2,4)]
meta_35$Freezer_Id <- meta_35$Freezer_Id %>% as.character
meta_35$taxon_oid <- meta_35$taxon_oid %>% as.character
meta_35$Freezer_Id[22] <- "MF47"
meta_35$taxon_oid[22] <- "2636416056"
meta_35 <- data.frame(taxon_oid = meta_35$taxon_oid, Cluster = "35Memb",value = rep(1,nrow(meta_35)))


#Prepare files for heatmap
df <- df_end[,c(1,11)]
df$value <- rep(1,nrow(df))
df$taxon_oid <- factor(df$taxon_oid,levels = mtips)
df$value <- df$value %>% factor
df <- df %>%
  #subset(Cluster != "A") %>%
  subset(Cluster != 35) %>% droplevels
df <- rbind(df,meta_35)
df$Cluster <- df$Cluster %>% factor(levels = c("A","B","C","D","35Memb"))

df$taxon_oid <- df$taxon_oid %>% factor(levels = mtips)
p <- ggplot(data = df,aes(Cluster,taxon_oid, fill = value)) + 
  facet_grid(.~Cluster,space = "free",scales = "free") + 
  geom_tile() + theme_ohchibi() + 
  theme(
    legend.position = "none",
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.major.y = element_blank(),
    strip.text.x = element_blank()
    ) +
  coord_cartesian(expand = F) +
  scale_fill_manual(values = "black")

#Create a bar 
df_end <- match(levels(df$taxon_oid) %>% rev,df_end$taxon_oid) %>%
  df_end[.,]
df_end$Genus
#Create a bar to align to the tree and put the genus name
bar<- c(rep("white",11),rep("black",25),rep("white",1),rep("black",7),
        rep("white",10),rep("black",1),rep("white",6),rep("black",1),
        rep("white",1),rep("black",13),rep("white",1),rep("black",3),
        #Leifsonia
        rep("white",11),rep("black",10),rep("white",1),rep("black",2),
        rep("white",2),rep("black",1),rep("white",2),rep("black",1),
        #Methylobacterium
        rep("white",10),rep("black",2),rep("white",3),rep("black",5),
        #Agrobacterium
        rep("white",6),rep("black",2),rep("white",6),rep("black",10),
        rep("white",1),rep("black",2),rep("white",6),rep("black",2),
        rep("white",10)
) 

paleta_blanco_negro <- c("white","black")
names(paleta_blanco_negro) <- c("white","black")
df_end$Bar <- bar
df_end$Type <- rep("Bar",nrow(df_end))
df_end$taxon_oid <- df_end$taxon_oid %>% factor(levels = df_end$taxon_oid %>% unique)
p_colorbar <- ggplot(df_end,aes(Type,taxon_oid %>% rev, fill = Bar,color = Bar)) +
  geom_tile(size = 0)  + theme_ohchibi(size_axis_text.y = 0,size_axis_title.x = 0,
                                       size_axis_title.y = 0,size_axis_text.x = 0,
                                       size_panel_border = 0.1,size_strip_text.x = 0) +
  theme(axis.ticks = element_blank(), 
        axis.text.x = element_blank(),
        panel.spacing = unit(0.25, "lines"),
        axis.text.y = element_blank(),
        strip.text.x = element_blank(),
        legend.position = "none",
        panel.grid.major.x = element_line(size = 0.4, color = "#D9D9D9"),
        panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank())+
  coord_cartesian(expand = FALSE) + scale_fill_manual(values = paleta_blanco_negro) +
  scale_color_manual(values = paleta_blanco_negro)
composition <- egg::ggarrange(p_colorbar,p,nrow = 1)
oh.save.pdf(p = composition,outname = "primaryroot_elongation_variovoraxdropout_biotic_membership.pdf",
            outdir = "../figures/",height = 15,width = 10)

