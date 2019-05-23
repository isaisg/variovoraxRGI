library(ohchibi)

setwd('/home/isai/Documents/results/github/variovoraxRGI/scripts/')

Map_norm <- read.table("../cleandata/primaryroot_elongation_monoassociation_normalized.tsv",header = T,
                       sep = "\t",quote = "",comment.char = "")

df_nb <- Map_norm %>% subset(Strain == "NB") %$% NormLength %>% quantile

#Read metadata information for 201 unc isolates
Map <- read.table(file = "../rawdata/metadata_201_isolates_unc.tsv",
                  header = F,sep = "\t",quote = "",comment.char = "")
Map <- Map[,c(1,4:8,3)] 
colnames(Map) <- c("taxon_oid","Phylum","Class","Order","Family","Genus","Name")

tree <- read.tree(file = "../rawdata/tree_201_uncisolates_47markers.rooted.newick")
tree <- ladderize(tree)


#Read the 185 Strains in the inoculum
Strains_185 <- read.table(file = "../rawdata/185_strains_innoculum.txt",header = F,sep = "\t") %$% 
  V1  %>% as.character
meta <- read.table(file = "../rawdata/metadata_97useq_to_ids.tsv",
                   header = T,sep = "\t",quote = "",comment.char = "")
df_found <- which((Strains_185 %in% meta$Freezer_Id)) %>% Strains_185[.] %>%
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


#which(!(Strains_185 %in% df_all$Freezer_Id)) %>% Strains_185[.]
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


meta_mapped <- match(Map$taxon_oid,meta$taxon_oid) %>% 
  meta[.,] %>% droplevels
Map$Freezer_Id <- meta_mapped$Freezer_Id %>% as.character
Map$Freezer_Id[which(is.na(Map$Freezer_Id))] <- "MF142"
Map$Freezer_Id <- Map$Freezer_Id %>% gsub(pattern = "_.*",replacement = "")
Map$NLabel <- paste(Map$Genus,Map$Freezer_Id,sep = " ")

#Intersect
chosen_ids <- intersect(Map$Freezer_Id %>% unique,Map_norm$Strain %>% unique)
Map <- which(Map$Freezer_Id %in% chosen_ids) %>% Map[.,] %>% droplevels
Map_norm <- (Map_norm$Strain %in% chosen_ids) %>% Map_norm[.,] %>% droplevels

#Drop tips
chosen <- Map$taxon_oid %>% unique %>% as.character
dropti <- which(!(tree$tip.label %in% chosen)) %>% tree$tip.label[.]
tree <- ape::drop.tip(phy = tree,tip = dropti)


plot(tree,align.tip.label = T,show.tip.label = T,
     no.margin = TRUE,cex = 1,rotate.tree = 160,adj = 0.5,
     edge.width = 2,font = 4,label.offset = 0.02)
tiplabels(pch = 21,col = "#414141",
          bg = Map$PhylaColor,cex = 2)

is_tip <- tree$edge[,2] <= length(tree$tip.label)
ordered_tips <- tree$edge[is_tip, 2]
mtips <- tree$tip.label[ordered_tips]


tree$tip.label <- match(tree$tip.label,Map$taxon_oid) %>% Map$NLabel[.] %>%
  as.character
cairo_pdf(filename = "../figures/primaryroot_elongation_monoassociation_phylogeny.pdf",onefile = FALSE, fallback_resolution = 1200,
          width = 5, height = 15, family = "Arial", antialias = "default",
          pointsize = 12)
plot(tree,align.tip.label = T,show.tip.label = T,
     no.margin = TRUE,cex = 0.6,rotate.tree = 160,adj = 0.5,
     edge.width = 2,font = 4,label.offset = 0.02)
tiplabels(pch = 21,col = "#414141",
          bg = Map$PhylaColor,cex = 1)
dev.off()

#Now create the boxplots
Map_norm$taxon_oid <- match(Map_norm$Strain,Map$Freezer_Id) %>% Map$taxon_oid[.]
Map_norm$taxon_oid <- Map_norm$taxon_oid %>% factor(levels = mtips)
df_seg <- data.frame(two = df_nb[2] %>% as.numeric, four = df_nb[4] %>% as.numeric)

#write.table(x = Map,file = "../cleandata/map_cleaned_174.tsv",append = F,quote = F,sep = "\t",row.names = F,col.names = T)

paleta <- palettesPM::pm.colors.phyla()
paleta <- which(names(paleta) %in% c("Actinobacteria","Bacteroidetes","Firmicutes")) %>%
  paleta[.]
paleta_proteo <- c("#CC873A","#FF6100")
names(paleta_proteo) <- c("Alphaproteobacteria","Gammaproteobacteria")
paleta <- c(paleta,paleta_proteo)
p <- ggplot(data = Map_norm,aes(taxon_oid,NormLength)) +
  geom_boxplot(aes(fill = Color),size = 0.5,outlier.size = NA,outlier.shape = NA) 
dat <- ggplot_build(p)$data[[1]]

p <- p + geom_segment(data=dat, aes(x=xmin, xend=xmax, 
                                    y=middle, yend=middle), colour="white", size=2) +
  geom_rect(data = df_seg,
            mapping = aes(ymin = two,ymax = four,xmin = 0,xmax = Inf),
            inherit.aes = F,color = NA,fill = "#56B4E9",alpha = 0.2) +
  geom_hline(yintercept = 3,size = 1,colour="black", linetype="dashed") +
  
  coord_flip() + theme_ohchibi() +
  theme(
    panel.grid.major.y = element_blank(),
    legend.position = "none",
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.y = element_blank()
  )   +
  scale_fill_manual(values = paleta,na.value = "#D9D9D9")


#Create_colorbar
Map$taxon_oid <- Map$taxon_oid %>% factor(levels = mtips %>% rev)
Map <- with(Map,order(taxon_oid)) %>%
  Map[.,]
bar<- c(rep("white",11),rep("black",25),rep("white",1),rep("black",7),
        rep("white",10),rep("black",1),rep("white",6),rep("black",1),
        rep("white",1),rep("black",13),rep("white",1),rep("black",3),
        rep("white",11),rep("black",10),rep("white",1),rep("black",2),
        rep("white",2),rep("black",1),rep("white",2),rep("black",1),
        rep("white",10),rep("black",2),rep("white",3),rep("black",5),
        rep("white",6),rep("black",2),rep("white",5),rep("black",10),
        rep("white",1),rep("black",2),rep("white",6),rep("black",2),
        rep("white",10)
) 

paleta_blanco_negro <- c("white","black")
names(paleta_blanco_negro) <- c("white","black")
Map$Bar <- bar
Map$Type <- rep("Bar",nrow(Map))
Map$taxon_oid <- Map$taxon_oid %>% factor(levels = Map$taxon_oid %>% as.character)

p_colorbar <- ggplot(Map,aes(Type,taxon_oid %>% rev, fill = Bar,color = Bar)) +
  geom_tile(size = 0)  + theme_ohchibi(size_axis_text.y = 0,size_axis_title.x = 0,
                                       size_axis_title.y = 0,size_axis_text.x = 0,
                                       size_panel_border = 0.1,size_strip_text.x = 0) +
  theme(axis.ticks = element_blank(), 
        axis.text.x = element_blank(),
        panel.spacing = unit(0.25, "lines"),
        axis.text.y = element_blank(),
        legend.position = "none",
        panel.grid.major.x = element_line(size = 0.4, color = "#D9D9D9"),
        panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank())+
  coord_cartesian(expand = FALSE) + scale_fill_manual(values = paleta_blanco_negro) +
  scale_color_manual(values = paleta_blanco_negro)



composition <- egg::ggarrange(p_colorbar,p,nrow = 1,widths = c(0.05,1))

### Test phylogenetic signal
df_sig <- aggregate(NormLength ~Strain,data = Map_norm,FUN = mean)
tree$tip.label <- tree$tip.label %>% gsub(pattern = ".* ",replacement = "")
df_sig <- match(tree$tip.label,df_sig$Strain) %>% df_sig[.,]

phylosig(tree = tree,df_sig$NormLength,method = "lambda",test = TRUE)

oh.save.pdf(p = composition,outname = "primaryroot_elongation_monoassociation_boxplots_phylogeny.pdf",outdir = "../figures/",
            height = 15,width = 15)
