library(ohchibi)
library(egg)

set.seed(130816)
setwd('/home/isai/Documents/results/github/variovoraxRGI/scripts/')

tree <- read.tree(file = "../rawdata/tree_124_isolates_burkholderia_56markers.rooted.newick")
meta <- read.table("../rawdata/metadata_124_isolates_burkholderia.tsv",header = T,
                   sep = "\t",quote = "",comment.char = "")


meta_vario <- meta$Genus %>% grep(pattern = "Variovorax") %>% 
  meta[.,] %>% droplevels

#Use acidovorax and cl11 to root
meta_outgroups <- meta$Genome.Name %>% grep(pattern = "CL11|Root219") %>%
  meta[.,] %>% droplevels

meta_sub <- rbind(meta_outgroups,meta_vario)

todrop <- which(!(tree$tip.label %in% (meta_sub$taxon_oid %>% as.character))) %>%
  tree$tip.label[.]
tree <- ape::drop.tip(phy = tree,todrop)
tree <- ape::ladderize(tree)


tree$tip.label <- match(tree$tip.label,meta_sub$taxon_oid) %>% meta_sub$Genome.Name[.] %>% 
  as.character
cairo_pdf(filename = "../figures/primaryroot_elongation_all_public_variovorax_phylogeny.pdf",onefile = FALSE, fallback_resolution = 1200,
           width = 5, height = 15, family = "Arial", antialias = "default",
           pointsize = 12)
plot(tree,align.tip.label = T,show.tip.label = ,
     no.margin = TRUE,cex = 0.6,
     edge.width = 2,font = 4,label.offset = 0.01)
tiplabels(pch = 21,col = "black",
          bg =meta_sub$Genus,cex = 1.5)
 dev.off()

is_tip <- tree$edge[,2] <= length(tree$tip.label)
ordered_tips <- tree$edge[is_tip, 2]
mtips <- tree$tip.label[ordered_tips]


#Create a structure to quantify which ones revert
meta_sub$CL28 <- rep("Untested",nrow(meta_sub))
meta_sub$CL28[meta_sub$Genome.Name %>% grep(pattern = "MF|CL14|110B",value = F)] <- "Yes"
meta_sub$CL28[meta_sub$Genome.Name %>% grep(pattern = "295")] <- "Untested"

#Read close relatives tested
reverters <- read.table(file = "../rawdata/data_primaryroot_elongation_cl28vsvariovoraxphylogeny_binary.txt") %$% V1 %>%
  as.character

pos <- meta_sub$Genome.Name %>% 
  grep(pattern = "YR266|PDC80|CF313|Root411|YR216|OK605|Root318D1|Root473|B4|YR750") 
meta_sub$CL28[pos] <- "Yes"
meta_sub$CL28 <- meta_sub$CL28 %>% factor(levels = c("Untested","No","Yes"))

meta_sub$CL28[meta_sub$Genome.Name %>% grep(pattern = "CL11|Root219")] <- "No"
meta_sub$Genome.Name <- meta_sub$Genome.Name  %>% factor(levels = mtips)
meta_sub$Type <- rep("Bar",nrow(meta_sub))
p_cl28 <- ggplot(data = meta_sub,aes(Type,Genome.Name, fill = CL28)) + geom_tile() +
  coord_cartesian(expand = F) +
  scale_fill_manual(values = c("grey","white","black")) +
  theme_ohchibi() +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "none"
  )

#Add another bar associated with presence of acc deaminase
positive_acc <- read.table(file = "../rawdata/data_K01505_distribution_burkholderia.txt") %$% V1
#Burkholderia CL11 and Acidovorax Root219 also have ACC deaminase
positive_acc <- c(2546825541,2643221611,positive_acc)
meta_sub$ACC <- rep("No",nrow(meta_sub))
meta_sub$ACC[which(meta_sub$taxon_oid %in% positive_acc)] <- "Yes"
meta_sub$ACC <- meta_sub$ACC %>% factor

p_acc <- ggplot(data = meta_sub,aes(Type,Genome.Name, fill = ACC)) + geom_tile() +
  coord_cartesian(expand = F) +
  scale_fill_manual(values = c("white","black")) +
  theme_ohchibi() +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "none"
  )

composition <- egg::ggarrange(p_cl28,p_acc, nrow = 1)
oh.save.pdf(p = composition,outname = "primaryroot_elongation_all_public_variovorax_bars.pdf",outdir = "../figures/",height = 15,width = 10)
