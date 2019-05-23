library(ohchibi)
library(multcomp)

set.seed(130816)
setwd('/home/isai/Documents/results/github/variovoraxRGI/scripts/')

tree <- read.tree(file = "../rawdata/tree_124_isolates_burkholderia_56markers.rooted.newick")
meta <- read.table("../rawdata/metadata_124_isolates_burkholderia.tsv",header = T,
                   sep = "\t",quote = "",comment.char = "")

#Read experiment of relatives
df <- read.table("../rawdata/dataS8_primaryroot_elongation_cl28vsvariovoraxphylogeny.tsv",sep = "\t",header = T)


#Read hormonde data
df_hormones <- read.table(file = "../rawdata/dataS13_primaryroot_elongation_hormone_MAMP.tsv",
                          header = T,sep = "\t",quote = "",comment.char = "",check.names = F)

df_hormones <- df_hormones %>% subset(Treatment == "CL28" | Treatment== "NB") %>% droplevels


df_hormones <- data.frame(Length = df_hormones$Length,Strain1=df_hormones$Strain,Strain2 = df_hormones$Treatment)
df_hormones$Exp <- rep("Hormones",nrow(df_hormones))
df$Exp <- rep("Relatives",nrow(df))

merged <- rbind(df,df_hormones)

df <- merged



meta$Id <- meta$Genome.Name %>% gsub(pattern = ".* ",replacement = "")
meta$Id <- meta$Id %>% gsub(pattern = "CF079",replacement = "CF79") %>%
  gsub(pattern = "160MFSha2.1",replacement = "MF160")

df <- df %>% subset(Strain1 != "OF79") %>%
  subset(Strain1 != "CF79") %>%
  droplevels
merged <- merged%>% subset(Strain1 != "OF79") %>% 
  subset(Strain1 != "CF79") %>%
  droplevels

ids <- df$Strain1 %>% levels 
#Add the unc variovorax isolates
unc_vario <- c(2517572245,2519899633,2643221492,2643221497,2643221498,2643221499,2643221501,2643221501)
which(!(ids %in% meta$Id)) %>% ids[.]

todrop <- which(!(meta$Id %in% ids)) %>%
  meta$taxon_oid[.] %>% as.character
todrop <- which(!(todrop %in% unc_vario)) %>% todrop[.]
tree <- ape::drop.tip(phy = tree,tip = todrop)
meta<- match(tree$tip.label,meta$taxon_oid) %>%
  meta[.,] %>% droplevels
tree <- ladderize(tree)
# cairo_pdf(filename = "../figures/phylogeny_boxplots_monoassociation.pdf",onefile = FALSE, fallback_resolution = 1200,
#           width = 5, height = 15, family = "Arial", antialias = "default",
#           pointsize = 12)
plot(tree,align.tip.label = T,show.tip.label = T,
     no.margin = TRUE,cex = 0.6,rotate.tree = 160,adj = 0.5,
     edge.width = 2,font = 4,label.offset = 0.02)
tiplabels(pch = 21,col = "#414141",
          bg =meta$Genus,cex = 1)
# dev.off()

is_tip <- tree$edge[,2] <= length(tree$tip.label)
ordered_tips <- tree$edge[is_tip, 2]
mtips <- tree$tip.label[ordered_tips]

 cairo_pdf(filename = "../figures/primaryroot_elongation_all_public_variovorax_phylogenyvscl28.pdf",onefile = FALSE, fallback_resolution = 1200,
           width = 5, height = 15, family = "Arial", antialias = "default",
           pointsize = 12)
meta$Nom <- paste0(meta$Genus," ",meta$Id)
tree$tip.label <- meta$Nom
plot(tree,align.tip.label = T,show.tip.label = T,
     no.margin = TRUE,cex = 0.6,
     edge.width = 2,font = 4,label.offset = 0.001)
tiplabels(pch = 21,col = "#414141",
          bg =meta$Genus,cex = 3)
dev.off()


df_nb <- df %>% subset(Strain1 == "NB") %>% droplevels

df <- df %>% subset(Strain1 != "NB") %>% droplevels
df$taxon_oid <- match(df$Strain1,meta$Id) %>% meta$taxon_oid[.] %>% as.character
df$Genus <- match(df$Strain1,meta$Id) %>% meta$Genus[.] %>% as.character
df$taxon_oid <- df$taxon_oid  %>% factor(levels = mtips)
df$Strain2 <- df$Strain2 %>% factor(levels = c("NB","CL28","MF48"))

df_segment <- NULL
for(lev in df_nb$Strain2 %>% levels){
  temp <- df_nb %>% subset(Strain2 == lev) %>% 
    droplevels %$% Length %>% quantile
  temp <- data.frame(Strain2 = lev,two = temp[2] %>% as.numeric, four = temp[4] %>% as.numeric)
  df_segment <- rbind(df_segment,temp)
}

#Do testing
merged$Count <- rep(1,nrow(merged))
temp <- acast(data = merged,formula = Exp~Strain1,fun.aggregate = sum,
      value.var = "Count") %>% 
  apply(X = .,MARGIN = 2,FUN = function(x)which(x>0) %>% length)
temp <- which(temp == 2) %>% temp[.] %>% names

#Test the one with two reps
df_two <- which(merged$Strain1 %in% temp) %>%
  merged[.,] %>%
  subset(Strain2 != "MF48") %>% droplevels
df_two$Strain1 <- df_two$Strain1 %>% relevel(ref = "NB")
df_stat <- NULL
for(st in levels(df_two$Strain2)){
  df_temp <- df_two %>% subset(Strain2 == st) %>%
    droplevels
  m1 <- aov(data = df_temp,formula = Length ~ Strain1 + Exp)
  Dunnet <- glht(m1, linfct=mcp(Strain1="Dunnett"))
  sum <- summary(Dunnet)
  noms <- coefficients(sum) %>% names %>% gsub(pattern = " .*",replacement = "")
  pvals <- sum$test$pvalues
  df_res <- data.frame(Strain1 = noms,Strain2 = rep(st,length(noms)),
             pvalue = pvals)
  df_stat <- rbind(df_stat,df_res)
}

#Test the ones with one rep in MF48
df_sub <- merged %>% subset(Strain2 == "MF48") %>% droplevels
df_sub$Strain1  <- df_sub$Strain1 %>% relevel(ref = "NB")
for(st in levels(df_sub$Strain2)){
  df_temp <- df_sub %>% subset(Strain2 == st) %>%
    droplevels
  m1 <- aov(data = df_temp,formula = Length ~ Strain1 )
  Dunnet <- glht(m1, linfct=mcp(Strain1="Dunnett"))
  sum <- summary(Dunnet)
  noms <- coefficients(sum) %>% names %>% gsub(pattern = " .*",replacement = "")
  pvals <- sum$test$pvalues
  df_res <- data.frame(Strain1 = noms,Strain2 = rep(st,length(noms)),
                       pvalue = pvals)
  df_stat <- rbind(df_stat,df_res)
}

df_sub <- merged %>% subset(Strain2 != "MF48") %>% droplevels
temp <- temp %>% grep(pattern = "NB",value = T,invert = T)
df_sub <- which(!(df_sub$Strain1 %in% temp)) %>% 
  df_sub[.,] %>% droplevels
df_sub$Strain1  <- df_sub$Strain1 %>% relevel(ref = "NB")
for(st in levels(df_sub$Strain2)){
  df_temp <- df_sub %>% subset(Strain2 == st) %>%
    droplevels
  m1 <- aov(data = df_temp,formula = Length ~ Strain1 )
  Dunnet <- glht(m1, linfct=mcp(Strain1="Dunnett"))
  sum <- summary(Dunnet)
  noms <- coefficients(sum) %>% names %>% gsub(pattern = " .*",replacement = "")
  pvals <- sum$test$pvalues
  df_res <- data.frame(Strain1 = noms,Strain2 = rep(st,length(noms)),
                       pvalue = pvals)
  df_stat <- rbind(df_stat,df_res)
}

#Adjust all the pvalus
df_stat$padj <- df_stat$pvalue %>% p.adjust(method = "fdr")
#df_stat$padj <- df_stat$pvalue 

df_stat$Label <- rep("NS",nrow(df_stat))
df_stat$Label[which(df_stat$padj < 0.05 )] <- "*"
df_stat$Label[which(df_stat$padj < 0.01 )] <- "**"
df_stat$Label[which(df_stat$padj < 0.001 )] <- "***"
df_stat$taxon_oid <- match(df_stat$Strain1,df$Strain1) %>%
  df$taxon_oid[.] %>% as.character

df <- df %>% subset(Strain2 == "CL28") %>% droplevels

#Load previous data from the tripartite experiments
Dat <- readRDS(file = "../cleandata/dataset_tripartite_moduleAvs4RGIs.RDS")
Dat <- Dat %>% subset(InhibitedStrain == "CL28") %>% droplevels
meta_original <- read.table("../rawdata/metadata_124_isolates_burkholderia.tsv",header = T,
                   sep = "\t",quote = "",comment.char = "")
meta_original <- meta_original$Genome.Name %>% grep(pattern = "MF|CL14|110B")  %>% meta_original[.,]  %>%
  droplevels
meta_original <- meta_original$Genome.Name %>% grep(pattern = "Variovorax") %>% meta_original[.,] %>% droplevels
meta_original <- meta_original$Genome.Name %>% grep(pattern = "CL14|160",invert = T) %>% meta_original[.,] %>%
  droplevels
meta_original$InhibitorStrain <- c("MF4","MF295","MF278","MF350","MF369","MF110","MF375","MF349")
meta_original <- meta_original %>% subset(InhibitorStrain != "MF295")
Dat_sub <- which(Dat$InhibitorStrain  %in% c("MF4","MF278","MF350","MF369","MF110","MF375","MF349")) %>%
  Dat[.,] %>% droplevels
Dat_sub <- merge(Dat_sub,meta_original, by = "InhibitorStrain")

df_tri <- data.frame(Length = Dat_sub$MainRoot,Strain1 = Dat_sub$InhibitorStrain,Strain2 = Dat_sub$InhibitedStrain,
           Exp = rep("Tripartite",nrow(Dat_sub)),taxon_oid = Dat_sub$taxon_oid, Genus = Dat_sub$Genus)
df$taxon_oid <- df$taxon_oid %>% as.character
df_all <- rbind(df,df_tri)

#Prepate the df stat
df_stat <- df_stat %>% subset(Strain2 == "CL28") %>% droplevels 
head(df_stat)
temp <- Dat_sub[,c(1,6,8,11)] %>% unique
df_stat_tri <- data.frame(Strain1 = temp$InhibitorStrain,Strain2 = temp$InhibitedStrain,pvalue = temp$p.adj,padj = temp$p.adj,
           Label = rep("***",nrow(temp)),taxon_oid = temp$taxon_oid)
df_stat_all  <- rbind(df_stat,df_stat_tri)

#Refactor
df_all$taxon_oid <- df_all$taxon_oid %>% factor(levels = mtips)
df_stat_all$taxon_oid <- df_stat_all$taxon_oid %>% factor(levels = mtips)

p <- ggplot(data = df_all,aes(taxon_oid,Length)) +
  facet_grid(.~Strain2) + 
  geom_boxplot(aes(color = Genus),outlier.size = NA,outlier.shape = NA) +
  geom_sina(size = 2,aes(color = Genus),alpha = 0.2)+
  geom_hline(yintercept = 3,size = 1,colour="black", linetype="dashed") +
  geom_text(data = df_stat_all  %>% subset(Strain2 == "CL28"),aes(x = taxon_oid,y = 8,label = Label),
            inherit.aes = F,size = 5,family ="Arial",angle = -90) +
geom_rect(data = df_segment %>% subset(Strain2 == "CL28"),
            mapping = aes(ymin = two,ymax = four,xmin = 0,xmax = Inf),
            inherit.aes = F,color = NA,fill = "#56B4E9",alpha = 0.2) +
  theme_ohchibi() + coord_flip() +
  theme(
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90,vjust = 0.5),
        legend.position = "none"
        )

oh.save.pdf(p = p,outname = "primaryroot_elongation_all_public_variovorax_boxplotsvscl28.pdf",
            outdir = "../figures/",height = 15,width = 10)


