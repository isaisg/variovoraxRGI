library(ohchibi)
library(emmeans)
library(RColorBrewer)
library(egg)
library(palettesPM)

setwd('/home/isai/Documents/results/github/variovoraxRGI/scripts/')
set.seed(130816)

dir.create("../cleandata/")
dir.create("../figures/")

Map <- read.table(file = "../rawdata/dataS4_primaryroot_elongation_monoassociation.tsv",header = T,
           sep = "\t",quote = "",comment.char = "")

#Calculate normalization factor using NB
nb_me <- subset(Map,Strain == "NB") %$% Length %>% mean
Map_nb <- subset(Map,Strain == "NB")  %>% droplevels
Map_norm <- NULL
for(exp in levels(Map_nb$Experiment)){
  me_exp  <- subset(Map_nb, Experiment == exp) %$% Length %>% mean
  nf <- nb_me /me_exp
  map_temp <- subset(Map, Experiment == exp) %>% droplevels
  map_temp$NormLength <- map_temp$Length*nf
  Map_norm <- rbind(Map_norm,map_temp)
}




#Aggregate samples
df_exp <- aggregate(Experiment ~ Strain,Map_norm,table)
rStrains <- df_exp$Strain
df_exp <- df_exp[,-1]
rownames(df_exp) <- rStrains
df_exp[df_exp >1 ] <- 1
df_exp <- rowSums(df_exp) %>% data.frame
all_Strains <- rownames(df_exp)[which(df_exp>1)] %>% grep(pattern = "NB",value = T,invert = T)

#Anova framework to compute significance
Res <- NULL
for(st in all_Strains){
  print(st)
  Map_st <- subset(Map_norm,Strain == st) %>% droplevels
  Map_nb <- subset(Map_norm, Strain == "NB")  %>% droplevels
  Map_nb <- Map_nb[Map_nb$Experiment %in% levels(Map_st$Experiment) %>% which,] %>% droplevels
  map <- rbind(Map_nb,Map_st)
  map$Strain <- relevel(x = map$Strain,ref = "NB")
  rdf <- aov(formula = NormLength ~ Strain + Experiment,data = map)  %>% emmeans(specs = "Strain") %>%
    pairs(reverse = TRUE,adjust = "none") %>% data.frame 
  rdf$contrast <- rdf %$% contrast %>% as.character %>% 
    gsub(pattern = " ",replacement = "") %>% factor
  rdf$Strains <- st
  Res <- rbind(Res,rdf)
}

#Correct for multiple testing
Res$q.value <- Res$p.value %>% p.adjust(method = "fdr")
Res$l10q.value <- -log10(Res$q.value)
Res$Significance <- rep("NoSignificant",nrow(Res))
Res$Significance[which(Res$q.value < 0.05)] <- "Significant"
Res$Significance <- Res$Significance %>% factor(levels = c("NoSignificant","Significant"))
paleta_sig <- c("#1A3431","#FF854E")

#Plot FDR pvalues versus coefficient of the anova
p <- ggplot(data = Res,aes(l10q.value,estimate,fill = Significance)) +
  geom_point(size=10,alpha=1,shape = 21) + scale_fill_manual(values = paleta_sig)+ 
  geom_hline(yintercept = c(-2),linetype="dashed",lwd=1.2,color = "#414141") + 
  geom_vline(xintercept = c(1),linetype="dashed",lwd=1.2 , color = "#414141") +
  xlab(label = "-log10 q.value ") + ylab(label = "Coefficient") + 
  theme(axis.line = element_blank(),
        panel.background = element_rect(fill = 'white'),
        panel.grid.major = element_line(colour =   "#D9D9D9"),
        panel.grid.minor = element_line(colour = "#D9D9D9"),
        panel.border = element_rect(fill=NA,color =  "#414141",size = 1),
        axis.ticks = element_line(colour = "black",size = 2.5),
        axis.text.x = element_text(family = "AvantGarde",face = "plain",size =30,colour="#414141"),
        axis.text.y = element_text(family = "AvantGarde",face="plain",size=30,colour="#414141"),
        axis.title.x = element_text(family = "AvantGarde",face="plain",size = 40,colour = "#414141"),
        axis.title.y = element_text(family = "AvantGarde",face="plain",size=40,colour="#414141"),
        legend.background = element_blank(),legend.key.size = unit(40,"point"),
        legend.title=element_blank(),legend.key = element_blank(), 
        legend.text = element_text(size=25,
                                   family = "AvantGarde",face = "plain",colour = "#414141"),
        legend.position ="right",strip.background = element_blank(),
        strip.text.x = element_text(family = "AvantGarde",face = "plain",size =30,colour="#414141")
  ) 


#Map norm
Map_norm$Strain <- Map_norm$Strain %>% gsub(pattern = " ",replacement = "")
Map_norm$Strain <- Map_norm$Strain %>% gsub(pattern = "a",replacement = "A")
Map_norm$Strain <- Map_norm$Strain %>% factor

meta <- read.table("../rawdata/metadata_97useq_to_ids.tsv",
                   header = T,sep = "\t",check.names = F,comment.char = "",quote = "")


meta$Freezer_Id <- meta$Freezer_Id %>%gsub(pattern = "_.*",replacement = "") %>% as.character
meta_2cluster <- read.table(file = "../rawdata/amplicon_4stresses_map_useq2cluster_abcd.tsv",header = T,sep = "\t")

meta <- merge(meta,meta_2cluster, by = "Useq")
meta <- meta[,c(1,2,6,7)]
colnames(meta)[2] <- "Strain"
Map_norm <- merge(Map_norm, meta, by = "Strain", all = TRUE)
Map_norm$Cluster <- Map_norm$Cluster %>% as.character

Map_norm$Cluster[Map_norm$Strain %>% grep(pattern = "NB")] <- "NB"
Map_norm <- Map_norm[-which(is.na(Map_norm$Cluster)),] %>% droplevels
Map_norm$Cluster <- Map_norm$Cluster %>% factor(levels = c("NB","A","B","C","D"))

#Order them 
df_mean <- aggregate(NormLength ~ Strain,data = Map_norm,FUN = "mean")
df_mean <- with(df_mean,order(-NormLength)) %>% df_mean[.,]
Map_norm$Strain <- Map_norm$Strain %>% factor(levels = df_mean$Strain)


#Create df to mape geom_rect
#Here use the combinatorial shade
df_combi <- read.table(file = "../rawdata/dataS3_primaryroot_elongation_modules.tsv",header = T,sep = "\t")
df_combi <- which(df_combi$SynCom %in% c("NB","A","B","C","D")) %>% df_combi[.,] %>% droplevels 
df_seg_module <- aggregate(Length ~ SynCom,df_combi,FUN = quantile) 
df_seg_module_end <- df_seg_module$Length %>% as.data.frame
df_seg_module_end$SynCom <- df_seg_module$SynCom
colnames(df_seg_module_end) <- c("one","two","three","four","five","Cluster")
df_seg <- df_seg_module_end %>% subset(Cluster == "NB") %>% droplevels
dg_seg <- df_seg[,1:5]
df_seg <- rbind(df_seg,df_seg,df_seg,df_seg)
df_seg$Cluster <- c("A","B","C","D")

##Append taxonomy
mphyla <- NULL
mclass <- NULL
for(ele in Map_norm$Genome_Classification){
  if(is.na(ele)){
    mphyla <- c(mphyla,NA)
    mclass <- c(mclass,NA)
  }else{
    temp_phyla <- ele %>% as.character %>% strsplit(split = "\\;") %>% unlist %>%
      grep(pattern = "p__",value = T) %>% gsub(pattern = "p__",replacement = "") %>%
      as.character
    mphyla <- c(mphyla,temp_phyla)
    temp_class <- ele %>% as.character %>% strsplit(split = "\\;") %>% unlist %>%
      grep(pattern = "c__",value = T) %>% gsub(pattern = "c__",replacement = "") %>%
      as.character
    mclass <- c(mclass,temp_class)
  }
  
}
Map_norm$Phylum <- mphyla
Map_norm$Class <- mclass
Map_norm$Color <- Map_norm$Phylum %>% as.character
Map_norm$Color[which(Map_norm$Class == "Alphaproteobacteria")] <- "Alphaproteobacteria"
Map_norm$Color[which(Map_norm$Color == "Proteobacteria")] <- "Gammaproteobacteria"
Map_norm$Color <- Map_norm$Color %>% factor
paleta <- palettesPM::pm.colors.phyla()
paleta <- which(names(paleta) %in% c("Actinobacteria","Bacteroidetes","Firmicutes")) %>%
  paleta[.]
paleta_proteo <- c("#CC873A","#FF6100")
names(paleta_proteo) <- c("Alphaproteobacteria","Gammaproteobacteria")
paleta <- c(paleta,paleta_proteo)

p_mono <- Map_norm %>% subset(Strain != "NB") %>% droplevels %>%
  ggplot(data = .,aes(Strain,NormLength,fill = Color)) +
  geom_point(size = 0,alpha = 0) + 
  geom_rect(data = df_seg,
            mapping = aes(ymin = two,ymax = four,xmin = 0,xmax = Inf),
            inherit.aes = F,color = NA,fill = "#56B4E9",alpha = 0.2) +
  geom_rect(data = df_seg_module_end %>% subset(Cluster != "NB") %>% droplevels,
            aes(ymin = two,ymax = four,xmin = 0,xmax = Inf),inherit.aes = F,color =NA,fill = "#00B233",alpha = 0.2) +
  geom_hline(yintercept = 3,size = 1,colour="black", linetype="dashed") +
  geom_boxplot(outlier.size = 0,outlier.stroke = 0,alpha = 1) + 
  facet_grid(~Cluster,scales = "free",space = "free") + theme_ohchibi(font_family = "Arial") +
  ylab(label = "Main Root Elongation (cms)") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none") +
  scale_fill_manual(values = paleta,na.value = "#D9D9D9")+
  scale_y_continuous(limits = c(0,8)) +
  coord_cartesian(expand = TRUE)


oh.save.pdf(p = p_mono,outname = "primaryroot_elongation_monoassociation_boxplots.pdf",outdir = "../figures/")
saveRDS(object = p_mono,file = "../figures/primaryroot_elongation_monoassociation_plots.RDS")

#Wrote clean dataset
write.table(x = Map_norm,file = "../cleandata/primaryroot_elongation_monoassociation_normalized.tsv",
            append = F,quote = F,sep = "\t",row.names = F,col.names = T)



