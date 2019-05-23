library(ohchibi)
setwd('/home/isai/Documents/results/github/variovoraxRGI/scripts/')
set.seed(130816)


#Read the clean matrix and create dataset 
Dat_amplicon<-readRDS(file = "../cleandata/dat_amplicon_4stresses_useq97.RDS")
Dat <- Dat_amplicon$RelativeAbundance

Tax <- Dat$Tax
vario_useqs <- Tax$Taxonomy %>% grep(pattern = "Variovorax") %>%
  Tax$ID[.] %>% as.character

melted <- Dat$Tab %>% as.matrix %>% melt
colnames(melted) <- c("USeq","ID_Matrix","value")

melted <- merge(Dat$Map,melted, by = "ID_Matrix")

melted_vario <- melted$USeq %in% vario_useqs %>% which %>% melted[.,] %>% droplevels
melted_vario$Genus <- rep("Variovorax",nrow(melted_vario))

melted_vario_sum <- dcast(data = melted_vario,formula = ID_Matrix ~ Genus,fun.aggregate = mean,value.var = "value")
melted_vario_sum <- merge(Dat$Map,melted_vario_sum,by = "ID_Matrix")




##Load the Dropout data
Dat_drop <- readRDS(file = "../cleandata/dat_amplicon_variovoraxdropout_useq97.RDS")
Dat_drop <- Dat_drop$RelativeAbundance

melted_drop <- Dat_drop$Tab %>% as.matrix %>% melt
colnames(melted_drop) <- c("USeq","ID_Matrix","value")

melted_drop <- merge(Dat_drop$Map,melted_drop, by = "ID_Matrix")

melted_drop_vario <- melted_drop$USeq %in% vario_useqs %>% which %>% melted_drop[.,] %>% droplevels
melted_drop_vario$Genus <- rep("Variovorax",nrow(melted_drop_vario))

melted_drop_vario_sum <- dcast(data = melted_drop_vario,formula = ID_Matrix ~ Genus,fun.aggregate = mean,value.var = "value")
melted_drop_vario_sum <- merge(Dat_drop$Map,melted_drop_vario_sum,by = "ID_Matrix")

#Subset
melted_drop_vario_sum <- melted_drop_vario_sum %>% subset( 
  syncom == "I_Full" | syncom ==  "Full") %>% droplevels


#Merge both data frames
a <- melted_vario_sum[,c(5,9,13)]
b <- data.frame(
  experiment = rep("Dropout",nrow(melted_drop_vario_sum)),
  typebyTissue = melted_drop_vario_sum$fraction,
  Variovorax = melted_drop_vario_sum$Variovorax
)

all <- rbind(a,b)

all <- all %>% subset(typebyTissue != "Agar_NoPlant") %>% droplevels
all$experiment <- all$experiment  %>% gsub(pattern = "Pi",replacement = "Phosphate") %>%
  gsub(pattern = "salinity",replacement = "Salinity") %>%
  gsub(pattern = "T",replacement = "Temperature") %>%
  factor(levels = c("Salinity","Phosphate","pH","Temperature","Dropout"))


all$typebyTissue <- all$typebyTissue %>% as.character %>%
  gsub(pattern = "Agar_Plant",replacement = "AgarPlant") %>%
  factor(levels = c("Inoculum","AgarNoPlant","AgarPlant","Root","Shoot"))

all <- all %>% subset(typebyTissue != "AgarNoPlant") %>% droplevels

#Define paleta
pal_frac <- palettesPM::pm.colors.fractions()
pal_frac <- match(all$typebyTissue %>% levels,names(pal_frac)) %>%
  pal_frac[.]

p <- chibi.boxplot(Map = all,x_val = "typebyTissue",y_val = "Variovorax",col_val = "typebyTissue",
                   facet_formula = "experiment",mpalette = pal_frac,size_boxplot = 1,
                   style = "open",size_point = 8,size_axis_text.y = 30,size_axis_text.x = 0,alpha_point = 0.8,
                   size_axis_title.x = 40,size_axis_title.y = 40,size_median = 6,median_color = "black",legend_proportion_size = 4,size_legend_text = 30,strip_text_size = 30) + 
  ylab(label = "Variovorax Relative Abundance")  + xlab(label = "Fraction") + 
  geom_vline(xintercept = c(1.5,2.5,3.5,4.5),size = 0.3, color = "#D9D9D9") +
  theme(axis.title.x = element_blank())


p <- chibi.boxplot(Map = all,x_val = "typebyTissue",y_val = "Variovorax",col_val = "typebyTissue",mpalette = pal_frac,size_boxplot = 1,
              style = "open",size_point = 8,size_axis_text.y = 30,size_axis_text.x = 0,alpha_point = 0.8,
              size_axis_title.x = 40,size_axis_title.y = 40,size_median = 6,median_color = "black",
              legend_proportion_size = 4,size_legend_text = 30,strip_text_size = 30) +
  theme(panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank(),axis.title.x = element_blank()) +
  ylab("Relative Abundance")


p <- chibi.boxplot(Map = all,x_val = "typebyTissue",y_val = "Variovorax",size_boxplot = 1,
                   style = "open",size_point = 8,size_axis_text.y = 30,size_axis_text.x = 0,alpha_point = 0.8,
                   size_axis_title.x = 40,size_axis_title.y = 40,size_median = 6,median_color = "red",
                   legend_proportion_size = 4,size_legend_text = 30,strip_text_size = 30) +
  theme(panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank(),axis.title.x = element_blank()) +
  ylab("Relative Abundance")

p <- chibi.boxplot(Map = all,x_val = "typebyTissue",y_val = "Variovorax",size_boxplot = 1,
              style = "open",size_point = 0,size_axis_text.y = 30,size_axis_text.x = 0,alpha_point = 0,
              size_axis_title.x = 40,size_axis_title.y = 40,size_median = 6,median_color = "red",
              legend_proportion_size = 4,size_legend_text = 30,strip_text_size = 30) +
  theme(panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank(),axis.title.x = element_blank()) +
  ylab("Relative Abundance") + geom_sina(alpha = 0.3, 
                                         stroke = 0,size = 3,color = "black") + scale_y_continuous(limits = c(0,0.02))


oh.save.pdf(p = p,outname = "amplicon_4stresses_relative_abundance_variovorax.pdf",outdir = "../figures/",width = 18,height = 9)
