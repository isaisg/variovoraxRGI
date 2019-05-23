library(ohchibi)
library(paletteer)
library(palettesPM)
library(egg)
library(RColorBrewer)
library(extrafont)
loadfonts(device = "pdf") 
library(emmeans)

### Optimized for the oh.ggsave function
size_legend_text <- 45
size_axis_title <- 20
size_axis_text <- 15
size_legend_title <- 55
legend_proportion_size <- 4
size_points <- 6
strip_text_size <- 25
size_median <- 3
alpha_point <- 0.7
stroke_point <- 0.5


setwd('/home/isai/Documents/results/github/variovoraxRGI/scripts/')

set.seed(seed = 130816)

#Read the dataset 
Dat_amplicon <- readRDS(file = "../cleandata/dat_amplicon_4stresses_useq97.RDS")
distfun <- function(x,method) vegan::vegdist(x = x, method = "bray")

#Define palettes
palette_variance <- paletteer_d(package = "dutchmasters",palette = "pearl_earring",11)[c(1:4,7)] %>%
  c(.,"white")
names(palette_variance) <- c("typebyTissue","Genotype","condition",
                             "typebyTissue:condition","typebyTissue:Genotype","Residual")
palette_variance <- palette_variance[c(1,3,4,6)]

colores_phosphate <- c("#ffffe5","#f7fcb9","#addd8e","#41ab5d","#238443","#005a32")
colores_salinity <- c("#fde0dd","#fa9fb5","#dd3497", "#7a0177")
colores_temp <- c("#FFC88C","#FFA33F","#7F6446")
colores_ph <- c("#C3C1EA","#7E78E5","#555466")



################# Alpha Diversity #############
Dat <- Dat_amplicon$Rarefied
Dat$Map$Shannon <- vegan::diversity(x = Dat$Tab, index = "shannon", MARGIN = 2 )

########### Salinity ##################
Dat_sub <- Dat %>% subset.Dataset(subset = experiment == "salinity"&
                                    typebyTissue != "Inoculum",drop = T,clean = T) %>%
  subset.Dataset(typebyTissue != "AgarNoPlant",drop = T,clean = T)
m1 <- Dat_sub %$% Map %>%
  aov(formula = Shannon ~ typebyTissue + condition + typebyTissue:condition + Rep,data = .)
Res_em <- m1 %>% emmeans(pairwise ~typebyTissue:condition) %$% emmeans %>% CLD %>% as.data.frame 
Res_em$.group <- Res_em$.group %>% gsub(pattern = " ",replacement = "")
colnames(Res_em)[which(colnames(Res_em) == ".group")] <- "group"

#Adjust the letters of tukey
Res_em$LLetters  <- Res_em$group %>%
  strsplit(split = "") %>%
  lapply(X = .,FUN = function(x)letters[as.numeric(x)])

#Loop
Res_em$Letters <- rep(0,nrow(Res_em))
for(i in 1:nrow(Res_em)){
  Res_em$Letters[i] <- Res_em$LLetters[[i]] %>%
    unlist %>% toupper  %>% paste0(collapse = "")
}

p <- chibi.boxplot(Map = Dat_sub$Map,facet_formula = "typebyTissue",x_val = "condition",
                   y_val = "Shannon",col_val = "condition",
                   size_median = size_median,style = "mix",
                   size_point = size_points,alpha_point = alpha_point,stroke_point = stroke_point,
                   median_colored_as_points = F,mpalette = colores_salinity,median_color = "black",
                   size_axis_text.x = size_axis_text,
                   size_axis_text.y = size_axis_text,
                   size_axis_title.x = size_axis_title,
                   size_axis_title.y = size_axis_title,
                   strip_text_size = strip_text_size,
                   legend_proportion_size = legend_proportion_size,font_family = "Arial") +
  ylab(label ="Shannon Diversity Index") + scale_y_continuous(limits = c(1.8,3.5)) 
  #geom_vline(xintercept = c(1.5,2.5,3.5),size = 0.3,color = "#D9D9D9")
p_alpha_salinity <- p + geom_text(data = Res_em,
                                  aes(x = condition,y = 3.5,label = Letters),
                                  inherit.aes = F,size = 4,family ="Arial",color = "black") +
  theme(legend.position = "none",axis.title.x = element_blank())

# p_anova <- chibi.anova(m1) + scale_fill_manual(values = palette_variance) + 
#   theme(legend.position = "none")

########### Phosphate ##################
Dat_sub <- Dat %>% subset.Dataset(subset = experiment == "Pi"&
                                    typebyTissue != "Inoculum",drop = T,clean = T) %>%
  subset.Dataset(typebyTissue != "AgarNoPlant",drop = T,clean = T)
m1 <- Dat_sub %$% Map %>%
  aov(formula = Shannon ~ typebyTissue + condition + typebyTissue:condition + Rep,data = .)
Res_em <- m1 %>% emmeans(pairwise ~typebyTissue:condition) %$% emmeans %>% CLD %>% as.data.frame 
Res_em$.group <- Res_em$.group %>% gsub(pattern = " ",replacement = "")
colnames(Res_em)[which(colnames(Res_em) == ".group")] <- "group"

#Adjust the letters of tukey
Res_em$LLetters  <- Res_em$group %>%
  strsplit(split = "") %>%
  lapply(X = .,FUN = function(x)letters[as.numeric(x)])

#Loop
Res_em$Letters <- rep(0,nrow(Res_em))
for(i in 1:nrow(Res_em)){
  Res_em$Letters[i] <- Res_em$LLetters[[i]] %>%
    unlist %>% toupper  %>% paste0(collapse = "")
}
p <- chibi.boxplot(Map = Dat_sub$Map,facet_formula = "typebyTissue",x_val = "condition",
                   y_val = "Shannon",col_val = "condition",
                   size_median = size_median,style = "mix",
                   size_point = size_points,alpha_point = alpha_point,stroke_point = stroke_point,
                   median_colored_as_points = F,mpalette = colores_phosphate,median_color = "black",
                   size_axis_text.x = size_axis_text,
                   size_axis_text.y = size_axis_text,
                   size_axis_title.x = size_axis_title,
                   size_axis_title.y = size_axis_title,
                   strip_text_size = strip_text_size,
                   legend_proportion_size = legend_proportion_size,font_family = "Arial") +
  ylab(label ="Shannon Diversity Index")  + scale_y_continuous(limits = c(1.8,3.5)) 
  #geom_vline(xintercept = c(1.5,2.5,3.5,4.5,5.5),size = 0.3,color = "#D9D9D9") 
p_alpha_phosphate <- p + geom_text(data = Res_em,
                                   aes(x = condition,y = 3.5,label = Letters),
                                   inherit.aes = F,size = 4,family ="Arial",color = "black") +
  theme(legend.position = "none",axis.title.x = element_blank())
# p_anova <- chibi.anova(m1) + scale_fill_manual(values = palette_variance) + 
#   theme(legend.position = "none")

########### pH ##################
Dat_sub <- Dat %>% subset.Dataset(subset = experiment == "pH"&
                                    typebyTissue != "Inoculum",drop = T,clean = T) %>%
  subset.Dataset(typebyTissue != "AgarNoPlant",drop = T,clean = T)
m1 <- Dat_sub %$% Map %>%
  aov(formula = Shannon ~ typebyTissue + condition + Rep,data = .)
Res_em <- m1 %>% emmeans(pairwise ~typebyTissue:condition) %$% emmeans %>% CLD %>% as.data.frame 
Res_em$.group <- Res_em$.group %>% gsub(pattern = " ",replacement = "")
colnames(Res_em)[which(colnames(Res_em) == ".group")] <- "group"

#Adjust the letters of tukey
Res_em$LLetters  <- Res_em$group %>%
  strsplit(split = "") %>%
  lapply(X = .,FUN = function(x)letters[as.numeric(x)])

#Loop
Res_em$Letters <- rep(0,nrow(Res_em))
for(i in 1:nrow(Res_em)){
  Res_em$Letters[i] <- Res_em$LLetters[[i]] %>%
    unlist %>% toupper  %>% paste0(collapse = "")
}
p <- chibi.boxplot(Map = Dat_sub$Map,facet_formula = "typebyTissue",x_val = "condition",
                   y_val = "Shannon",col_val = "condition",
                   size_median = size_median,style = "mix",
                   size_point = size_points,alpha_point = alpha_point,stroke_point = stroke_point,
                   median_colored_as_points = F,mpalette = colores_ph,median_color = "black",
                   size_axis_text.x = size_axis_text,
                   size_axis_text.y = size_axis_text,
                   size_axis_title.x = size_axis_title,
                   size_axis_title.y = size_axis_title,
                   strip_text_size = strip_text_size,
                   legend_proportion_size = legend_proportion_size,font_family = "Arial") +
  ylab(label ="Shannon Diversity Index") + scale_y_continuous(limits = c(1.8,3.5)) 
  #geom_vline(xintercept = c(1.5,2.5),size = 0.3,color = "#D9D9D9")

p_alpha_ph<- p + geom_text(data = Res_em,
                           aes(x = condition,y = 3.5,label = Letters),
                           inherit.aes = F,size = 4,family ="Arial",color = "black") +
  theme(legend.position = "none",axis.title.x = element_blank())

########### Temperature ##################
Dat_sub <- Dat %>% subset.Dataset(subset = experiment == "T"&
                                    typebyTissue != "Inoculum",drop = T,clean = T) %>%
  subset.Dataset(typebyTissue != "AgarNoPlant",drop = T,clean = T)
m1 <- Dat_sub %$% Map %>%
  aov(formula = Shannon ~ typebyTissue + condition + typebyTissue:condition + Rep,data = .)
Res_em <- m1 %>% emmeans(pairwise ~typebyTissue:condition) %$% emmeans %>% CLD %>% as.data.frame 
Res_em$.group <- Res_em$.group %>% gsub(pattern = " ",replacement = "")
colnames(Res_em)[which(colnames(Res_em) == ".group")] <- "group"

#Adjust the letters of tukey
Res_em$LLetters  <- Res_em$group %>%
  strsplit(split = "") %>%
  lapply(X = .,FUN = function(x)letters[as.numeric(x)])

#Loop
Res_em$Letters <- rep(0,nrow(Res_em))
for(i in 1:nrow(Res_em)){
  Res_em$Letters[i] <- Res_em$LLetters[[i]] %>%
    unlist %>% toupper  %>% paste0(collapse = "")
}
p <- chibi.boxplot(Map = Dat_sub$Map,facet_formula = "typebyTissue",x_val = "condition",
                   y_val = "Shannon",col_val = "condition",
                   size_median = size_median,style = "mix",
                   size_point = size_points,alpha_point = alpha_point,stroke_point = stroke_point,
                   median_colored_as_points = F,mpalette = colores_temp,median_color = "black",
                   size_axis_text.x = size_axis_text,
                   size_axis_text.y = size_axis_text,
                   size_axis_title.x = size_axis_title,
                   size_axis_title.y = size_axis_title,
                   strip_text_size = strip_text_size,
                   legend_proportion_size = legend_proportion_size,font_family = "Arial") +
  ylab(label ="Shannon Diversity Index") + scale_y_continuous(limits = c(1.8,3.5)) 
  #geom_vline(xintercept = c(1.5,2.5),size = 0.3,color = "#D9D9D9")

p_alpha_temp<- p + geom_text(data = Res_em,
                             aes(x = condition,y = 3.5,label = Letters),
                             inherit.aes = F,size = 4,family ="Arial",color = "black") +
  theme(legend.position = "none",axis.title.x = element_blank())

#Create composition 
p_alpha_phosphate <- p_alpha_phosphate + theme(panel.grid.minor.y = element_blank(),panel.grid.major.y = element_blank())
p_alpha_salinity <- p_alpha_salinity + theme(panel.grid.minor.y = element_blank(),panel.grid.major.y = element_blank())
p_alpha_ph <- p_alpha_ph + theme(panel.grid.minor.y = element_blank(),panel.grid.major.y = element_blank())
p_alpha_temp <- p_alpha_temp + theme(panel.grid.minor.y = element_blank(),panel.grid.major.y = element_blank())

composition <- egg::ggarrange(nrow = 2,ncol =2,
                              p_alpha_phosphate,p_alpha_salinity,p_alpha_ph,p_alpha_temp,
                              widths = c(1,1),heights = c(1,1),debug = F,byrow = TRUE
)

oh.save.pdf(p = composition,outname = "amplicon_4stresses_composition_alpha_diversity.pdf",outdir = "../figures/",
            height = 15,width = 20)

