library(ohchibi)
library(emmeans)
library(multcompView)
library(tidyr)
library(paletteer)

set.seed(130816)

setwd('/home/isai/Documents/results/github/variovoraxRGI/scripts/')

Map <- read.table(file = "../rawdata/dataS13_primaryroot_elongation_hormone_MAMP.tsv",
                  header = T,sep = "\t",quote = "",comment.char = "")

Map$Strain <- Map$Strain  %>%  
  factor(levels = c("NB","CL11","CL14","MF160","B4","YR216"))


Res_em <- NULL
for (treat in levels(Map$Treatment)){
  df <- Map %>% subset(Treatment == treat) %>% droplevels
  m1 <- aov(formula = Length ~ Strain ,data = df)
  m1_em <- emmeans(m1, specs = "Strain") %>% CLD() 
  m1_em$.group <- m1_em$.group %>% gsub(pattern = " ",replacement = "")
  colnames(m1_em)[which(colnames(m1_em) == ".group")] <- "group"
  m1_em$Treatment <- rep(treat,nrow(m1_em))
  Res_em <- rbind(Res_em,m1_em)
}


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

#Subset the treatments for main figure
g1 <- c("CL28","IAA_100nM","2_4-D_100nM","ACC_100nM","flg22_1000nM","BAP_100nM","Zeatin_100nM")
df_sel <- which(Map$Treatment %in% g1) %>% Map[.,] %>% droplevels
df_em <- which(Res_em$Treatment %in% g1) %>% Res_em[.,] %>%
  droplevels
df_sel$Treatment <- df_sel$Treatment %>% 
  factor(levels = c("CL28","IAA_100nM","2_4-D_100nM","ACC_100nM","flg22_1000nM","BAP_100nM","Zeatin_100nM"))
df_em$Treatment <- df_em$Treatment %>% 
  factor(levels = c("CL28","IAA_100nM","2_4-D_100nM","ACC_100nM","flg22_1000nM","BAP_100nM","Zeatin_100nM"))

p <- chibi.sina(Map = df_sel,x_val = "Strain",y_val = "Length",
                facet_formula = "Treatment",size_point = 11,size_bar = 3,alpha_point = 0.15,
                df_stats = df_em,mean_var = "emmean",ymin_var = "lower.CL",ymax_var = "upper.CL",
                size_axis_text.x = 0,show_points = T,font_family = "Arial",
                size_axis_text.y = 30,size_axis_title.x = 35,size_axis_title.y = 35,
                size_legend_text = 35,strip_text_size = 35,legend_proportion_size = 4) 

p <- p + geom_text(data = df_em,
                   aes(x = Strain,y = 7.8,label = Letters),
                   inherit.aes = F,size = 8,family ="Arial") +
  ylab(label = "Main root elongation (cm)") +
  theme(axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none",
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank()
        )

#Create checkerboard panel
#Create matrix to manuallay adjust
mat <- rep(c(c(1,0,0,0,0,0),c(0,1,0,0,0,0),c(0,0,1,0,0,0),c(0,0,0,1,0,0),c(0,0,0,0,1,0),c(0,0,0,0,0,1)),7) %>%
  matrix(data = .,nrow = 6,byrow = T)
rownames(mat) <- c("NB","CL11","CL14","MF160","B4","YR216")  %>% rev

colnames(mat) <- levels(df_sel$Treatment %>% droplevels) %>% sapply(FUN = function(x)rep(x,6)) %>%
  as.vector

mat <- mat %>% melt
mat$value <- factor(mat$value)

mat$Var1 <- factor(mat$Var1,levels = mat$Var1 %>% levels)
p_mat <- ggplot(data = mat,mapping = aes(Var1,rev(Var1), fill = value)) + 
  geom_tile()  + theme_ohchibi() +facet_grid(.~Var2,scales = "free",space = "free")+
  scale_fill_manual(values = c("#ffffff00","black")) +
  theme(axis.title.x = element_blank(),
        strip.text.x = element_blank(),
        axis.ticks.x = element_blank(),axis.text.x = element_blank(),
        panel.grid  =  element_blank(),legend.position = "none",
        panel.background = element_blank(),
        panel.grid.major.y = element_blank())  +
  ylab(label = "Strain") +
  geom_hline(yintercept = c(1.5,2.5,3.5,4.5,5.5),size = 0.6, color = "#D9D9D9") + 
  geom_vline(xintercept = c(1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5,10.5,11.5,12.5), size = 0.6, color = "#D9D9D9") 


composition <- egg::ggarrange(p,p_mat, ncol = 1, heights = c(1,0.2),debug = F)


oh.save.pdf(p = composition,outname = "primaryroot_elongation_hormone_MAMP.pdf",outdir = "../figures/",width = 30,height = 20)


