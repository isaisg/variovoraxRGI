library(ohchibi)
library(emmeans)
library(paletteer)
library(egg)

setwd('/home/isai/Documents/results/github/variovoraxRGI/scripts/')
set.seed(130816)
dir.create("../cleandata/")
dir.create("../figures/")

#Read dataset
df <- read.table(file = "../rawdata/dataS14_primaryroot_elongation_DR5_tripartitecl28cl14.tsv",
                 header = T,sep = "\t",comment.char = "") %>%
  droplevels
df$Strain <- df$Strain %>% factor(levels = c("NB","CL28","CL28+CL14"))

### Run anova per Timepoint
Res_em <- NULL
Res_tukey <- NULL
for(tp in levels(df$Time)){
  m1 <- df %>% subset(Time == tp) %>% droplevels %>%
    aov(GFP_Intensity ~ Strain ,data = .)
  df_tukey <- TukeyHSD(m1) %$% Strain %>% as.data.frame
  df_tukey$Time <- rep(tp,nrow(df_tukey))
  df_tukey$Contrast <- rownames(df_tukey)
  rownames(df_tukey) <- NULL
  Res_tukey <- rbind(Res_tukey,df_tukey)
  temp <- m1 %>%  
    emmeans(object = .,specs = "Strain")  %>% CLD
  temp$Time <- rep(tp,nrow(temp))
  Res_em <- rbind(Res_em ,temp)
}

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


##Compare pair of Strains
Res_tukey$FDR <- p.adjust(p = Res_tukey$`p adj` ,method = "fdr")
which(Res_tukey$Contrast %in% c("CL28+CL14-CL28")) %>%
  Res_tukey[.,]



p<- ggplot(data = Res_em,mapping = aes(Time,emmean,group = Strain)) +
  geom_line(aes(y = emmean,linetype = Strain),size =3) +
  geom_point(aes(y = emmean, 
                 group = Strain,
                 shape = Strain), 
             ,size = 12,stroke = 3) +
  theme_ohchibi(legend_proportion_size = 3) + 
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank()) +
  ylab(label = "GFP intensity")  + 
  scale_shape_manual(values = c(1,2,0))+
  scale_linetype_manual(values = c("solid","dashed","dotdash")) +
  theme(
    panel.background = element_rect(fill = "#fff2ff"),
    legend.position = "none"
    
  )

oh.save.pdf(p = p,outname = "primaryroot_elongation_DR5_tripartitecl28cl14.pdf",outdir = "../figures/")
