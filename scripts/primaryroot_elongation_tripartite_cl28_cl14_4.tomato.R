library(ohchibi)
library(emmeans)
setwd('/home/isai/Documents/results/github/variovoraxRGI/scripts/')

set.seed(130816)

###Tomato 
#Read tomato data
df_tomato <- read.table(file = "../rawdata/dataS9_primaryroot_elongation_tripartite_cl28cl14tomato.tsv",
                        header = T,sep = "\t")

df_tomato$SynCom <- df_tomato$SynCom %>% factor(levels = c("NB","CL28","CL28+CL14"))

m1 <- aov(Length ~ SynCom + Rep ,data = df_tomato)
Res_em_tomato <- m1 %>%
  emmeans(object = .,specs = "SynCom") %>% CLD %>%
  as.data.frame
Res_em_tomato$.group <- Res_em_tomato$.group %>% gsub(pattern = " ",replacement = "")
colnames(Res_em_tomato)[which(colnames(Res_em_tomato) == ".group")] <- "group"


#Adjust the letters of tukey
Res_em_tomato$LLetters  <- Res_em_tomato$group %>% 
  strsplit(split = "") %>% 
  lapply(X = .,FUN = function(x)letters[as.numeric(x)])

#Loop 
Res_em_tomato$Letters <- rep(0,nrow(Res_em_tomato))
for(i in 1:nrow(Res_em_tomato)){
  Res_em_tomato$Letters[i] <- Res_em_tomato$LLetters[[i]] %>% 
    unlist %>% toupper  %>% paste0(collapse = "")
}
Res_em_tomato$condition <- rep("Tomato",nrow(Res_em_tomato))



p <- chibi.boxplot(Map = df_tomato,x_val = "SynCom",
                   y_val = "Length",
                   #facet_formula = "Rep",
                   style = "open",
                   size_median = 8,
                   size_point = 12,size_axis_text.y = 35,
                   size_axis_text.x = 30,median_color = "red",
                   size_axis_title.y = 40,strip_text_size = 30,
                   size_legend_text = 30,legend_proportion_size = 4) +
  theme(
    axis.title.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank()
  ) + 
  ylab(label = "Primary root elongation (cm)")  +
  geom_text(data = Res_em_tomato,aes(y =12, label = Letters),size = 10)

oh.save.pdf(p = p,outname = "primatyroot_elongation_tripartite_cl28_cl14_tomato.pdf",
            outdir = "../figures",height = 20,width = 18)
