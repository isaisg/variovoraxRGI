library(ohchibi)
library(emmeans)


setwd('/home/isai/Documents/results/github/variovoraxRGI/scripts/')
set.seed(130816)
dir.create("../figures/")
dir.create("../cleandata/")

#Read the rawdata
df <- read.table("../rawdata/dataS7_primaryroot_elongation_totalrootnetwork_variovoraxdropout.tsv",
                 header = T,sep = "\t",quote = "",comment.char = "")

df <- df %>% subset(Experiment == "Original") %>% droplevels
df$SynCom <- df$SynCom %>% factor(levels = c("NB","Full","DOB","DOV","DOVB"))
#Test using anova framewotk
m1_em <- df %>% aov(formula = Length ~ SynCom,data = .) %>% 
  emmeans(object = .,specs = "SynCom") %>% CLD
m1_em$.group <- m1_em$.group %>% gsub(pattern = " ",replacement = "")
colnames(m1_em)[which(colnames(m1_em) == ".group")] <- "group"
Res_em <- m1_em

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

#Create boxplot
nb_quantiles <- df %>% subset(SynCom == "NB") %$% Length %>% quantile
df_nb <- data.frame(lower = nb_quantiles[2],upper = nb_quantiles[4])
full_quantiles <- df %>% subset(SynCom == "Full") %$% Length %>% quantile
df_full <- data.frame(lower = full_quantiles[2],upper = full_quantiles[4])

palette_dropout_duo <- c("#56B4E9","#000000","#E69F00")
names(palette_dropout_duo) <- c("NB","Full","DOV")




#Alternative figure for main figure
palette_dropout_duo <- c("#56B4E9","#000000","#E69F00","#7F5800","#7F7F7F")
names(palette_dropout_duo) <- c("NB","Full","DOV","DOVB","DOB")
palette_dropout_duo <- palette_dropout_duo[c(1,2,5,3,4)]
p <- chibi.boxplot(Map = df,x_val = "SynCom",y_val = "Length",col_val = "SynCom",
                   style = "mix",mpalette = palette_dropout_duo,median_colored_as_points = T,
                   size_median = 4,
                   size_point = 0,stroke_point = 0,size_axis_text.y = 35,
                   size_axis_text.x = 30,median_color = "black",
                   size_axis_title.y = 40,strip_text_size = 30,
                   size_legend_text = 30,legend_proportion_size = 4) +
  scale_shape_manual(values = 21:22) + ylab(label = "Main root elongation (cm)") + 
  xlab(label = "SynCom") +
  geom_sina(alpha = 0.3, 
            stroke = 0,size = 8,aes(color = SynCom)) +
  scale_color_manual(values = palette_dropout_duo) +
  #geom_text(data = Res_em,
  #          aes(x = SynCom,y = 8,label = Letters),
  #          inherit.aes = F,size = 5,family ="Arial")  +
  #geom_vline(xintercept = c(1.5,2.5,3.5,4.5,5.5),size = 0.6, color = "#D9D9D9") +
  theme(legend.position = "none",
        axis.text.x = element_blank(),axis.title.x = element_blank()) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank()) +
  scale_y_continuous(limits = c(0,7.5))
#p <- p + geom_rect(data = df_nb,
#          mapping = aes(ymin = lower,ymax = upper,xmin = 0,xmax = Inf),
#          inherit.aes = F,alpha = 0.2,color = NA,fill = "#56B4E9") +
#  geom_rect(data = df_full,
#            mapping = aes(ymin = lower,ymax = upper,xmin = 0,xmax = Inf),
#            inherit.aes = F,alpha = 0.1,color = NA,fill = "#000000")
oh.save.pdf(p = p,outname = "primaryroot_elongation_variovoraxdropout_original.pdf",
            outdir = "../figures",height = 20,width = 18)

