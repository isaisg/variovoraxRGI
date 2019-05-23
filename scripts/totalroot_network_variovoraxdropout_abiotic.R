library(ohchibi)
library(emmeans)

#The raw data were taken from supplementary datasets 2 and 3 from Nature 2012 Paper
setwd('/home/isai/Documents/results/github/variovoraxRGI/scripts/')
set.seed(130816)

dir.create("../cleandata")
dir.create("../figures")

#Read data
df <- read.table("../rawdata/dataS7_primaryroot_elongation_totalrootnetwork_variovoraxdropout.tsv",
                 header = T,sep = "\t",quote = "",comment.char = "")

df <- df %>% subset(Experiment == "TotalRootNetwork") %>% droplevels

#Reorder the stress levels
df$Condition <- df$Condition %>% 
  factor(c("Control","-P","Na","pH8.2","31deg","MS"))
df$SynCom <- df$SynCom %>% factor(levels = c("NB","Full","DOV"))


#TEst using ANOVaAframework
#Run anova indide each condition tested
Res_em <- NULL
for(ct in levels(df$Condition)){
  temp <- df %>% subset(Condition == ct) %>% droplevels %>%
    aov(Length ~ SynCom + Rep,data = .)  %>%
    emmeans(object = .,specs = "SynCom") %>% CLD %>%
    as.data.frame
  temp$Condition <- rep(ct,nrow(temp))
  Res_em <- rbind(Res_em,temp)
  
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

Res_em$Condition <- Res_em$Condition %>% factor(levels = levels(df$Condition))




palette_dropout_duo <- c("#56B4E9","#000000","#E69F00")
names(palette_dropout_duo) <- c("NB","Full","DOV")
p <- chibi.boxplot(Map = df,x_val = "SynCom",y_val = "Length",col_val = "SynCom",
                   facet_formula = "Condition",
                   style = "mix",mpalette = palette_dropout_duo,median_colored_as_points = T,
                   size_median = 4,
                   size_point = 0,stroke_point = 0,size_axis_text.y = 35,
                   size_axis_text.x = 30,median_color = "black",
                   size_axis_title.y = 40,strip_text_size = 30,
                   size_legend_text = 30,legend_proportion_size = 4) +
  scale_shape_manual(values = 21:22) + ylab(label = "Main root elongation (cm)") + 
  xlab(label = "Syncom") +
  geom_sina(alpha = 0.3, 
            stroke = 0,size = 3,aes(color = SynCom)) +
  scale_color_manual(values = palette_dropout_duo) +
  geom_text(data = Res_em,
            aes(x = SynCom,y = 9.8,label = Letters),
            inherit.aes = F,size = 5,family ="Arial")  +
  geom_vline(xintercept = c(1.5,2.5,3.5,4.5,5.5),size = 0.6, color = "#D9D9D9") +
  theme(legend.position = "none",
        axis.text.x = element_blank(),axis.title.x = element_blank())

oh.save.pdf(p = p,outname = "totalroot_network_variovoraxdropout_abiotic.pdf",outdir = "../figures/",
            width = 20,height = 10)

