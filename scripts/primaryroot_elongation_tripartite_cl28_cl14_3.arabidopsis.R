library(ohchibi)
library(paletteer)
library(palettesPM)
library(emmeans)


setwd('/home/isai/Documents/results/github/variovoraxRGI/scripts/')
set.seed(130816)

dir.create("../cleandata/")
dir.create("../figures//")

#Read the datasets of monoassociation
dat1 <- readRDS("../cleandata/dataset_tripartite_moduleAvs4RGIs.RDS")
long1 <- readRDS(file = "../cleandata/dataset_inhibitors_moduleAvs4RGIs.RDS")
dat2 <- readRDS("../cleandata/dataset_tripartite_variovoraxburkholderiavsRGIs.RDS")
long2 <- readRDS(file = "../cleandata/dataset_inhibitors_variovoraxburkholderiavsRGIs.RDS")

mono1 <- dat1 %>% subset(InhibitedStrain == "CL28" & 
                           InhibitorStrain == "NB")

long1 <- long1 %>% subset(Treatment == "CL14") %>% droplevels
duo1 <- dat1 %>% subset(InhibitedStrain == "CL28" & 
                          InhibitorStrain == "CL14")

df1 <- rbind(mono1,duo1)
df1$Rep <- rep("Rep1",nrow(df1))

mono2 <- dat2 %>% subset(InhibitedStrain == "CL28" & 
                           InhibitorStrain == "NB")
long2cl14 <- long2 %>%  subset(InhibitedStrain == "NB" & 
                                 InhibitorStrain == "CL14")
long2nb <- long2 %>%  subset(InhibitedStrain == "NB" & 
                               InhibitorStrain == "NB")

duo2 <- dat2 %>% subset(InhibitedStrain == "CL28" & 
                          InhibitorStrain == "CL14")
df2 <- rbind(mono2,duo2)
df2$Rep <- rep("Re2",nrow(df2))

#Create a united structure
df1 <- df1[,c(2,7,11)]
cl14_1 <- data.frame(MainRoot = long1$MainRoot,
                     Comparison = rep("CL14-NB",nrow(long1)),
                     Rep = rep("Rep1",nrow(long1))
)
df1 <- rbind(df1,cl14_1)


df2 <- df2[,c(3,5,4)]
cl14_2 <- data.frame(MainRoot = long2cl14$MainRoot,
                     Comparison = rep("CL14-NB",nrow(long2cl14)),
                     Rep = rep("Rep2",nrow(long2cl14))
)
nb_2 <- data.frame(MainRoot = long2nb$MainRoot,
                   Comparison = rep("NB-NB",nrow(long2nb)),
                   Rep = rep("Rep2",nrow(long2nb))
)
df2 <- rbind(df2,cl14_2,nb_2)

df <- rbind(df1,df2)
df$Comparison <- df$Comparison %>%
  factor(levels = c("NB-NB","CL14-NB","CL28-NB","CL28-CL14"))

#Compare
res_em <- aov(data = df,formula = MainRoot ~ Comparison) %>%
  emmeans(specs = "Comparison") %>% CLD

nb_quantiles <- df %>% subset(Comparison == "NB-NB") %$% MainRoot %>% quantile
df_nb <- data.frame(lower = nb_quantiles[2],upper = nb_quantiles[4])
full_quantiles <- df %>% subset(Comparison == "CL14-NB") %$% MainRoot %>% quantile
df_full <- data.frame(lower = full_quantiles[2],upper = full_quantiles[4])


p <- chibi.boxplot(Map = df,x_val = "Comparison",y_val = "MainRoot",
                   style = "open",
                   size_median = 8,
                   size_point = 12,size_axis_text.y = 35,
                   size_axis_text.x = 30,median_color = "red",
                   size_axis_title.y = 40,strip_text_size = 30,
                   size_legend_text = 30,legend_proportion_size = 4) +
  geom_vline(xintercept = c(1.5,2.5,3.5,4.5,5.5),size = 0.6, color = "#D9D9D9")  +
  theme(axis.title.x = element_blank()) + ylab(label = "Main root elongation (cm)") +
  geom_rect(data = df_nb,
            mapping = aes(ymin = lower,ymax = upper,xmin = 0,xmax = Inf),
            inherit.aes = F,alpha = 0.2,color = NA,fill = "#56B4E9") +
  geom_rect(data = df_full,
            mapping = aes(ymin = lower,ymax = upper,xmin = 0,xmax = Inf),
            inherit.aes = F,alpha = 0.1,color = NA,fill = "#5E5E5E")

saveRDS(object = df,file = "../cleandata/dataset_tripartite_cl28_cl14.RDS")
oh.save.pdf(p = p,outname = "primaryroot_elongation_tripartite_cl28_cl14.pdf",
            outdir = "../figures",height = 20,width = 18)
