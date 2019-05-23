library(ohchibi)
library(egg)
setwd('/home/isai/Documents/results/github/variovoraxRGI/scripts/')
set.seed(130816)

df <- read.table(file = "../rawdata/data_auxindegradation_variovoraxcl14_growth.tsv",header = T,sep = "\t")
df$Strain <- df$Strain %>% factor(levels = c("none","CL14"))

melted <- melt(data = df,measure.vars = c("OD600","Auxin"))

paleta <- c("#56B4E9","darkgreen")
p_od <- melted %>% subset(variable == "OD600") %>%
  ggplot(data =.,aes(Time,value)) +
  geom_point(aes(color = Strain),size = 7) +
  stat_summary(mapping = aes(group = Strain,color = Strain),
               fun.y = "mean",geom = "line",size = 4)+
  scale_color_manual(values = paleta) +
  theme_ohchibi() +
  facet_grid(.~Media,space = "free",scales = "free") +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),axis.ticks.x = element_blank()
  ) + ylab("OD")

p_auxin <- melted %>% subset(variable == "Auxin") %>%
  ggplot(data =.,aes(Time,value)) +
  geom_point(aes(color = Strain),size = 7) +
  stat_summary(mapping = aes(group = Strain,color = Strain),
               fun.y = "mean",geom = "line",size = 4)+
  scale_color_manual(values = paleta) +
  theme_ohchibi() +
  facet_grid(.~Media,space = "free",scales = "free") +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.title.x = element_blank(),
    strip.text.x = element_blank()
  )+ ylab("Auxin")


composition <- egg::ggarrange(p_od,p_auxin,ncol =1)
oh.save.pdf(p = composition,outname = "auxindegradation_variovoraxcl14growth.pdf",outdir = "../figures/",width = 20,height = 15)
