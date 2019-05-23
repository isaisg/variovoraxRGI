library(ohchibi)
library(emmeans)

setwd('/home/isai/Documents/results/github/variovoraxRGI/scripts/')
set.seed(130816)

Map <- read.table(file = "../rawdata/dataS10_bacterialabsoluteabundance_cl28variovorax.tsv",header = T,
                  sep = "\t",quote = "",comment.char = "")
Map$Partner <- Map$Partner %>%
  factor(levels = c("No bacteria","Variovorax CL14","Variovorax B4"))

p <- chibi.boxplot(Map = Map,x_val = "Fraction",y_val = "LCountperg",col_val = "Fraction",
                   facet_formula = "Partner",style = "mix",mpalette = c("saddlebrown"),
                   median_color = "#414141",size_point = 15,size_median = 6,alpha_point = 0.8,
                   size_axis_text.x = 0,size_boxplot = 1,size_axis_text.y = 30,
                   size_axis_title.x = 0,size_axis_title.y = 35,
                   size_legend_text = 35,legend_proportion_size = 4,strip_text_size = 35) 
oh.save.pdf(p = p,outname = "bacterialabsoluteabundance_tripartite_cl28variovorax.pdf",outdir = "../figures/",height = 20,width = 18)

aov(formula = LCountperg ~ Partner,data = Map) %>% 
  emmeans::emmeans(specs = "Partner") %>% emmeans::CLD()
