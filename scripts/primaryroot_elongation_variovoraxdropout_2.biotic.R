library(ohchibi)
library(emmeans)
library(paletteer)

setwd('/home/isai/Documents/results/github/variovoraxRGI/scripts/')
set.seed(130816)


#Read raw data 
df <- read.table("../rawdata/dataS7_primaryroot_elongation_totalrootnetwork_variovoraxdropout.tsv",
                 header = T,sep = "\t",quote = "",comment.char = "")

df <- df %>% subset(Experiment == "Biotic") %>% droplevels
df$SynCom <- df$SynCom %>% factor(levels = c("Full","DOV"))

#Test each pairwise comparison using t.text
Res <- NULL
for(sync in levels(df$Condition)){
  temp <- df %>% subset(Condition == sync) %>% droplevels
  m1 <- temp %>% subset(SynCom == "DOV") %>% droplevels %$%
    Length %>% as.vector
  m2 <- temp %>% subset(SynCom == "Full") %>% droplevels %$%
    Length %>% as.vector
  model <- t.test(x = m1,y = m2)
  temp <- data.frame(SynCom = sync,pval =model$p.value)
  Res <- rbind(Res,temp)
}
Res$padj <- Res$pval %>% p.adjust(method = "fdr")
Res$Type <- rep("Test",nrow(Res))
colnames(Res)[1] <- "Condition"

palette_dropout_duo <- c("#000000","#E69F00")
names(palette_dropout_duo) <- c("Full","DOV")
df$Condition <- df$Condition %>% factor(levels = c("a","c","d","35-member"))
p <- chibi.boxplot(Map = df,x_val = "SynCom",y_val = "Length",col_val = "SynCom",
                   facet_formula = "Condition",
                   style = "mix",mpalette = palette_dropout_duo,
                   median_colored_as_points = T,
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
  #geom_vline(xintercept = c(1.5,2.5,3.5,4.5,5.5),size = 0.6, color = "#D9D9D9") +
  theme(legend.position = "none",
        axis.text.x = element_blank(),axis.title.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank()
  ) 

Res$Text <- Res$padj %>% format.pval(digits = 3)
Res$Condition <- Res$Condition %>% factor(levels = c("a","c","d","35-member"))
p <- p + geom_text(data = Res,aes(x = 1.5, y = 7.5,label = Text))+
  scale_y_continuous(limits = c(0,7.5))


oh.save.pdf(p = p,outname = "primaryroot_elongation_variovoraxdropout_biotic.pdf",outdir = "../figures/",
            width = 18,height = 20)
