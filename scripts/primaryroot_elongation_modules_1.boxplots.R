library(ohchibi)
library(emmeans)
library(egg)

set.seed(130816)

setwd('/home/isai/Documents/results/github/variovoraxRGI/scripts/')
dir.create("../cleandata/")
dir.create("../figures/")

#Read rawdata
df <- read.table(file = "../rawdata/dataS3_primaryroot_elongation_modules.tsv",
                 header = T,sep = "\t")

#Adjust factors

df$SynCom <- df$SynCom %>% 
  factor(levels = c("NB","Full","A","B","C","D",
                    "AB","AC","AD","BC","BD","CD",
                    "ABC","ABD","ACD","BCD",
                    "R1","R2","R3","R4"))

df$Diversity <- df$Diversity %>% 
  factor(levels = c("NB","Full","Single","Pairs","Trios","Random"))

m1_em <- df %>% aov(formula = Length ~ SynCom + Rep,data = .) %>% 
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

Res_em$Diversity <- match(Res_em$SynCom, df$SynCom) %>% df$Diversity[.]

#Create boxplot
p <- chibi.boxplot(Map = df,x_val = "SynCom",y_val = "Length",
                   style = "open",
                   facet_formula = "Diversity",size_median = 6,
                   size_point = 0,stroke_point = 0,size_axis_text.y = 35,
                   size_axis_text.x = 30,median_color = "red",
                   size_axis_title.y = 40,strip_text_size = 30,
                   size_legend_text = 30,legend_proportion_size = 4) +
  scale_shape_manual(values = 21:22) + ylab(label = "Main Root Elongation (cm)") + 
  xlab(label = "SynCom") +
  geom_sina(alpha = 0.15, color = "#414141",
            stroke = 0,size = 4) +
  geom_text(data = Res_em,
            aes(x = SynCom,y = 8,label = Letters),
            inherit.aes = F,size = 7.5,family ="Arial") +
  coord_cartesian(expand = TRUE) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank())  +
  scale_y_continuous(limits = c(0,8)) 


#Add some rects for the NB and the module A
nb_quantiles <- df %>% subset(SynCom == "NB") %$% Length %>% quantile
df_seg <- data.frame(lower = nb_quantiles[2],upper = nb_quantiles[4])
p_quant <- p +  geom_rect(data = df_seg,
                          mapping = aes(ymin = lower,ymax = upper,xmin = 0,xmax = Inf),
                          inherit.aes = F,alpha = 0.2,color = NA,fill = "#56B4E9")

a_quantiles <- df %>% subset(SynCom == "Full") %$% Length %>% quantile
df_a <- data.frame(lower = a_quantiles[2],upper = a_quantiles[4])

p_quant <- p_quant +  geom_rect(data = df_a,
                                mapping = aes(ymin = lower,ymax = upper,xmin = 0,xmax = Inf),
                                inherit.aes = F,alpha = 0.1,color = NA,fill = "#000000",alpha =0.3 )



#Create checkerboard panel
mat <- c(rep(0,4),rep(1,4),c(1,0,0,0),c(0,1,0,0),c(0,0,1,0),c(0,0,0,1),
         c(1,1,0,0),c(1,0,1,0),c(1,0,0,1),c(0,1,1,0),c(0,1,0,1),c(0,0,1,1)) %>%
  matrix(data = .,nrow = 4,byrow = F) 
rownames(mat) <- c("A","B","C","D")  
colnames(mat) <- c("NB","Full","A","B","C","D","AB","AC","AD","BC","BD","CD")
mat <- mat %>% melt
mat$value <- factor(mat$value)
mat$Var1 <- mat$Var1 %>% factor(levels = c("D","C","B","A"))
colnames(mat) <- c("Module","SynCom","value")

mat$Diversity <- rep("NB",nrow(mat))
mat$Diversity[mat$SynCom %>% grep(pattern = "^[A-Z]$")] <- "Single"
mat$Diversity[mat$SynCom %>% grep(pattern = "^[A-Z][A-Z]$")] <- "Pairs"
mat$Diversity[mat$SynCom %>% grep(pattern = "^[A-Z][A-Z][A-Z]$")] <- "Trios"
mat$Diversity[mat$SynCom %>% grep(pattern = "Full")] <-"Full"
mat$Diversity[mat$SynCom %>% grep(pattern = "NB")] <-"NB"

mat$Diversity <- mat$Diversity %>% factor(levels = c("NB","Full","Single","Pairs","Trios","Random"))


p_mat <- ggplot(data = mat,mapping = aes(SynCom,Module, fill = value)) + 
  geom_tile()  + facet_grid(.~Diversity,space = "free",scales = "free") +
  theme_ohchibi() +
  scale_fill_manual(values = c("#ffffff00","black")) +
  theme(axis.title.x = element_blank(),
        strip.text.x = element_blank(),
        axis.ticks.x = element_blank(),axis.text.x = element_blank(),
        panel.grid  =  element_blank(),legend.position = "none",
        panel.background = element_blank(),
        panel.grid.major.y = element_blank())  +
  ylab(label = "Module") +
  geom_hline(yintercept = c(1.5,2.5,3.5),size = 0.6, color = "#D9D9D9") + 
  geom_vline(xintercept = c(1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5,10.5,11.5,12.5),
             size = 0.6, color = "#D9D9D9")  +
  coord_cartesian(expand = TRUE)  



composition <- egg::ggarrange(ncol = 2,nrow = 2,
                              p,p_quant,
                              p_mat,
                              heights = c(1,0.2),debug = F)

oh.save.pdf(p = composition,outname = "primaryroot_elongation_modules.pdf",
            outdir = "../figures/")

#Save the figures to arrange a panel figure afterwards
mplots <- list(p_checker = p_mat,p_boxplot = p, p_quant = p_quant)
saveRDS(object = mplots,file  = "../figures/primaryroot_elongation_modules_plots.RDS")






