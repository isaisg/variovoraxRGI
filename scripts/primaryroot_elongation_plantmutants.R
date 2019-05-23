library(ohchibi)
library(emmeans)
library(paletteer)
library(Rmisc)

setwd('/home/isai/Documents/results/github/variovoraxRGI/scripts/')
set.seed(130816)
dir.create("../cleandata/")
dir.create("../figures/")

#Define palette
paleta <- paletteer_d(package = "RSkittleBrewer",palette = "original",n = 5,direction = 1)[2:5]
paleta[4] <- "#CCBF14"


df <- read.table("../rawdata/dataS15_primaryroot_elongation_plantmutants.tsv",
                 header=T,sep = "\t",comment.char = "",quote = "")


#Split Full and duo
df_controls <- df %>% subset(Strain == "NB" | Strain == "Full" |Strain == "CL28_CL14") %>% droplevels


df$Strain <- factor(df$Strain,levels=  c("NB","Full","CL28_CL14","IAA","ACC","CL28","DOV","MF48","CL52","MF6"))
df$Genotype_MCP<-factor(df$Genotype_MCP,levels=c("col","col_mcp","axr","axr_mcp"))

###Standarize against NB ##### 
res <- NULL
for (geno in levels(df$Genotype_MCP)){
  df_sub <- df %>% subset(Genotype_MCP == geno)
  nf <- df %>% subset(Genotype_MCP == geno & Strain == "NB") %$% Length %>% mean
  df_sub$nLength <- (df_sub %>% subset(Genotype_MCP == geno) %$% Length) / nf
  res <- rbind(res,df_sub)
}

res_NB <- res %>% 
  subset(Strain == "NB") %>% droplevels

res <- res %>% 
  subset(Strain != "NB") %>% droplevels

res_Full_duo <- res %>%
  subset(Strain == "Full" | Strain == "CL28_CL14") %>%
  droplevels

res <- res %>% subset(Strain != "Full") %>%
  subset(Strain != "CL28_CL14") %>% 
  droplevels

### Test using anova framework
Res_em <- NULL
for(st in levels(res$Strain)){
  df_sub <- res %>% subset(Strain == st)
  m1 <-  df_sub %>% aov(formula = nLength ~ Genotype_MCP + Rep)  
  m1_em <- emmeans(m1, specs = "Genotype_MCP") %>% CLD() 
  m1_em$.group <- m1_em$.group %>% gsub(pattern = " ",replacement = "")
  colnames(m1_em)[which(colnames(m1_em) == ".group")] <- "group"
  m1_em$Strain <- rep(st,nrow(m1_em))
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

#Adjust levels Genotype
res$Genotype_MCP <- res$Genotype_MCP %>%
  factor(levels = c("col","axr","col_mcp","axr_mcp"))
Res_em$Genotype_MCP <- Res_em$Genotype_MCP %>%
  factor(levels = c("col","axr","col_mcp","axr_mcp"))

#Subset only CL28 and droput
res_sub <- res %>% subset(Strain == "CL28" | Strain == "DOV") %>%
  droplevels
Res_em_sub <- Res_em %>% subset(Strain == "CL28" | Strain == "DOV") %>%
  droplevels


paleta <- rep("red",4)
p <- chibi.sina(Map = res_sub,x_val = "Genotype_MCP",y_val = " nLength",
                facet_formula = "Strain",size_point = 10,size_bar = 2,
                df_stats = Res_em_sub,mean_var = "emmean",
                ymin_var = "lower.CL",ymax_var = "upper.CL",
                size_axis_text.x = 0,show_points = T,font_family = "Arial",
                size_axis_text.y = 30,size_axis_title.x = 35,size_axis_title.y = 35,
                size_legend_text = 35,strip_text_size = 35,legend_proportion_size = 4,
                mpalette = paleta) + 
  #geom_vline(xintercept = c(1.5,2.5,3.5,4.5,5.5),size = 0.3,color = "#D9D9D9") +
  theme(legend.position = "none",axis.title.x = element_blank())

p <- p + geom_text(data = Res_em_sub,
                   aes(x = Genotype_MCP,y = 1.2,label = Letters),
                   inherit.aes = F,size = 8,family ="Arial") +
  ylab(label = "Standarized root Length") 


res_controls_plot <- rbind(res_NB,res_Full_duo)
df_segments <- summarySE(data = res_controls_plot,measurevar = "nLength",groupvars = "Strain")
df_segments$lower <- df_segments$nLength - df_segments$ci
df_segments$upper <- df_segments$nLength + df_segments$ci

#Add no bacteria segment
df_segment_NB <- df_segments %>% subset(Strain == "NB")
df_segment_NB <- rbind(df_segment_NB,df_segment_NB)
df_segment_NB$Strain <- df_segment_NB$Strain %>% as.character
df_segment_NB$Strain <- c("CL28","DOV") %>% factor


p <- p  + geom_hline(yintercept = c(1),size = 0.6, color = "#56B4E9", linetype="dashed") 


#Add the presece of Full segment only
df_segment_vario <- df_segments %>% subset(Strain == "Full")
df_segment_vario$Strain <- df_segment_vario$Strain %>% as.character
df_segment_vario$Strain <- c("DOV") %>%
  factor

p <- p +  geom_rect(data = df_segment_vario,
                    mapping = aes(ymin = lower,ymax = upper,xmin = 0,xmax = Inf),
                    inherit.aes = F,alpha = 0.2,color = NA,fill = "#000000")

#Add the presece of Full segment only
df_segment_CL28 <- df_segments %>% subset(Strain == "CL28_CL14")
df_segment_CL28$Strain <- df_segment_CL28$Strain %>% as.character
df_segment_CL28$Strain <- c("CL28") %>%
  factor

p <- p +  geom_rect(data = df_segment_CL28,
                    mapping = aes(ymin = lower,ymax = upper,xmin = 0,xmax = Inf),
                    inherit.aes = F,alpha = 0.2,color = NA,fill = "#000000")
p <- p + theme(
  panel.grid.major.y = element_blank(),
  panel.grid.minor.y = element_blank()
)

oh.save.pdf(p = p,outdir = "../figures/",outname = "primaryroot_elongation_plantmutants.pdf",width = 20,height = 20)
