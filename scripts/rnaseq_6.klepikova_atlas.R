library(ohchibi)
library(paletteer)
library(corrplot)

set.seed(130816)
setwd('/home/isai/Documents/results/github/variovoraxRGI/scripts/')
dir.create("../cleandata")
dir.create("../figures")

files <- list.files(path = "../rawdata/",
                    pattern = "klepikova.*",full.names = T)

#Read the files
df <- NULL
for(file in files){
  temp <- read.table(file = file,header = T,sep = "\t",
             quote = "",comment.char = "")
  temp <- temp[,2:4]
  colnames(temp) <- c("Tissue","Expression","SD")
  gene_id <- file %>% strsplit(split = "\\/") %>% unlist %>%
    grep(pattern = ".txt",value = T) %>%
    gsub(pattern = ".txt",replacement = "")
  temp$Gene <- rep(gene_id,nrow(temp))
  df <- rbind(df,temp)
  
}


df$Gene <- df$Gene %>% gsub(pattern = "klepikova_",replacement = "")
###Compute correlation
Tab <- acast(data = df,formula = Tissue~Gene,
      fun.aggregate = mean,value.var = "Expression") 
Tab_cor <- cor(Tab)



dendrogram_Tab_cor <- as.dist(Tab_cor) %>% hclust(method = "complete")
order_genes <- dendrogram_Tab_cor$order %>% dendrogram_Tab_cor$labels[.]

# Get upper triangle of the correlation matrix
#Taken from http://www.sthda.com/english/wiki/ggplot2-quick-correlation-matrix-heatmap-r-software-and-data-visualization

get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}

#Read list of gene identifier
df_names <- read.table(file = "../rawdata/rnaseq_18genes_intersection_dropout_tripartite_names.txt",
                       header = T,sep = "\t",check.names = F,quote = "",comment.char = "",row.names = 1)
df_names <- match(rownames(Tab_cor),rownames(df_names)) %>%
  df_names[.,]

found <- df_names$`Primary Gene Symbol` %>% as.character
found <- found %>% gsub(pattern = ".*\\(",replacement = "") %>%
  gsub(pattern = "\\)",replacement = "")
found[7] <- rownames(df_names)[7]
found[8] <- rownames(df_names)[8]
found[11] <- rownames(df_names)[11]
found[15] <- rownames(df_names)[15]

rownames(Tab_cor) <- found
colnames(Tab_cor) <- found
Tab_cor <- Tab_cor[dendrogram_Tab_cor$order,dendrogram_Tab_cor$order]

upper_tri <- get_upper_tri(Tab_cor)
melted <- upper_tri %>% melt(na.rm = T)

colnames(melted) <- c("Gene2","Gene1","value")
p <- ggplot(melted, aes(Gene1, Gene2, fill = value))+
  geom_tile(color = "white") +
  scale_fill_paletteer_c(package = "viridis",
                         palette = "cividis",
                         name = "Pearson\nCorrelation") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 1, 
                                   size = 12, hjust = 1),
        axis.text.y = element_text(size = 12))+
  coord_fixed() +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    legend.justification = c(1, 0),
    legend.position = c(0.6, 0.7),
    legend.direction = "horizontal")+
  guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                               title.position = "top", title.hjust = 0.5))
oh.save.pdf(p = p,outname = "rnaseq_18genes_intersection_tripartitecl28_dov_correlation_klepikova.pdf",
            outdir = "../figures/")
