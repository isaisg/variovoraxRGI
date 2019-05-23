library(ohchibi)


setwd('/home/isai/Documents/results/github/variovoraxRGI/scripts/')
set.seed(130816)

df <- read.table(file = "../cleandata/primaryroot_elongation_monoassociation_normalized.tsv",header = T,quote = "",comment.char = "")

res <- NULL            
for(st in levels(df$Strain)){
  med <- df %>% subset(Strain == st) %>% droplevels %$% 
    NormLength %>% median
  res <- data.frame(Strain = st, median = med) %>% 
    rbind(res,.)
}

srti <- res %>% subset(median  < 3) %$% Strain %>% as.character

meta <- read.table(file = "../rawdata/metadata_97useq_to_ids.tsv",header = T,sep = "\t")
meta$Freezer_Id <- meta$Freezer_Id %>% gsub(pattern = "_.*",replacement = "")
which(!(srti %in% meta$Freezer_Id)) %>% srti[.]


meta <- match(srti,meta$Freezer_Id) %>% meta[.,] %>% droplevels

meta <- meta[,c(2,4,5,6)]
write.table(x = meta,file = "../cleandata/primaryroot_elongation_monoassociation_34_RGIs.tsv",append = F,quote = F,sep = "\t",row.names = F,col.names = T)
