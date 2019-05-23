library(ohchibi)
library(DESeq2)

setwd('/home/isai/Documents/results/github/variovoraxRGI/scripts/')
set.seed(130816)

dir.create("../cleandata/")
dir.create("../figures/")


#Read the Matrix files obtained directly from Ombnibus
Tab_1 <- read.table(file = "../rawdata/SynCom_mRNA_counts_Tripartite.txt",skip = 1,header = T,sep = "\t",row.names = 1)
Tab_2 <- read.table(file = "../rawdata/SynCom_mRNA_counts_Dropout.txt",skip = 1,header = T,sep = "\t",row.names = 1)

#Put the ids in relationship with Tab_1
Tab_2 <- Tab_2[match(rownames(Tab_1),rownames(Tab_2)),]

#Merge all the Tabs
Tab <- Tab_1 %>% cbind(Tab_2) 


#Read the metadata
Map <- read.table(file = "../rawdata/rnaseq_metadata.tsv",
                  header = T,sep = "\t",comment.char = "",quote = "")


#Check the intersection between Tab and Map
Map <- Map[which(Map$Sample_Id %in% colnames(Tab)),] %>% droplevels
Tab <- match(Map$Sample_Id,colnames(Tab)) %>% Tab[,.]

#Create the dataset
rownames(Map) <- Map$Sample_Id
Dat_rnaseq_comb_do_raw <- create_dataset(Tab = Tab,Map = Map)

#Read the length of the genes
df_length_genes <- read.table(file = "../rawdata/rnaseq_genelength.tsv",header = T,sep = "\t")

Dat_rnaseq_comb_do <- list(Dat_rnaseq = Dat_rnaseq_comb_do_raw, 
                           df_length_genes = df_length_genes)
saveRDS(object = Dat_rnaseq_comb_do,file = "../cleandata/dat_rnaseq_plants_tripartitecl28cl14_variovoraxdropout.RDS")
