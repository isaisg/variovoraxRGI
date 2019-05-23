library(ohchibi)
library(gridExtra)
setwd('/home/isai/Documents/results/github/variovoraxRGI/scripts/')


outdir <- "../cleandata/"
figures <- "../figures/"

dir.create(figures)
dir.create(outdir)

#Read mapping table
Tab_usearch <- read.table("../rawdata/amplicon_4stresses_matrix_useq97.tsv",header=T,row.names = 1,sep="\t",check.names = F)
#Use read.am to read the otu table
Tab_otus <- read.am(file = "../rawdata/amplicon_4stresses_matrix_otu.tsv",format = "qiime",taxonomy = "taxonomy")

Tax_otus<-Tab_otus$Tax

#Compute proportion of read that mapped to strains 
usearch_counts <- data.frame((colSums(Tab_usearch)/sum(Tab_usearch))*100)
colnames(usearch_counts) <- "Global_abundance"
usearch_counts$Strain <- rownames(usearch_counts)
percentage_mapped_to_strains <- sum(usearch_counts[grep(pattern = "Sequence_",x = usearch_counts$Strain),][,1])
percentage_mapped_to_contaminants <- sum(usearch_counts[grep(pattern = "contaminant",x = usearch_counts$Strain),][,1])
percentage_unmapped <- sum(usearch_counts[grep(pattern = "Unmapped",x = usearch_counts$Strain),][,1])

###Create a table with the mapping data results
df <- data.frame(Category=c("Reads mapped to strains","Reads mapped to contaminants","Reads Unmapped"),Percentage=c(percentage_mapped_to_strains,percentage_mapped_to_contaminants,percentage_unmapped))
grid.table(df,rows=NULL)

#
total <- 81.023833+13.959454
#Percentage of read mapping to strains
(81.023833/total)*100


#Remove the ones that hit to contaminants
Tab_usearch <- Tab_usearch[,-grep(pattern = "contaminant",colnames(Tab_usearch))]


#Read metadata file
metadata_file <- "../rawdata/amplicon_4stresses_metadata.tsv"
Map <- read.table(metadata_file,header=T,sep="\t")

##Perform the merging 
Tab_usearch <- Tab_usearch[match(Map$ID_Matrix,rownames(Tab_usearch)),]

#Remove the unmapped
Tab_usearch <- Tab_usearch[,-which(colnames(Tab_usearch)=="Unmapped")]

#Now Map the OTUS
Tab_otus <- t(Tab_otus$Tab)
Tab_otus <- Tab_otus[match(Map$OTU_Header,rownames(Tab_otus)),]


Tab_otus <- Tab_otus[match(Map$OTU_Header,rownames(Tab_otus)),]

#Perform the merging
rownames(Tab_otus) <- Map$ID_Matrix
Tab_merged <- cbind(Tab_usearch,Tab_otus)
rownames(Map) <- Map$ID_Matrix
Tab_merged <- t(Tab_merged)
Tab_merged <- Tab_merged[,match(Map$ID_Matrix,colnames(Tab_merged))]

#Read the Taxonomic profile based on Mothur
Tax_usearch <- read.table(file = "../rawdata/df_97useq_taxonomy.tsv",header=F,sep="\t")

#Modify the structure to contain the proper qiime format needed by AMOR
temp_new<-NULL
for(element in Tax_usearch$V2){
  temp_vec<-unlist(strsplit(x = element,split = ";"))
  qiime_format<-paste("Root; k__",temp_vec[1],"; p__",temp_vec[2],"; c__",temp_vec[3],"; o__",temp_vec[4],"; f__",temp_vec[5],"; g__",temp_vec[6],sep="")
  temp_new<-c(temp_new,qiime_format)
}

Tax_usearch <- data.frame(ID=Tax_usearch$V1,Taxonomy=temp_new)
rownames(Tax_usearch) <- Tax_usearch$ID

#Merge taxonomies
Tax_merged <- rbind(Tax_usearch,Tax_otus)

#Adjust the taxonomy
Tax_merged <- Tax_merged[match(rownames(Tab_merged),rownames(Tax_merged)),]

##Change the values of metadata to make it more comprehensible
Map$type <-  factor(gsub(pattern = "NP",replacement = "NoPlant",x = gsub(pattern = "SC",replacement = "Plant",x = Map$type)))
Map$typebyTissue <- Map$typebyTissue %>% gsub(pattern = "Agar_NP",replacement = "AgarNoPlant") %>% 
  gsub(pattern = "Agar_SC",replacement = "AgarPlant") %>%
  gsub(pattern = "Inoculum_Inoculum",replacement = "Inoculum") %>%
  gsub(pattern = "Root_SC",replacement = "Root") %>%
  gsub(pattern = "Shoot_SC",replacement = "Shoot") %>%
  factor(levels = c("Inoculum","AgarNoPlant","AgarPlant","Root","Shoot"))


#Create the dataset object
Dat <- create_dataset(Tab = Tab_merged,Map = Map,Tax = Tax_merged)


#Compute the Depth per Sample
Dat$Map$Depth <- colSums(Dat$Tab)

#Compute the Log Depth per Sample
Dat$Map$LDepth <-log(colSums(Dat$Tab))


#Readjust the levels
Dat$Map$condition <-   Dat$Map$condition %>% gsub(pattern = "5.5",replacement = "5.5pH") %>%
  gsub(pattern = "7",replacement = "7pH") %>%
  gsub(pattern = "8.2",replacement = "8.2pH") %>%  
  gsub(pattern = "P$",replacement = "Pi") %>%
  gsub(pattern = "S$",replacement = "NaCl") %>%
  gsub(pattern = "T$",replacement = "C") %>%
  factor(
    levels=c("Inoculum","0Pi","10Pi","30Pi","50Pi","100Pi","1000Pi","50NaCl","100NaCl","150NaCl","200NaCl",
             "5.5pH","7pH","8.2pH","10C","21C","31C")
  )


##Plot usable reads
Dat <- clean(Dat)

p <- chibi.boxplot(Map = Dat$Map,x_val = "Tissue",y_val = "Depth",col_val ="type")+ scale_y_log10() 

p <- chibi.boxplot(Map = Dat$Map,x_val = "Tissue",y_val = "Depth",col_val ="experiment")+ scale_y_log10()

#Remove samples with less than 1000 reads
Dat_raw <- remove_samples(Dat,row.names(Dat$Map)[ Dat$Map$Depth < 1000 ])

#Rarefaction
set.seed(130817)
#Remove the OTUs
#Dat_temp <- remove_taxons(Dat = Dat,taxons = as.character(Dat$Tax$ID[ grep(pattern = "OTU",x = Dat$Tax$ID) ]))
Dat_temp <- Dat
Dat_temp$Map$Counts_Rar <- colSums(Dat_temp$Tab)
Dat_temp <- subset(Dat_temp, Counts_Rar > 1000, drop = TRUE, clean = TRUE)
Dat_rar <- rarefaction(x = Dat_temp,sample = 1000)
Dat_rar <- clean(Dat_rar)


###Relative abundance 
Dat_ra <- remove_samples(Dat,row.names(Dat$Map)[ Dat$Map$Depth < 1000 ])
Dat_ra$Tab<-scale(x = Dat_ra$Tab,center = F,scale = colSums(Dat_ra$Tab))


##Save the raw , rarefied and relative abundance datasets into a global structure
Dat_amplicon <- list(RawCounts=Dat_raw,Rarefied=Dat_rar,RelativeAbundance=Dat_ra)
filename <- paste(outdir,"../cleandata/dat_amplicon_4stresses_useq97.RDS",sep = "")
saveRDS(object = Dat_amplicon,file = filename)

