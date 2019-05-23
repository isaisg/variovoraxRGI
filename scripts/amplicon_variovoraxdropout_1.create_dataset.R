library(ohchibi)
library(gridExtra)

setwd('/home/isai/Documents/results/github/variovoraxRGI/scripts/')
set.seed(130816)
dir.create("../figures")
dir.create("../cleandata")


#Read the matrices files
matrices_files <- list.files(path = "../rawdata/",
                             pattern = "amplicons_variovoraxdropout_matrix*",full.names = T) 
#Use only the full and v
matrices_files <- matrices_files %>% grep(pattern = "otu",value = T,invert = T) 
  #grep(pattern = "matrix_i",value = T,invert = T)

melted <- NULL
for(mat in matrices_files){
  melt1 <- read.table(file = mat,
                      header = T,sep = "\t",quote = "",row.names = 1,check.names = F)%>% 
    as.matrix %>%  melt
  melted <- rbind(melted,melt1)
}
Tab_usearch <- acast(data = melted,Var1~Var2,value.var = "value",fill =T)


#Read the otus
#Use read.am to read the otu table
Tab_otus <- read.am(file = "../rawdata/amplicons_variovoraxdropout_matrix_otu.tsv",format = "qiime",taxonomy = "taxonomy")
Tax_otus<-Tab_otus$Tax

#Create contaminants to remove later
contam_otus <- c(grep(pattern = "chloroplast", as.character(Tax_otus$Taxonomy),ignore.case = TRUE),
                 grep(pattern = "mitochondri", as.character(Tax_otus$Taxonomy),ignore.case = TRUE),
                 grep(pattern = "oomycete", as.character(Tax_otus$Taxonomy),ignore.case = TRUE),
                 grep(pattern = "Root; k__Bacteria; p__NA", as.character(Tax_otus$Taxonomy),ignore.case = TRUE),
                 grep(pattern = "Root; k__NA;", as.character(Tax_otus$Taxonomy),ignore.case = TRUE))

contam_otus <- row.names(Tax_otus)[contam_otus]

# Filter Dataset
Tab_otus <- remove_taxons(Dat = Tab_otus, taxons = contam_otus)
Tab_otus <- clean(Dat = Tab_otus,verbose = TRUE)


#Compute proportion of Mapped to strain 
usearch_counts <- data.frame((colSums(Tab_usearch)/sum(Tab_usearch))*100)
colnames(usearch_counts) <- "Global_abundance"
usearch_counts$Strain <- rownames(usearch_counts)
percentage_mapped_to_strains <- sum(usearch_counts[grep(pattern = "Sequence_",x = usearch_counts$Strain),][,1])
percentage_mapped_to_contaminants <- sum(usearch_counts[grep(pattern = "contaminant",x = usearch_counts$Strain),][,1])
percentage_unmapped <- sum(usearch_counts[grep(pattern = "Unmapped",x = usearch_counts$Strain),][,1])

###Create a table with the mapping data results
df <- data.frame(Category=c("Reads mapped to strains","Reads mapped to contaminants","Reads Unmapped"),Percentage=c(percentage_mapped_to_strains,percentage_mapped_to_contaminants,percentage_unmapped))
grid.table(df,rows=NULL)


#Remove the ones that hit to contaminants
Tab_usearch <- Tab_usearch[,-grep(pattern = "contaminant",colnames(Tab_usearch))]

#Read the metadata
Map <- read.table(file = "../rawdata/amplicons_variovoraxdropout_metadata.tsv",
                  header = T,sep = "\t")

##Perform the merging r
shared_samples <- intersect(Map$ID_Matrix,rownames(Tab_usearch))
Map <- match(shared_samples,Map$ID_Matrix) %>% Map[.,] %>% droplevels
Tab_usearch <- match(shared_samples,rownames(Tab_usearch)) %>% Tab_usearch[.,]

#Create columns with percentage of mappability
mappability <- data.frame(rowSums(Tab_usearch[,-grep(pattern = "Unmapped",x = colnames(Tab_usearch))]),
                        Ummapped=Tab_usearch[,grep(pattern = "Unmapped",x = colnames(Tab_usearch))],
                        Total=rowSums(Tab_usearch))
colnames(mappability)[1] <- "Mapped"

mappability <- mappability[match(Map$ID_Matrix,rownames(mappability)),]
Map <- cbind(Map,mappability)
Map$Unmapped_prop <- Map$Ummapped/Map$Total

Map$syncom <- Map$syncom %>% as.character %>% 
  gsub(pattern = "I_FULL$",replacement = "I_Full") %>%
  gsub(pattern = "I_MINUS_V",replacement = "I_DOV") %>%
  gsub(pattern = "^FULL$",replacement = "Full") %>%
  gsub(pattern = "MINUS_V",replacement = "DOV") %>% 
  factor(levels = c("I_Full","I_DOV","Full","DOV"))


#Plot level of unmappability
#chibi.boxplot(Map = Map,x_val = "syncom",y_val = "Unmapped_prop",col_val = "fraction",mypch = 21)

#Remove the unmapped
Tab_usearch <- Tab_usearch[,-which(colnames(Tab_usearch)=="Unmapped")]

#Now Map the OTUS
Tab_otus <- t(Tab_otus$Tab)
Tab_otus <- Tab_otus[match(Map$OTU_Header,rownames(Tab_otus)),]


#Perform the merging
rownames(Tab_otus) <- Map$ID_Matrix

Tab_merged<-cbind(Tab_usearch,Tab_otus)
Tab_merged<-t(Tab_merged)
Tab_merged<-Tab_merged[,match(Map$ID_Matrix,colnames(Tab_merged))]

Tab_merged[is.na(Tab_merged)]<-0

rownames(Map) <- Map$ID_Matrix
#Tab_merged<-Tab_merged[,-which(is.na(colnames(Tab_merged)))]
Map <- Map[match(colnames(Tab_merged),Map$ID_Matrix),]
#Read the Taxonomic profile based on Mothur
Tax_usearch<-read.table(file = "../rawdata/df_97useq_taxonomy.tsv",header=F,sep="\t")

#Modify the structure to contain the proper qiime format needed by AMOR
temp_new<-NULL
for(element in Tax_usearch$V2){
  temp_vec <- unlist(strsplit(x = element,split = ";"))
  qiime_format <- paste("Root; k__",temp_vec[1],"; p__",temp_vec[2],"; c__",temp_vec[3],"; o__",temp_vec[4],"; f__",temp_vec[5],"; g__",temp_vec[6],sep="")
  temp_new <- c(temp_new,qiime_format)
}

Tax_usearch <- data.frame(ID=Tax_usearch$V1,Taxonomy=temp_new)
rownames(Tax_usearch) <- Tax_usearch$ID

#Merge taxonomies
Tax_merged <- rbind(Tax_usearch,Tax_otus)

#Adjust the taxonomy
Tax_merged<- Tax_merged[match(rownames(Tab_merged),rownames(Tax_merged)),]


#Create the dataset object
Dat <- create_dataset(Tab = Tab_merged,Map = Map,Tax = Tax_merged)


#Compute the Depth per Sample
Dat$Map$Depth <- colSums(Dat$Tab)

#Compute the Log Depth per Sample
Dat$Map$LDepth <-log(colSums(Dat$Tab))

##Plot usable reads
Dat <- clean(Dat)

#Remove sapmles with less than 1000 reads
Dat_raw <- remove_samples(Dat,row.names(Dat$Map)[ Dat$Map$Depth < 1000 ])

#Rarefaction
set.seed(130816)
#Remove the OTUs
Dat_temp <- remove_taxons(Dat = Dat,taxons = as.character(Dat$Tax$ID[ grep(pattern = "OTU",x = Dat$Tax$ID) ]))
Dat_temp$Map$Counts_Rar<-colSums(Dat_temp$Tab)
Dat_temp <- subset(Dat_temp, Counts_Rar > 1000, drop = TRUE, clean = TRUE)
Dat_rar <- rarefaction(x = Dat_temp,sample = 1000)
Dat_rar <- clean(Dat_rar)

###Relative abundance 
Dat_ra <- remove_samples(Dat,row.names(Dat$Map)[ Dat$Map$Depth < 1000 ])
Dat_ra$Tab <- scale(x = Dat_ra$Tab,center = F,scale = colSums(Dat_ra$Tab))

##Save the raw , rarefied and relative abundance datasets into a global structure
Dat_amplicon_combi <- list(RawCounts=Dat_raw,Rarefied=Dat_rar,RelativeAbundance=Dat_ra)
saveRDS(object = Dat_amplicon_combi,file =,"../cleandata/dat_amplicon_variovoraxdropout_useq97.RDS" )


