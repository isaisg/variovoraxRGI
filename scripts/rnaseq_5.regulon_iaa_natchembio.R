library(ohchibi)
library(DESeq2)
library(ggrepel)
library(egg)

set.seed(seed = 130816)
setwd('/home/isai/Documents/results/github/variovoraxRGI/scripts/')


## Tri partite dataset ##
Dat_rnaseq_comb_do <- readRDS(file = "../cleandata/dat_rnaseq_plants_tripartitecl28cl14_variovoraxdropout.RDS")

Dat_raw <- Dat_rnaseq_comb_do$Dat_rnaseq

Dat_sub <- Dat_raw %>% subset.Dataset(subset = Experiment == "Fernanda",drop = T,clean = T)

#Create design for the transformed object
dds<-DESeqDataSetFromMatrix(countData = Dat_sub$Tab,colData = Dat_sub$Map,
                            design =~ Syncom)


dds <- DESeq(dds)
vsd <- vst(dds, blind=FALSE)

Tab_vsd_z <- assay(vsd) %>% t %>% scale %>% t %>% as.matrix

## Load the regulon from the natchembio article
df_auxin_regulon <- read.table(file = "../rawdata/rnaseq_regulon_iaa_natchembio.tsv",header = T,sep = "\t")
df_auxin_regulon <- with(df_auxin_regulon,order(-logFC)) %>%
  df_auxin_regulon[.,]
plot(df_auxin_regulon$logFC)
df_auxin_regulon %>% subset(logFC >=2) %>% dim
sig_genes <- df_auxin_regulon %>% head(n = 13) %$% AGI.code %>%
  as.character

Tab_vsd_z_sub <- match(sig_genes,rownames(Tab_vsd_z)) %>%
  Tab_vsd_z[.,] %>% na.omit

### Scale the data
melted_Tab_vsd_Z <- Tab_vsd_z_sub %>% melt
colnames(melted_Tab_vsd_Z) <- c("genes","Sample_Id","value")
melted_Tab_vsd_Z <- merge(melted_Tab_vsd_Z,Dat_sub$Map,by = "Sample_Id")

mdf <- dcast(melted_Tab_vsd_Z,genes ~Syncom,value.var = "value",fun.aggregate = mean) %>%
  melt
mdf$Ranking <- 1:12 %>% factor
mdf$Names <- df_auxin_regulon$Symbol[1:12] %>% as.character
mdf$Names <- mdf$Names %>% gsub(pattern = "YDK1",replacement = "GH3.2")
mdf$Names <- mdf$Names %>% factor(levels = mdf$Names %>% unique) 

mdf <- mdf %>% subset(variable == "CL28" | variable == "CL28_CL14") %>% droplevels


#Compute median values to compare
temp1 <- mdf %>% subset(variable == "CL28") %$% value %>%
  quantile 
done <- data.frame(variable = "CL28",
                   median = temp1[3],
                   q1 = temp1[2],q3 = temp1[4],row.names = NULL
)
temp2 <- mdf %>% subset(variable == "CL28_CL14") %$% value %>%
  quantile 
dtwo <- data.frame(variable = "CL28_CL14",
                   median = temp2[3],
                   q1 = temp2[2],q3 = temp2[4],row.names = NULL
)

dquant <- rbind(done,dtwo)
mdf_tri <- mdf
dquant_tri <- dquant

#### Compute significance of the trends ####
#Original comparison
Tab_vsd_z_sub <- match(sig_genes,rownames(Tab_vsd_z)) %>%
  Tab_vsd_z[.,] %>% na.omit
melted_Tab_vsd_Z <- Tab_vsd_z_sub %>% melt
colnames(melted_Tab_vsd_Z) <- c("genes","Sample_Id","value")
melted_Tab_vsd_Z <- merge(melted_Tab_vsd_Z,Dat_sub$Map,by = "Sample_Id")
mdf <- dcast(melted_Tab_vsd_Z,genes ~Syncom,value.var = "value",fun.aggregate = mean) %>%
  melt
mdf <- mdf %>% subset(variable == "CL28" | variable == "CL28_CL14") %>% droplevels
v1 <- mdf %>% subset(variable == "CL28") %$% value %>% mean
v2 <- mdf %>% subset(variable == "CL28_CL14") %$% value %>% mean
original_diff <- v1-v2

#Boot 1000 times
boot <- 10000
genes <- rownames(Tab_vsd_z)
df_boot <- NULL
for(i in 1:boot){
  cat("Boot: ",i,"\n",sep="")
  sub_genes <- sample(x = genes,size = 12,replace = F)
  Tab_vsd_z_sub <- match(sub_genes,rownames(Tab_vsd_z)) %>%
    Tab_vsd_z[.,]
  melted_Tab_vsd_Z <- Tab_vsd_z_sub %>% melt
  colnames(melted_Tab_vsd_Z) <- c("genes","Sample_Id","value")
  melted_Tab_vsd_Z <- merge(melted_Tab_vsd_Z,Dat_sub$Map,by = "Sample_Id")
  mdf <- dcast(melted_Tab_vsd_Z,genes ~Syncom,value.var = "value",fun.aggregate = mean) %>%
    melt(id.vars = "genes")
  mdf <- mdf %>% subset(variable == "CL28" | variable == "CL28_CL14") %>% droplevels
  v1 <- mdf %>% subset(variable == "CL28") %$% value %>% mean
  v2 <- mdf %>% subset(variable == "CL28_CL14") %$% value %>% mean
  temp_df <- data.frame(DOV = v1,Full = v2, Diff = v1-v2, boot = i)
  df_boot <- rbind(df_boot,temp_df)
}

ggplot(data = df_boot,aes(Diff)) + geom_histogram() + theme_ohchibi() +
  geom_vline(xintercept =original_diff ,color = "red")

pval_tri <- (which(abs(df_boot$Diff) > abs(original_diff)) %>% length)/boot


# Dropout experiment
Dat_rnaseq_comb_do <- readRDS(file = "../cleandata/dat_rnaseq_plants_tripartitecl28cl14_variovoraxdropout.RDS")

Dat_raw <- Dat_rnaseq_comb_do$Dat_rnaseq

Dat_sub <- Dat_raw %>% subset.Dataset(subset = Experiment == "VarioDO",drop = T,clean = T)

Dat_sub$Map$Syncom <- Dat_sub$Map$Syncom %>% factor(levels = c("NB","Full","DOB","DOV","DOVB"))

#Create design for the transformed object
dds<-DESeqDataSetFromMatrix(countData = Dat_sub$Tab,colData = Dat_sub$Map,
                            design =~ Syncom)


dds <- DESeq(dds)
vsd <- vst(dds, blind=FALSE)

Tab_vsd_z <- assay(vsd) %>% t %>% scale %>% t %>% as.matrix

## Extract the top 
df_auxin_regulon <- read.table(file = "../rawdata/rnaseq_regulon_iaa_natchembio.tsv",header = T,sep = "\t")
df_auxin_regulon <- with(df_auxin_regulon,order(-logFC)) %>%
  df_auxin_regulon[.,]
plot(df_auxin_regulon$logFC)
df_auxin_regulon %>% subset(logFC >=2) %>% dim
sig_genes <- df_auxin_regulon %>% head(n = 13) %$% AGI.code %>%
  as.character

Tab_vsd_z_sub <- match(sig_genes,rownames(Tab_vsd_z)) %>%
  Tab_vsd_z[.,] %>% na.omit

### Scale the data
melted_Tab_vsd_Z <- Tab_vsd_z_sub %>% melt
colnames(melted_Tab_vsd_Z) <- c("genes","Sample_Id","value")
melted_Tab_vsd_Z <- merge(melted_Tab_vsd_Z,Dat_sub$Map,by = "Sample_Id")

mdf <- dcast(melted_Tab_vsd_Z,genes ~Syncom,value.var = "value",fun.aggregate = mean) %>%
  melt
mdf$Ranking <- 1:12 %>% factor
mdf$Names <- df_auxin_regulon$Symbol[1:12] %>% as.character
mdf$Names <- mdf$Names %>% gsub(pattern = "YDK1",replacement = "GH3.2")
mdf$Names <- mdf$Names %>% factor(levels = mdf$Names %>% unique)


mdf <- mdf %>% subset(variable == "Full" | variable == "DOV") %>% droplevels
#Compute median values to compare
temp1 <- mdf %>% subset(variable == "DOV") %$% value %>%
  quantile 
done <- data.frame(variable = "DOV",
                   median = temp1[3],
                   q1 = temp1[2],q3 = temp1[4],row.names = NULL
)
temp2 <- mdf %>% subset(variable == "Full") %$% value %>%
  quantile 
dtwo <- data.frame(variable = "Full",
                   median = temp2[3],
                   q1 = temp2[2],q3 = temp2[4],row.names = NULL
)

dquant <- rbind(done,dtwo)
dquant_drop <- dquant
mdf_drop <- mdf


#### Compute significance of the trends ####
#Original comparison
Tab_vsd_z_sub <- match(sig_genes,rownames(Tab_vsd_z)) %>%
  Tab_vsd_z[.,] %>% na.omit
melted_Tab_vsd_Z <- Tab_vsd_z_sub %>% melt
colnames(melted_Tab_vsd_Z) <- c("genes","Sample_Id","value")
melted_Tab_vsd_Z <- merge(melted_Tab_vsd_Z,Dat_sub$Map,by = "Sample_Id")
mdf <- dcast(melted_Tab_vsd_Z,genes ~Syncom,value.var = "value",fun.aggregate = mean) %>%
  melt
mdf <- mdf %>% subset(variable == "Full" | variable == "DOV") %>% droplevels
v1 <- mdf %>% subset(variable == "DOV") %$% value %>% mean
v2 <- mdf %>% subset(variable == "Full") %$% value %>% mean
original_diff <- v1-v2

#Boot 1000 times
boot <- 10000
genes <- rownames(Tab_vsd_z)
df_boot <- NULL
for(i in 1:boot){
  cat("Boot: ",i,"\n",sep="")
  sub_genes <- sample(x = genes,size = 12,replace = F)
  Tab_vsd_z_sub <- match(sub_genes,rownames(Tab_vsd_z)) %>%
    Tab_vsd_z[.,]
  melted_Tab_vsd_Z <- Tab_vsd_z_sub %>% melt
  colnames(melted_Tab_vsd_Z) <- c("genes","Sample_Id","value")
  melted_Tab_vsd_Z <- merge(melted_Tab_vsd_Z,Dat_sub$Map,by = "Sample_Id")
  mdf <- dcast(melted_Tab_vsd_Z,genes ~Syncom,value.var = "value",fun.aggregate = mean) %>%
    melt(id.vars = "genes")
  mdf <- mdf %>% subset(variable == "Full" | variable == "DOV") %>% droplevels
  v1 <- mdf %>% subset(variable == "DOV") %$% value %>% mean
  v2 <- mdf %>% subset(variable == "Full") %$% value %>% mean
  temp_df <- data.frame(DOV = v1,Full = v2, Diff = v1-v2, boot = i)
  df_boot <- rbind(df_boot,temp_df)
}

ggplot(data = df_boot,aes(Diff)) + geom_histogram() + theme_ohchibi() +
  geom_vline(xintercept =original_diff ,color = "red")

pval_drop <- (which(abs(df_boot$Diff) > abs(original_diff)) %>% length)/boot



### Merge both structures ###
mdf_drop$Experiment <- rep("Dropout",nrow(mdf_drop))
mdf_tri$Experiment <- rep("Tripartite",nrow(mdf_tri))
mdf <- rbind(mdf_drop,mdf_tri)

dquant_drop$Experiment <- rep("Dropout",nrow(dquant_drop))
dquant_tri$Experiment <- rep("Tripartite",nrow(dquant_tri))
dquant <- rbind(dquant_drop,dquant_tri)


###Construct two separate 
p_tri <- ggplot(mdf_tri,aes(variable,value)) +  
  facet_grid(.~Experiment,space = "free",scales = "free") +
  geom_pointrange(data = dquant_tri,
                  aes(x = variable ,y = median,ymin =q1,ymax = q3),inherit.aes = F,
                  size = 3,color = "red") +
  stat_summary(aes(group = 1),geom = 'smooth', alpha = 0.2, fill = 'red', color = 'red',
               fun.data = median_hilow, fun.args = list(conf.int = 1)) +
  geom_line(aes(group = (Ranking)),alpha = 1, colour = "black") +
  geom_point(size = 12,fill = "#414141",shape = 21,alpha = 0.1) + 
  geom_text_repel(data = subset(mdf,variable == "CL28"),aes(label = Names),nudge_x = -0.2,segment.size = 0.1,size = 10) +
  geom_text_repel(data = subset(mdf,variable == "CL28_CL14"),aes(label = Names),nudge_x = 0.2,segment.size = 0.1,size = 10) +
  theme_ohchibi() +
  #geom_point(data = dquant,aes(x = variable, y= median,fill = variable),shape = 21) +
  scale_y_continuous(limits = c(-0.8,1.7)) +
  ylab(label = "z-score") +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.title.x = element_blank(),
        panel.background  =element_rect(fill = "#fee5cb"),
        strip.text.x = element_blank()
        
  ) 


p_drop <- ggplot(mdf_drop,aes(variable,value)) +  
  facet_grid(.~Experiment,space = "free",scales = "free") +
  geom_pointrange(data = dquant_drop,
                  aes(x = variable ,y = median,ymin =q1,ymax = q3),inherit.aes = F,
                  size = 3,color = "red") +
  stat_summary(aes(group = 1),geom = 'smooth', alpha = 0.2, fill = 'red', color = 'red',
               fun.data = median_hilow, fun.args = list(conf.int = 1)) +
  geom_line(aes(group = (Ranking)),alpha = 1, colour = "black") +
  geom_point(size = 12,fill = "#414141",shape = 21,alpha = 0.1) + 
  geom_text_repel(data = subset(mdf,variable == "Full"),aes(label = Names),nudge_x = 0.2,segment.size = 0.1,size = 10) +
  geom_text_repel(data = subset(mdf,variable == "DOV"),aes(label = Names),nudge_x = -0.2,segment.size = 0.1,size = 10) +
  theme_ohchibi() +
  #geom_point(data = dquant,aes(x = variable, y= median,fill = variable),shape = 21) +
  scale_y_continuous(limits = c(-0.8,1.7)) +
  ylab(label = "z-score") +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.title.x = element_blank(),
        panel.background  =element_rect(fill = "#f8e9e6"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        strip.text.x = element_blank()
  ) 

### Create composition of both figures ###
composition <- egg::ggarrange(p_tri,p_drop,nrow =1)

#Save figure
oh.save.pdf(p = composition,outname = "rnaseq_iaaregulon_tripartitecl28cl14_variovoraxdropout.pdf",
            outdir = "../figures/",height = 20,width = 18)
