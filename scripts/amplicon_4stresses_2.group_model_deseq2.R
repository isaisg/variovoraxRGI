library(ohchibi)
library(DESeq2)


setwd('/home/isai/Documents/results/github/variovoraxRGI/scripts/')


set.seed(seed = 130816)
##Phosphate
Dat <- readRDS(file = "../cleandata/dat_amplicon_4stresses_useq97.RDS")
Dat_raw <- Dat$RawCounts

#Create the dds object
Dat_raw$Map$group <- paste(Dat_raw$Map$typebyTissue,Dat_raw$Map$condition,sep = "_") %>%
  factor
dds <- DESeqDataSetFromMatrix(countData = Dat_raw$Tab,colData = Dat_raw$Map,
                              design =~ Rep + experiment + group)

#Run the model
dds <- DESeq(object = dds)

############# Root vs Agar ####################
#Phosphate
res <- dds %>% results(contrast = c("group","Root_0Pi","AgarPlant_0Pi")) %>%
  as.data.frame
res$Contrast <- rep("Root_0Pi_vs_AgarPlant_0Pi",nrow(res))
res$Id <- rownames(res)
res_pi0 <- res

res <- dds %>% results(contrast = c("group","Root_10Pi","AgarPlant_10Pi")) %>%
  as.data.frame
res$Contrast <- rep("Root_10Pi_vs_AgarPlant_10Pi",nrow(res))
res$Id <- rownames(res)
res_pi10 <- res

res <- dds %>% results(contrast = c("group","Root_30Pi","AgarPlant_30Pi")) %>%
  as.data.frame
res$Contrast <- rep("Root_30Pi_vs_AgarPlant_30Pi",nrow(res))
res$Id <- rownames(res)
res_pi30 <- res

res <- dds %>% results(contrast = c("group","Root_50Pi","AgarPlant_50Pi")) %>%
  as.data.frame
res$Contrast <- rep("Root_50Pi_vs_AgarPlant_50Pi",nrow(res))
res$Id <- rownames(res)
res_pi50 <- res

res <- dds %>% results(contrast = c("group","Root_100Pi","AgarPlant_100Pi")) %>%
  as.data.frame
res$Contrast <- rep("Root_100Pi_vs_AgarPlant_100Pi",nrow(res))
res$Id <- rownames(res)
res_pi100 <- res

res <- dds %>% results(contrast = c("group","Root_1000Pi","AgarPlant_1000Pi")) %>%
  as.data.frame
res$Contrast <- rep("Root_1000Pi_vs_AgarPlant_1000Pi",nrow(res))
res$Id <- rownames(res)
res_pi1000 <- res

#Salinity
res <- dds %>% results(contrast = c("group","Root_50NaCl","AgarPlant_50NaCl")) %>%
  as.data.frame
res$Contrast <- rep("Root_50NaCl_vs_AgarPlant_50NaCl",nrow(res))
res$Id <- rownames(res)
res_nacl50 <- res

res <- dds %>% results(contrast = c("group","Root_100NaCl","AgarPlant_100NaCl")) %>%
  as.data.frame
res$Contrast <- rep("Root_100NaCl_vs_AgarPlant_100NaCl",nrow(res))
res$Id <- rownames(res)
res_nacl100 <- res

res <- dds %>% results(contrast = c("group","Root_150NaCl","AgarPlant_150NaCl")) %>%
  as.data.frame
res$Contrast <- rep("Root_150NaCl_vs_AgarPlant_150NaCl",nrow(res))
res$Id <- rownames(res)
res_nacl150 <- res

res <- dds %>% results(contrast = c("group","Root_200NaCl","AgarPlant_200NaCl")) %>%
  as.data.frame
res$Contrast <- rep("Root_200NaCl_vs_AgarPlant_200NaCl",nrow(res))
res$Id <- rownames(res)
res_nacl200 <- res

#pH
res <- dds %>% results(contrast = c("group","Root_5.5pH","AgarPlant_5.5pH")) %>%
  as.data.frame
res$Contrast <- rep("Root_5.5pH_vs_AgarPlant_5.5pH",nrow(res))
res$Id <- rownames(res)
res_ph5 <- res

res <- dds %>% results(contrast = c("group","Root_7pH","AgarPlant_7pH")) %>%
  as.data.frame
res$Contrast <- rep("Root_7pH_vs_AgarPlant_7pH",nrow(res))
res$Id <- rownames(res)
res_ph7 <- res

res <- dds %>% results(contrast = c("group","Root_8.2pH","AgarPlant_8.2pH")) %>%
  as.data.frame
res$Contrast <- rep("Root_8.2pH_vs_AgarPlant_8.2pH",nrow(res))
res$Id <- rownames(res)
res_ph8 <- res

#Temperature
res <- dds %>% results(contrast = c("group","Root_10C","AgarPlant_10C")) %>%
  as.data.frame
res$Contrast <- rep("Root_10C_vs_AgarPlant_10C",nrow(res))
res$Id <- rownames(res)
res_t10 <- res

res <- dds %>% results(contrast = c("group","Root_21C","AgarPlant_21C")) %>%
  as.data.frame
res$Contrast <- rep("Root_21C_vs_AgarPlant_21C",nrow(res))
res$Id <- rownames(res)
res_t21 <- res

res <- dds %>% results(contrast = c("group","Root_31C","AgarPlant_31C")) %>%
  as.data.frame
res$Contrast <- rep("Root_31C_vs_AgarPlant_31C",nrow(res))
res$Id <- rownames(res)
res_t31 <- res


Res_root <- rbind(res_pi0,res_pi10,res_pi30,res_pi50,res_pi100,res_pi1000,
                  res_nacl50,res_nacl100,res_nacl150,res_nacl200,
                  res_ph5,res_ph7,res_ph8,
                  res_t10,res_t21,res_t31)

############# Shoot vs Agar ####################
#Phosphate
res <- dds %>% results(contrast = c("group","Shoot_0Pi","AgarPlant_0Pi")) %>%
  as.data.frame
res$Contrast <- rep("Shoot_0Pi_vs_AgarPlant_0Pi",nrow(res))
res$Id <- rownames(res)
res_pi0 <- res

res <- dds %>% results(contrast = c("group","Shoot_10Pi","AgarPlant_10Pi")) %>%
  as.data.frame
res$Contrast <- rep("Shoot_10Pi_vs_AgarPlant_10Pi",nrow(res))
res$Id <- rownames(res)
res_pi10 <- res

res <- dds %>% results(contrast = c("group","Shoot_30Pi","AgarPlant_30Pi")) %>%
  as.data.frame
res$Contrast <- rep("Shoot_30Pi_vs_AgarPlant_30Pi",nrow(res))
res$Id <- rownames(res)
res_pi30 <- res

res <- dds %>% results(contrast = c("group","Shoot_50Pi","AgarPlant_50Pi")) %>%
  as.data.frame
res$Contrast <- rep("Shoot_50Pi_vs_AgarPlant_50Pi",nrow(res))
res$Id <- rownames(res)
res_pi50 <- res

res <- dds %>% results(contrast = c("group","Shoot_100Pi","AgarPlant_100Pi")) %>%
  as.data.frame
res$Contrast <- rep("Shoot_100Pi_vs_AgarPlant_100Pi",nrow(res))
res$Id <- rownames(res)
res_pi100 <- res

res <- dds %>% results(contrast = c("group","Shoot_1000Pi","AgarPlant_1000Pi")) %>%
  as.data.frame
res$Contrast <- rep("Shoot_1000Pi_vs_AgarPlant_1000Pi",nrow(res))
res$Id <- rownames(res)
res_pi1000 <- res

#Salinity
res <- dds %>% results(contrast = c("group","Shoot_50NaCl","AgarPlant_50NaCl")) %>%
  as.data.frame
res$Contrast <- rep("Shoot_50NaCl_vs_AgarPlant_50NaCl",nrow(res))
res$Id <- rownames(res)
res_nacl50 <- res

res <- dds %>% results(contrast = c("group","Shoot_100NaCl","AgarPlant_100NaCl")) %>%
  as.data.frame
res$Contrast <- rep("Shoot_100NaCl_vs_AgarPlant_100NaCl",nrow(res))
res$Id <- rownames(res)
res_nacl100 <- res

res <- dds %>% results(contrast = c("group","Shoot_150NaCl","AgarPlant_150NaCl")) %>%
  as.data.frame
res$Contrast <- rep("Shoot_150NaCl_vs_AgarPlant_150NaCl",nrow(res))
res$Id <- rownames(res)
res_nacl150 <- res

res <- dds %>% results(contrast = c("group","Shoot_200NaCl","AgarPlant_200NaCl")) %>%
  as.data.frame
res$Contrast <- rep("Shoot_200NaCl_vs_AgarPlant_200NaCl",nrow(res))
res$Id <- rownames(res)
res_nacl200 <- res

#pH
res <- dds %>% results(contrast = c("group","Shoot_5.5pH","AgarPlant_5.5pH")) %>%
  as.data.frame
res$Contrast <- rep("Shoot_5.5pH_vs_AgarPlant_5.5pH",nrow(res))
res$Id <- rownames(res)
res_ph5 <- res

res <- dds %>% results(contrast = c("group","Shoot_7pH","AgarPlant_7pH")) %>%
  as.data.frame
res$Contrast <- rep("Shoot_7pH_vs_AgarPlant_7pH",nrow(res))
res$Id <- rownames(res)
res_ph7 <- res

res <- dds %>% results(contrast = c("group","Shoot_8.2pH","AgarPlant_8.2pH")) %>%
  as.data.frame
res$Contrast <- rep("Shoot_8.2pH_vs_AgarPlant_8.2pH",nrow(res))
res$Id <- rownames(res)
res_ph8 <- res

#Temperature
res <- dds %>% results(contrast = c("group","Shoot_10C","AgarPlant_10C")) %>%
  as.data.frame
res$Contrast <- rep("Shoot_10C_vs_AgarPlant_10C",nrow(res))
res$Id <- rownames(res)
res_t10 <- res

res <- dds %>% results(contrast = c("group","Shoot_21C","AgarPlant_21C")) %>%
  as.data.frame
res$Contrast <- rep("Shoot_21C_vs_AgarPlant_21C",nrow(res))
res$Id <- rownames(res)
res_t21 <- res

res <- dds %>% results(contrast = c("group","Shoot_31C","AgarPlant_31C")) %>%
  as.data.frame
res$Contrast <- rep("Shoot_31C_vs_AgarPlant_31C",nrow(res))
res$Id <- rownames(res)
res_t31 <- res

Res_shoot <- rbind(res_pi0,res_pi10,res_pi30,res_pi50,res_pi100,res_pi1000,
                   res_nacl50,res_nacl100,res_nacl150,res_nacl200,
                   res_ph5,res_ph7,res_ph8,
                   res_t10,res_t21,res_t31)

############# Root vs Shoot ####################
#Phosphate
res <- dds %>% results(contrast = c("group","Root_0Pi","Shoot_0Pi")) %>%
  as.data.frame
res$Contrast <- rep("Root_0Pi_vs_Shoot_0Pi",nrow(res))
res$Id <- rownames(res)
res_pi0 <- res

res <- dds %>% results(contrast = c("group","Root_10Pi","Shoot_10Pi")) %>%
  as.data.frame
res$Contrast <- rep("Root_10Pi_vs_Shoot_10Pi",nrow(res))
res$Id <- rownames(res)
res_pi10 <- res

res <- dds %>% results(contrast = c("group","Root_30Pi","Shoot_30Pi")) %>%
  as.data.frame
res$Contrast <- rep("Root_30Pi_vs_Shoot_30Pi",nrow(res))
res$Id <- rownames(res)
res_pi30 <- res

res <- dds %>% results(contrast = c("group","Root_50Pi","Shoot_50Pi")) %>%
  as.data.frame
res$Contrast <- rep("Root_50Pi_vs_Shoot_50Pi",nrow(res))
res$Id <- rownames(res)
res_pi50 <- res

res <- dds %>% results(contrast = c("group","Root_100Pi","Shoot_100Pi")) %>%
  as.data.frame
res$Contrast <- rep("Root_100Pi_vs_Shoot_100Pi",nrow(res))
res$Id <- rownames(res)
res_pi100 <- res

res <- dds %>% results(contrast = c("group","Root_1000Pi","Shoot_1000Pi")) %>%
  as.data.frame
res$Contrast <- rep("Root_1000Pi_vs_Shoot_1000Pi",nrow(res))
res$Id <- rownames(res)
res_pi1000 <- res

#Salinity
res <- dds %>% results(contrast = c("group","Root_50NaCl","Shoot_50NaCl")) %>%
  as.data.frame
res$Contrast <- rep("Root_50NaCl_vs_Shoot_50NaCl",nrow(res))
res$Id <- rownames(res)
res_nacl50 <- res

res <- dds %>% results(contrast = c("group","Root_100NaCl","Shoot_100NaCl")) %>%
  as.data.frame
res$Contrast <- rep("Root_100NaCl_vs_Shoot_100NaCl",nrow(res))
res$Id <- rownames(res)
res_nacl100 <- res

res <- dds %>% results(contrast = c("group","Root_150NaCl","Shoot_150NaCl")) %>%
  as.data.frame
res$Contrast <- rep("Root_150NaCl_vs_Shoot_150NaCl",nrow(res))
res$Id <- rownames(res)
res_nacl150 <- res

res <- dds %>% results(contrast = c("group","Root_200NaCl","Shoot_200NaCl")) %>%
  as.data.frame
res$Contrast <- rep("Root_200NaCl_vs_Shoot_200NaCl",nrow(res))
res$Id <- rownames(res)
res_nacl200 <- res

#pH
res <- dds %>% results(contrast = c("group","Root_5.5pH","Shoot_5.5pH")) %>%
  as.data.frame
res$Contrast <- rep("Root_5.5pH_vs_Shoot_5.5pH",nrow(res))
res$Id <- rownames(res)
res_ph5 <- res

res <- dds %>% results(contrast = c("group","Root_7pH","Shoot_7pH")) %>%
  as.data.frame
res$Contrast <- rep("Root_7pH_vs_Shoot_7pH",nrow(res))
res$Id <- rownames(res)
res_ph7 <- res

res <- dds %>% results(contrast = c("group","Root_8.2pH","Shoot_8.2pH")) %>%
  as.data.frame
res$Contrast <- rep("Root_8.2pH_vs_Shoot_8.2pH",nrow(res))
res$Id <- rownames(res)
res_ph8 <- res

#Temperature
res <- dds %>% results(contrast = c("group","Root_10C","Shoot_10C")) %>%
  as.data.frame
res$Contrast <- rep("Root_10C_vs_Shoot_10C",nrow(res))
res$Id <- rownames(res)
res_t10 <- res

res <- dds %>% results(contrast = c("group","Root_21C","Shoot_21C")) %>%
  as.data.frame
res$Contrast <- rep("Root_21C_vs_Shoot_21C",nrow(res))
res$Id <- rownames(res)
res_t21 <- res

res <- dds %>% results(contrast = c("group","Root_31C","Shoot_31C")) %>%
  as.data.frame
res$Contrast <- rep("Root_31C_vs_Shoot_31C",nrow(res))
res$Id <- rownames(res)
res_t31 <- res

Res_rs <- rbind(res_pi0,res_pi10,res_pi30,res_pi50,res_pi100,res_pi1000,
                res_nacl50,res_nacl100,res_nacl150,res_nacl200,
                res_ph5,res_ph7,res_ph8,
                res_t10,res_t21,res_t31)

Res_root$OverallContrast <- rep("RootvsAgar",nrow(Res_root))
Res_shoot$OverallContrast <- rep("ShootvsAgar",nrow(Res_shoot))
Res_rs$OverallContrast <- rep("RootvsShoot",nrow(Res_rs))

Res <- rbind(Res_root,Res_shoot,Res_rs)



saveRDS(object = dds,file = "../cleandata/dds_4stresses_group_model.RDS")

Res <- Res[,c(8,1:7,9)]

#Read the 4map table
df_modules <- read.table(file = "../rawdata/amplicon_4stresses_map_useq2cluster_abcd.tsv",header = T,sep = "\t")
colnames(df_modules)[1] <- "Id"
Res <- merge(Res,df_modules, by = "Id",all.x = TRUE)
Res <- Res$Id %>% grep(pattern = "OTU",invert = T) %>%
  Res[.,]
colnames(Res)[1] <- "USeq"
Res$USeq <- Res$USeq %>% gsub(pattern = "Sequence_",replacement = "USeq")
colnames(Res)[10] <-  "Module"
write.table(x = Res,file = "../rawdata/dataS2_fraction_enrichments.tsv",
            append = F,quote = F,sep = "\t",row.names = F,col.names = T)
