library(ohchibi)

setwd('/home/isai/Documents/results/github/variovoraxRGI/scripts/')
set.seed(130816)
dir.create("../figures")
dir.create("../cleandata")

#Read dataset
Dat_amplicon<-readRDS(file = "../cleandata/dat_amplicon_variovoraxdropout_useq97.RDS")

#Remove OTUs
Dat <- Dat_amplicon$RawCounts
Dat <- Dat %>% subset.Dataset(syncom == "Full" | syncom == "DOV",drop = T,clean = T)
Dat <- Dat$Tab %>% rownames %>% grep(pattern = "OTU",value = T) %>% remove_taxons(Dat = Dat,taxons = .)

Dat <- Dat %>% subset.Dataset(condition == "control",drop = T,clean = T)

#Remove variovorax from the full
vario_useqs <- c("Sequence_1","Sequence_61","Sequence_78","Sequence_81")

##ALl conditions
Dat_full <- Dat %>% subset.Dataset(subset = syncom == "Full",drop = T,clean = T)
Dat_full$Tab <- Dat_full$Tab[which(!(Dat_full$Tab %>% rownames %in% vario_useqs)),]

Dat_minusv <- Dat %>% subset.Dataset(subset = syncom == "DOV",drop = T,clean = T)
Dat_minusv$Tab <- Dat_minusv$Tab[which(!(Dat_minusv$Tab %>% rownames %in% vario_useqs)),]

Tab <- Dat_full$Tab %>% merge(Dat_minusv$Tab,by = "row.names",all = TRUE)
Tab[is.na(Tab)] <- 0
rownames(Tab) <- Tab$Row.names
Tab <- Tab[,-1] %>% as.matrix
Map <- rbind(Dat_full$Map,Dat_minusv$Map)
Tab <- sweep(Tab,2,colSums(Tab),`/`)
dat <- create_dataset(Tab = Tab,Map = Map)


distfun <- function(x,method) vegan::vegdist(x = x, method = "bray")

paleta <- c("#56B4E9","#000000","#E69F00")
names(paleta) <- c("NB","Full","DOV")
mycap <- oh.cap(Tab = dat$Tab,Map = dat$Map,
                formula = "fraction + syncom +fraction:syncom +condition +Condition(rep)",
                distfun = distfun,perms = 5000,sqrt = T)

p <-chibi.cap(list_ohpco = mycap,comp_a = "CAP1",comp_b = "CAP2",
              col_val = "syncom",shape = "fraction",mypch = 21,size = 14,
              alpha = 1,alpha_shape_background = 0,lines_zero = T) +
  scale_color_manual(values = paleta) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),legend.position = "none") +
  stat_ellipse()


p <- ggplot(data = mycap$Map_cap,
            aes(CAP1,CAP2,fill = syncom,shape = fraction)) +
  geom_hline(yintercept = 0, size = 0.5, 
             color = "#D9D9D9", linetype = "longdash") +
  geom_vline(xintercept = 0,size = 0.5, color = "#D9D9D9",linetype = "longdash") +
  geom_point(size = 30) + scale_shape_manual(values = c(22,21,24)) +
  theme_ohchibi() + scale_fill_manual(values = paleta) +
  theme(panel.grid.major.y   = element_blank(),
        panel.grid.minor.y     = element_blank(),legend.position = "none")  
p <-  p+ stat_ellipse(aes(group = fraction)) +
  xlab(label = "CAP1(45.86%)") + ylab(label = "CAP2(19.45%)")

oh.save.pdf(p = p,outname = "amplicon_variovoraxdropout_cap_dropout_control.pdf",outdir = "../figures/",width = 20,height = 20)

mycap <- oh.cap(Tab = dat$Tab,Map = dat$Map,
                formula = "Condition(fraction) + syncom +  +Condition(rep)",
                distfun = distfun,perms = 5000,sqrt = T)

##PERMANOVA model ####
#Permanova
Tab_bray <- distfun(t(dat$Tab))
mypermanova <- adonis(Tab_bray ~  fraction + syncom  + fraction:syncom,
                      data = dat$Map,strata = dat$Map$rep,permutations = 5000)

mypermanova

capture.output(file = "../figures/amplicon_variovoraxdropout_summary_permanova_models.doc",
               append = F,print(mypermanova))


#Take the control root only
Dat_full <- Dat %>% subset.Dataset(subset = syncom == "Full" & condition == "control" & fraction == "Root",drop = T,clean = T)
Dat_full$Tab <- Dat_full$Tab[which(!(Dat_full$Tab %>% rownames %in% vario_useqs)),]

Dat_minusv <- Dat %>% subset.Dataset(subset = syncom == "DOV" & condition == "control"  & fraction == "Root",drop = T,clean = T)
Dat_minusv$Tab <- Dat_minusv$Tab[which(!(Dat_minusv$Tab %>% rownames %in% vario_useqs)),]


Tab <- Dat_full$Tab %>% merge(Dat_minusv$Tab,by = "row.names",all = TRUE)
Tab[is.na(Tab)] <- 0
rownames(Tab) <- Tab$Row.names
Tab <- Tab[,-1] %>% as.matrix
Map <- rbind(Dat_full$Map,Dat_minusv$Map)
Tab <- sweep(Tab,2,colSums(Tab),`/`)
dat <- create_dataset(Tab = Tab,Map = Map)

distfun <- function(x,method) vegan::vegdist(x = x, method = "bray")


mycap <- oh.cap(Tab = dat$Tab,Map = dat$Map,
                formula = "syncom  +Condition(rep)",
                distfun = distfun,perms = 5000,sqrt = T)

p <-chibi.cap(list_ohpco = mycap,comp_a = "CAP1",comp_b = "MDS1",
              col_val = "syncom",mypch = 21,size = 14,
              alpha = 1,alpha_shape_background = 0) 
##PERMANOVA model ####
#Permanova
Tab_bray <- distfun(t(dat$Tab))
mypermanova <- adonis(Tab_bray ~   syncom ,
                      data = dat$Map,strata = dat$Map$rep,permutations = 5000)

mypermanova

