library(ohchibi)
library(paletteer)
library(palettesPM)
library(egg)
library(RColorBrewer)
library(extrafont)
loadfonts(device = "pdf") 
library(emmeans)

### Optimized for the oh.ggsave function
size_legend_text <- 45
size_axis_title <- 20
size_axis_text <- 15
size_legend_title <- 55
legend_proportion_size <- 4
size_points <- 10
strip_text_size <- 25
size_median <- 4


ylim_min <- -1.5
ylim_max <- 1.5
xlim_min <- -1.5
xlim_max <- 1.5

setwd('/home/isai/Documents/results/github/variovoraxRGI/scripts/')

outfile_permanova <- "../figures/amplicon_4stresses_summary_permanova_models.doc"

set.seed(seed = 130816)

#Read the dataset 
Dat_amplicon <- readRDS(file = "../cleandata/dat_amplicon_4stresses_useq97.RDS")
Dat <- Dat_amplicon$RelativeAbundance
distfun <- function(x,method) vegan::vegdist(x = x, method = "bray")

palette_variance <- paletteer_d(package = "dutchmasters",palette = "pearl_earring",11)[c(1:4,7)] %>%
  c(.,"white")
names(palette_variance) <- c("typebyTissue","Genotype","condition",
                             "typebyTissue:condition","typebyTissue:Genotype","Residual")
palette_variance <- palette_variance[c(1,3,4,6)]

colores_phosphate <- c("#ffffe5","#f7fcb9","#addd8e","#41ab5d","#238443","#005a32")
colores_salinity <- c("#fde0dd","#fa9fb5","#dd3497", "#7a0177")
colores_temp <- c("#FFC88C","#FFA33F","#7F6446")
colores_ph <- c("#C3C1EA","#7E78E5","#555466")

### Full dataset for sup figure2
Dat_sub <- Dat %>% 
  subset.Dataset(subset = typebyTissue != "Inoculum",drop = T,clean = T)%>%
  subset.Dataset(typebyTissue != "AgarNoPlant",drop = T,clean = T)
Dat_sub$Map$Technical <- paste(Dat_sub$Map$experiment,Dat_sub$Map$Rep, sep = "_") %>% factor
mcap <- oh.cap(Tab = Dat_sub$Tab,Map = Dat_sub$Map,
               formula = "typebyTissue  + condition + Condition(Technical)",
               distfun = distfun,perms = 5000,sqrt = T)

p_all <- chibi.cap(list_ohpco = mcap,col_val = "typebyTissue",shape_val = "experiment",
                   comp_a = "CAP1",comp_b = "CAP2" ,size =size_points ,alpha=1,
                   size_legend_text = size_legend_text,size_axis_title = size_axis_title,size_axis_line = 2,
                   size_title_text = size_legend_title,size_axis_text =  size_axis_text,
                   font_family = "Arial",legend_proportion_size = 4,lines_zero = T) +
  scale_fill_fraction()+ theme(legend.position = "none") +
  scale_shape_manual(values = 21:24) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  theme(axis.ticks = element_blank(),axis.text.x = element_blank(),axis.text.y = element_blank())


p_all <- ggplot(data = mcap$Map_cap,mapping = aes(CAP1,CAP2)) +
  geom_vline(xintercept = 0,linetype = "dashed",size = 2 , color = "#D9D9D9") +
  geom_hline(yintercept = 0,linetype = "dashed",size = 2 , color = "#D9D9D9") +
  stat_ellipse(aes(group = typebyTissue, fill = typebyTissue),geom = "polygon",alpha = 0.4) +
  geom_point(aes(fill = typebyTissue,shape = experiment),size = 12) +
  scale_fill_fraction() + scale_shape_manual(values = 21:24) +
  theme_ohchibi() +  theme(panel.grid.major.y = element_blank(),
                           panel.grid.minor.y = element_blank()) +
  theme(axis.ticks = element_blank(),axis.text.x = element_blank(),axis.text.y = element_blank()) +
  xlab(label = "CAP1(40.87%)") + ylab(label = "CAP2(20.94%)")

  

Tab_bray <- distfun(t(Dat_sub$Tab))
mypermanova <- adonis(Tab_bray ~  typebyTissue + condition  + typebyTissue:condition,
                      data = Dat_sub$Map,strata = Dat_sub$Map$Technical,permutations = 5000)

mypermanova

Tab_bray <- distfun(t(Dat_sub$Tab))
mypermanova <- adonis(Tab_bray ~  typebyTissue + condition  + typebyTissue:condition + Technical,
                      data = Dat_sub$Map,permutations = 5000)

mypermanova


oh.save.pdf(p = p_all,outname = "amplicon_4stresses_supfigure_cap_fraction_4experiments.pdf",outdir = "../figures/",
            height = 20,width = 20)


####### Salinity #########
Dat_sub <- Dat %>% subset.Dataset(subset = experiment == "salinity"& 
                                    typebyTissue != "Inoculum",drop = T,clean = T) %>%
  subset.Dataset(typebyTissue != "AgarNoPlant",drop = T,clean = T)
mcap <- oh.cap(Tab = Dat_sub$Tab,Map = Dat_sub$Map,
               formula = "typebyTissue  + condition + Condition(Rep)",
               distfun = distfun,perms = 5000,sqrt = T)
#Remove the ticks and the axis labels given the fact all are in the same scale
p_salinity_cap_all <- chibi.cap(list_ohpco = mcap,col_val = "typebyTissue",
                                comp_a = "CAP1",comp_b = "CAP2" ,size =size_points ,alpha=1,
                                size_legend_text = size_legend_text,size_axis_title = size_axis_title,size_axis_line = 2,
                                size_title_text = size_legend_title,size_axis_text =  size_axis_text,
                                font_family = "Arial",legend_proportion_size = 4,lines_zero = T) +
  scale_fill_fraction()+ theme(legend.position = "none") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) + scale_y_continuous(limits = c(-1.5,1.5)) +
  scale_x_continuous(limits = c(-1.5,1.5)) +
  theme(axis.ticks = element_blank(),axis.text.x = element_blank(),axis.text.y = element_blank())

Tab_bray <- distfun(t(Dat_sub$Tab))
mypermanova <- adonis(Tab_bray ~  typebyTissue + condition  + typebyTissue:condition,
                      data = Dat_sub$Map,strata = Dat_sub$Map$Rep,permutations = 5000)

mypermanova
# p_salinity_perm_all <- chibi.permanova(mypermanova = mypermanova,size_legend_text = size_legend_text,size_axis_title = 15,
#                                        size_axis_text = size_axis_text,size_title_text = size_legend_text,legend_proportion_size = 4) + 
#   scale_fill_manual(values = palette_variance) + xlab(label = "Term Model") + 
#   theme(legend.position = "none")  +
#   theme(axis.ticks.x = element_blank(),axis.title.x  = element_blank(),
#         axis.ticks.y = element_line(size = 0.5),axis.text.y = element_text(vjust = 0.5))

### Agar ###
Dat_sub <- Dat %>% subset.Dataset(subset = experiment == "salinity"& 
                                    typebyTissue == "AgarPlant",drop = T,clean = T)
mcap <- oh.cap(Tab = Dat_sub$Tab,Map = Dat_sub$Map,
               formula = "condition + Condition(Rep)",
               distfun = distfun,perms = 5000,sqrt = T)
p_salinity_cap_agar <- chibi.cap(list_ohpco = mcap,col_val = "condition",
                                 comp_a = "CAP1",comp_b = "CAP2" ,size = size_points,alpha=1,
                                 size_legend_text = size_legend_text,size_axis_title = size_axis_title,
                                 size_axis_text =  size_axis_text,
                                 size_axis_line = 2,
                                 size_title_text = size_legend_title,
                                 font_family = "Arial",legend_proportion_size = 4,lines_zero = T) +
  theme(legend.position = "none") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_fill_manual(values = colores_salinity) + scale_y_continuous(limits = c(ylim_min,ylim_max)) +
  scale_x_continuous(limits = c(xlim_min,xlim_max))+
  theme(axis.ticks = element_blank(),axis.text.x = element_blank(),axis.text.y = element_blank())

Tab_bray <- distfun(t(Dat_sub$Tab))
mypermanova <- adonis(Tab_bray ~  condition  ,
                      data = Dat_sub$Map,strata = Dat_sub$Map$Rep,permutations = 5000)

mypermanova
capture.output(file = outfile_permanova,
               append = F,print(mypermanova))

# p_salinity_perm_agar <- chibi.permanova(mypermanova = mypermanova,size_legend_text = size_legend_text,size_axis_title = 15,
#                                         size_axis_text = size_axis_text,size_title_text = size_legend_text,legend_proportion_size = 4) + 
#   scale_fill_manual(values = palette_variance) + xlab(label = "Term Model") + 
#   theme(legend.position = "none")+
#   theme(axis.ticks.x = element_blank(),axis.title.x  = element_blank(),
#         axis.ticks.y = element_line(size = 0.5),axis.text.y = element_text(vjust = 0.5))


### Root ###
Dat_sub <- Dat %>% subset.Dataset(subset = experiment == "salinity"& 
                                    typebyTissue == "Root",drop = T,clean = T)
mcap <- oh.cap(Tab = Dat_sub$Tab,Map = Dat_sub$Map,
               formula = "condition + Condition(Rep)",
               distfun = distfun,perms = 5000,sqrt = T)
p_salinity_cap_root <- chibi.cap(list_ohpco = mcap,col_val = "condition",
                                 comp_a = "CAP1",comp_b = "CAP2" ,size = size_points,alpha=1,
                                 size_legend_text = size_legend_text,size_axis_title = size_axis_title,
                                 size_axis_text =  size_axis_text,size_axis_line = 2,
                                 size_title_text = size_legend_title,
                                 font_family = "Arial",legend_proportion_size = 4,lines_zero = T) +
  theme(legend.position = "none") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_fill_manual(values = colores_salinity)  + scale_y_continuous(limits = c(ylim_min,ylim_max)) +
  scale_x_continuous(limits = c(xlim_min,xlim_max))+
  theme(axis.ticks = element_blank(),axis.text.x = element_blank(),axis.text.y = element_blank())

Tab_bray <- distfun(t(Dat_sub$Tab))
mypermanova <- adonis(Tab_bray ~  condition  ,
                      data = Dat_sub$Map,strata = Dat_sub$Map$Rep,permutations = 5000)

mypermanova
capture.output(file = outfile_permanova ,
               append = T,print(mypermanova))

# p_salinity_perm_root <- chibi.permanova(mypermanova = mypermanova,size_legend_text = size_legend_text,size_axis_title = 15,
#                                         size_axis_text = size_axis_text,size_title_text = size_legend_text,legend_proportion_size = 4) + 
#   scale_fill_manual(values = palette_variance) + xlab(label = "Term Model") + 
#   theme(legend.position = "none")+
#   theme(axis.ticks.x = element_blank(),axis.title.x  = element_blank(),
#         axis.ticks.y = element_line(size = 0.5),axis.text.y = element_text(vjust = 0.5))


### Shoot ###
Dat_sub <- Dat %>% subset.Dataset(subset = experiment == "salinity"& 
                                    typebyTissue == "Shoot",drop = T,clean = T)
mcap <- oh.cap(Tab = Dat_sub$Tab,Map = Dat_sub$Map,
               formula = "condition + Condition(Rep)",
               distfun = distfun,perms = 5000,sqrt = T)
p_salinity_cap_shoot <- chibi.cap(list_ohpco = mcap,col_val = "condition",
                                  comp_a = "CAP1",comp_b = "CAP2" ,size = size_points,alpha=1,
                                  size_legend_text = size_legend_text,size_axis_title = size_axis_title,
                                  size_axis_text =  size_axis_text,size_axis_line = 2,
                                  size_title_text = size_legend_title,
                                  font_family = "Arial",legend_proportion_size = 4,lines_zero = T) +
  theme(legend.position = "none") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_fill_manual(values = colores_salinity)  + scale_y_continuous(limits = c(ylim_min,ylim_max)) +
  scale_x_continuous(limits = c(xlim_min,xlim_max))+
  theme(axis.ticks = element_blank(),axis.text.x = element_blank(),axis.text.y = element_blank())

Tab_bray <- distfun(t(Dat_sub$Tab))
mypermanova <- adonis(Tab_bray ~  condition  ,
                      data = Dat_sub$Map,strata = Dat_sub$Map$Rep,permutations = 5000)

mypermanova
capture.output(file = outfile_permanova,
               append = T,print(mypermanova))

# p_salinity_perm_shoot <- chibi.permanova(mypermanova = mypermanova,size_legend_text = size_legend_text,size_axis_title = 15,
#                                          size_axis_text = size_axis_text,size_title_text = size_legend_text,legend_proportion_size = 4) + 
#   scale_fill_manual(values = palette_variance) + xlab(label = "Term Model") + 
#   theme(legend.position = "none")+
#   theme(axis.ticks.x = element_blank(),axis.title.x  = element_blank(),
#         axis.ticks.y = element_line(size = 0.5),axis.text.y = element_text(vjust = 0.5))
# 

#############################
####### Phosphate #########
Dat_sub <- Dat %>% subset.Dataset(subset = experiment == "Pi"& 
                                    typebyTissue != "Inoculum",drop = T,clean = T) %>%
  subset.Dataset(typebyTissue != "AgarNoPlant",drop = T,clean = T)
mcap <- oh.cap(Tab = Dat_sub$Tab,Map = Dat_sub$Map,
               formula = "typebyTissue  + condition + Condition(Rep)",
               distfun = distfun,perms = 5000,sqrt = T)
p_phosphate_cap_all <- chibi.cap(list_ohpco = mcap,col_val = "typebyTissue",
                                 comp_a = "CAP1",comp_b = "CAP2" ,size = size_points,alpha=1,
                                 size_legend_text = size_legend_text,size_axis_title = size_axis_title,
                                 size_axis_text =  size_axis_text,size_axis_line = 2,
                                 size_title_text = size_legend_title,
                                 font_family = "Arial",legend_proportion_size = 4,lines_zero = T) +
  scale_fill_fraction()+ theme(legend.position = "none") + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+ scale_y_continuous(limits = c(-1.5,1.5)) +
  scale_x_continuous(limits = c(-1.5,1.5))+
  theme(axis.ticks = element_blank(),axis.text.x = element_blank(),axis.text.y = element_blank())

Tab_bray <- distfun(t(Dat_sub$Tab))
mypermanova <- adonis(Tab_bray ~  typebyTissue + condition  + typebyTissue:condition,
                      data = Dat_sub$Map,strata = Dat_sub$Map$Rep,permutations = 5000)

mypermanova
# p_phosphate_perm_all <- chibi.permanova(mypermanova = mypermanova,size_legend_text = size_legend_text,size_axis_title = 15,
#                                         size_axis_text = size_axis_text,size_title_text = size_legend_text,legend_proportion_size = 4) + 
#   scale_fill_manual(values = palette_variance) + xlab(label = "Term Model") + 
#   theme(legend.position = "none")+
#   theme(axis.ticks.x = element_blank(),axis.title.x  = element_blank(),
#         axis.ticks.y = element_line(size = 0.5),axis.text.y = element_text(vjust = 0.5))
# 

### Agar ###
Dat_sub <- Dat %>% subset.Dataset(subset = experiment == "Pi"& 
                                    typebyTissue == "AgarPlant",drop = T,clean = T)
mcap <- oh.cap(Tab = Dat_sub$Tab,Map = Dat_sub$Map,
               formula = "condition + Condition(Rep)",
               distfun = distfun,perms = 5000,sqrt = T)
p_phosphate_cap_agar <- chibi.cap(list_ohpco = mcap,col_val = "condition",
                                  comp_a = "CAP1",comp_b = "CAP2" ,size = size_points,alpha=1,
                                  size_legend_text = size_legend_text,size_axis_title = size_axis_title,
                                  size_axis_text =  size_axis_text,size_axis_line = 2,
                                  size_title_text = size_legend_title,
                                  font_family = "Arial",legend_proportion_size = 4,lines_zero = T) +
  theme(legend.position = "none") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_fill_manual(values = colores_phosphate) + scale_y_continuous(limits = c(ylim_min,ylim_max)) +
  scale_x_continuous(limits = c(xlim_min,xlim_max))+
  theme(axis.ticks = element_blank(),axis.text.x = element_blank(),axis.text.y = element_blank())

Tab_bray <- distfun(t(Dat_sub$Tab))
mypermanova <- adonis(Tab_bray ~  condition  ,
                      data = Dat_sub$Map,strata = Dat_sub$Map$Rep,permutations = 5000)

mypermanova
capture.output(file = outfile_permanova ,
               append = T,print(mypermanova))

# p_phosphate_perm_agar <- chibi.permanova(mypermanova = mypermanova,size_legend_text = size_legend_text,size_axis_title = 15,
#                                          size_axis_text = size_axis_text,size_title_text = size_legend_text,legend_proportion_size = 4) + 
#   scale_fill_manual(values = palette_variance) + xlab(label = "Term Model") + 
#   theme(legend.position = "none")+
#   theme(axis.ticks.x = element_blank(),axis.title.x  = element_blank(),
#         axis.ticks.y = element_line(size = 0.5),axis.text.y = element_text(vjust = 0.5))


### Root ###
Dat_sub <- Dat %>% subset.Dataset(subset = experiment == "Pi"& 
                                    typebyTissue == "Root",drop = T,clean = T)
mcap <- oh.cap(Tab = Dat_sub$Tab,Map = Dat_sub$Map,
               formula = "condition + Condition(Rep)",
               distfun = distfun,perms = 5000,sqrt = T)
p_phosphate_cap_root <- chibi.cap(list_ohpco = mcap,col_val = "condition",
                                  comp_a = "CAP1",comp_b = "CAP2" ,size = size_points,alpha=1,
                                  size_legend_text = size_legend_text,size_axis_title = size_axis_title,
                                  size_axis_text =  size_axis_text,size_axis_line = 2,
                                  size_title_text = size_legend_title,
                                  font_family = "Arial",legend_proportion_size = 4,lines_zero = T) +
  theme(legend.position = "none") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_fill_manual(values = colores_phosphate)  + scale_y_continuous(limits = c(ylim_min,ylim_max)) +
  scale_x_continuous(limits = c(xlim_min,xlim_max))+
  theme(axis.ticks = element_blank(),axis.text.x = element_blank(),axis.text.y = element_blank())

Tab_bray <- distfun(t(Dat_sub$Tab))
mypermanova <- adonis(Tab_bray ~  condition  ,
                      data = Dat_sub$Map,strata = Dat_sub$Map$Rep,permutations = 5000)

mypermanova
capture.output(file = outfile_permanova ,
               append = T,print(mypermanova))

# p_phosphate_perm_root <- chibi.permanova(mypermanova = mypermanova,size_legend_text = size_legend_text,size_axis_title = 15,
#                                          size_axis_text = size_axis_text,size_title_text = size_legend_text,legend_proportion_size = 4) + 
#   scale_fill_manual(values = palette_variance) + xlab(label = "Term Model") + 
#   theme(legend.position = "none")+
#   theme(axis.ticks.x = element_blank(),axis.title.x  = element_blank(),
#         axis.ticks.y = element_line(size = 0.5),axis.text.y = element_text(vjust = 0.5))


### Shoot ###
Dat_sub <- Dat %>% subset.Dataset(subset = experiment == "Pi"& 
                                    typebyTissue == "Shoot",drop = T,clean = T)
mcap <- oh.cap(Tab = Dat_sub$Tab,Map = Dat_sub$Map,
               formula = "condition + Condition(Rep)",
               distfun = distfun,perms = 5000,sqrt = T)
p_phosphate_cap_shoot <- chibi.cap(list_ohpco = mcap,col_val = "condition",
                                   comp_a = "CAP1",comp_b = "CAP2" ,size = size_points,alpha=1,
                                   size_legend_text = size_legend_text,size_axis_title = size_axis_title,
                                   size_axis_text =  size_axis_text,size_axis_line = 2,
                                   size_title_text = size_legend_title,
                                   font_family = "Arial",legend_proportion_size = 4,lines_zero = T) +
  theme(legend.position = "none") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_fill_manual(values = colores_phosphate)  + scale_y_continuous(limits = c(ylim_min,ylim_max)) +
  scale_x_continuous(limits = c(xlim_min,xlim_max))+
  theme(axis.ticks = element_blank(),axis.text.x = element_blank(),axis.text.y = element_blank())

Tab_bray <- distfun(t(Dat_sub$Tab))
mypermanova <- adonis(Tab_bray ~  condition  ,
                      data = Dat_sub$Map,strata = Dat_sub$Map$Rep,permutations = 5000)

mypermanova
capture.output(file = outfile_permanova ,
               append = T,print(mypermanova))

# p_phosphate_perm_shoot <- chibi.permanova(mypermanova = mypermanova,size_legend_text = size_legend_text,size_axis_title = 15,
#                                           size_axis_text = size_axis_text,size_title_text = size_legend_text,legend_proportion_size = 4) + 
#   scale_fill_manual(values = palette_variance) + xlab(label = "Term Model") + 
#   theme(legend.position = "none")+
#   theme(axis.ticks.x = element_blank(),axis.title.x  = element_blank(),
#         axis.ticks.y = element_line(size = 0.5),axis.text.y = element_text(vjust = 0.5))


##########################
####### pH #########
Dat_sub <- Dat %>% subset.Dataset(subset = experiment == "pH"& 
                                    typebyTissue != "Inoculum",drop = T,clean = T) %>%
  subset.Dataset(typebyTissue != "AgarNoPlant",drop = T,clean = T)
mcap <- oh.cap(Tab = Dat_sub$Tab,Map = Dat_sub$Map,
               formula = "typebyTissue  + condition + Condition(Rep)",
               distfun = distfun,perms = 5000,sqrt = T)
p_ph_cap_all <- chibi.cap(list_ohpco = mcap,col_val = "typebyTissue",
                          comp_a = "CAP1",comp_b = "CAP2" ,size = size_points,alpha=1,
                          size_legend_text = size_legend_text,size_axis_title = size_axis_title,
                          size_axis_text =  size_axis_text,size_axis_line = 2,
                          size_title_text = size_legend_title,
                          font_family = "Arial",legend_proportion_size = 4,lines_zero = T) +
  scale_fill_fraction()+ theme(legend.position = "none") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+ scale_y_continuous(limits = c(-1.5,1.5)) +
  scale_x_continuous(limits = c(-1.5,1.5))+
  theme(axis.ticks = element_blank(),axis.text.x = element_blank(),axis.text.y = element_blank())

Tab_bray <- distfun(t(Dat_sub$Tab))
mypermanova <- adonis(Tab_bray ~  typebyTissue + condition  + typebyTissue:condition,
                      data = Dat_sub$Map,strata = Dat_sub$Map$Rep,permutations = 5000)

mypermanova
# p_ph_perm_all <- chibi.permanova(mypermanova = mypermanova,size_legend_text = size_legend_text,size_axis_title = 15,
#                                  size_axis_text = size_axis_text,size_title_text = size_legend_text,legend_proportion_size = 4) + 
#   scale_fill_manual(values = palette_variance) + xlab(label = "Term Model") + 
#   theme(legend.position = "none")+
#   theme(axis.ticks.x = element_blank(),axis.title.x  = element_blank(),
#         axis.ticks.y = element_line(size = 0.5),axis.text.y = element_text(vjust = 0.5))


### Agar ###
Dat_sub <- Dat %>% subset.Dataset(subset = experiment == "pH"& 
                                    typebyTissue == "AgarPlant",drop = T,clean = T)
mcap <- oh.cap(Tab = Dat_sub$Tab,Map = Dat_sub$Map,
               formula = "condition + Condition(Rep)",
               distfun = distfun,perms = 5000,sqrt = T)
p_ph_cap_agar <- chibi.cap(list_ohpco = mcap,col_val = "condition",
                           comp_a = "CAP1",comp_b = "CAP2" ,size = size_points,alpha=1,
                           size_legend_text = size_legend_text,size_axis_title = size_axis_title,
                           size_axis_text =  size_axis_text,size_axis_line = 2,
                           size_title_text = size_legend_title,
                           font_family = "Arial",legend_proportion_size = 4,lines_zero = T) +
  theme(legend.position = "none") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_fill_manual(values = colores_ph)  + scale_y_continuous(limits = c(ylim_min,ylim_max)) +
  scale_x_continuous(limits = c(xlim_min,xlim_max))+
  theme(axis.ticks = element_blank(),axis.text.x = element_blank(),axis.text.y = element_blank())

Tab_bray <- distfun(t(Dat_sub$Tab))
mypermanova <- adonis(Tab_bray ~  condition  ,
                      data = Dat_sub$Map,strata = Dat_sub$Map$Rep,permutations = 5000)

mypermanova
capture.output(file = outfile_permanova ,
               append = T,print(mypermanova))

# p_ph_perm_agar <- chibi.permanova(mypermanova = mypermanova,size_legend_text = size_legend_text,size_axis_title = 15,
#                                   size_axis_text = size_axis_text,size_title_text = size_legend_text,legend_proportion_size = 4) + 
#   scale_fill_manual(values = palette_variance) + xlab(label = "Term Model") + 
#   theme(legend.position = "none")+
#   theme(axis.ticks.x = element_blank(),axis.title.x  = element_blank(),
#         axis.ticks.y = element_line(size = 0.5),axis.text.y = element_text(vjust = 0.5))


### Root ###
Dat_sub <- Dat %>% subset.Dataset(subset = experiment == "pH"& 
                                    typebyTissue == "Root",drop = T,clean = T)
mcap <- oh.cap(Tab = Dat_sub$Tab,Map = Dat_sub$Map,
               formula = "condition + Condition(Rep)",
               distfun = distfun,perms = 5000,sqrt = T)
p_ph_cap_root <- chibi.cap(list_ohpco = mcap,col_val = "condition",
                           comp_a = "CAP1",comp_b = "CAP2" ,size = size_points,alpha=1,
                           size_legend_text = size_legend_text,size_axis_title = size_axis_title,
                           size_axis_text =  size_axis_text,size_axis_line = 2,
                           size_title_text = size_legend_title,
                           font_family = "Arial",legend_proportion_size = 4,lines_zero = T) +
  theme(legend.position = "none") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_fill_manual(values = colores_ph)  + scale_y_continuous(limits = c(ylim_min,ylim_max)) +
  scale_x_continuous(limits = c(xlim_min,xlim_max))+
  theme(axis.ticks = element_blank(),axis.text.x = element_blank(),axis.text.y = element_blank())

Tab_bray <- distfun(t(Dat_sub$Tab))
mypermanova <- adonis(Tab_bray ~  condition  ,
                      data = Dat_sub$Map,strata = Dat_sub$Map$Rep,permutations = 5000)

mypermanova
capture.output(file = outfile_permanova ,
               append = T,print(mypermanova))

# p_ph_perm_root <- chibi.permanova(mypermanova = mypermanova,size_legend_text = size_legend_text,size_axis_title = 15,
#                                   size_axis_text = size_axis_text,size_title_text = size_legend_text,legend_proportion_size = 4) + 
#   scale_fill_manual(values = palette_variance) + xlab(label = "Term Model") + 
#   theme(legend.position = "none")+
#   theme(axis.ticks.x = element_blank(),axis.title.x  = element_blank(),
#         axis.ticks.y = element_line(size = 0.5),axis.text.y = element_text(vjust = 0.5))
# 

### Shoot ###
Dat_sub <- Dat %>% subset.Dataset(subset = experiment == "pH"& 
                                    typebyTissue == "Shoot",drop = T,clean = T)
mcap <- oh.cap(Tab = Dat_sub$Tab,Map = Dat_sub$Map,
               formula = "condition + Condition(Rep)",
               distfun = distfun,perms = 5000,sqrt = T)
p_ph_cap_shoot <- chibi.cap(list_ohpco = mcap,col_val = "condition",
                            comp_a = "CAP1",comp_b = "CAP2" ,size = size_points,alpha=1,
                            size_legend_text = size_legend_text,size_axis_title = size_axis_title,
                            size_axis_text =  size_axis_text,size_axis_line = 2,
                            size_title_text = size_legend_title,
                            font_family = "Arial",legend_proportion_size = 4,lines_zero = T) +
  theme(legend.position = "none") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_fill_manual(values = colores_ph)  + scale_y_continuous(limits = c(ylim_min,ylim_max)) +
  scale_x_continuous(limits = c(xlim_min,xlim_max))+
  theme(axis.ticks = element_blank(),axis.text.x = element_blank(),axis.text.y = element_blank())

Tab_bray <- distfun(t(Dat_sub$Tab))
mypermanova <- adonis(Tab_bray ~  condition  ,
                      data = Dat_sub$Map,strata = Dat_sub$Map$Rep,permutations = 5000)

mypermanova
capture.output(file = outfile_permanova ,
               append = T,print(mypermanova))

# p_ph_perm_shoot <- chibi.permanova(mypermanova = mypermanova,size_legend_text = size_legend_text,size_axis_title = 15,
#                                    size_axis_text = size_axis_text,size_title_text = size_legend_text,legend_proportion_size = 4) + 
#   scale_fill_manual(values = palette_variance) + xlab(label = "Term Model") + 
#   theme(legend.position = "none")+
#   theme(axis.ticks.x = element_blank(),axis.title.x  = element_blank(),
#         axis.ticks.y = element_line(size = 0.5),axis.text.y = element_text(vjust = 0.5))
# 


##########################
####### Temperature #########
Dat_sub <- Dat %>% subset.Dataset(subset = experiment == "T"& 
                                    typebyTissue != "Inoculum",drop = T,clean = T) %>%
  subset.Dataset(typebyTissue != "AgarNoPlant",drop = T,clean = T)
mcap <- oh.cap(Tab = Dat_sub$Tab,Map = Dat_sub$Map,
               formula = "typebyTissue  + condition + Condition(Rep)",
               distfun = distfun,perms = 5000,sqrt = T)
p_t_cap_all <- chibi.cap(list_ohpco = mcap,col_val = "typebyTissue",
                         comp_a = "CAP1",comp_b = "CAP2" ,size = size_points,alpha=1,
                         size_legend_text = size_legend_text,size_axis_title = size_axis_title,
                         size_axis_text =  size_axis_text,size_axis_line = 2,
                         size_title_text = size_legend_title,
                         font_family = "Arial",legend_proportion_size = 4,lines_zero = T) +
  scale_fill_fraction()+ theme(legend.position = "none") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+ scale_y_continuous(limits = c(-1.5,1.5)) +
  scale_x_continuous(limits = c(-1.5,1.5))+
  theme(axis.ticks = element_blank(),axis.text.x = element_blank(),axis.text.y = element_blank())

Tab_bray <- distfun(t(Dat_sub$Tab))
mypermanova <- adonis(Tab_bray ~  typebyTissue + condition  + typebyTissue:condition,
                      data = Dat_sub$Map,strata = Dat_sub$Map$Rep,permutations = 5000)

mypermanova
# p_t_perm_all <- chibi.permanova(mypermanova = mypermanova,size_legend_text = size_legend_text,size_axis_title = 15,
#                                 size_axis_text = size_axis_text,size_title_text = size_legend_text,legend_proportion_size = 4) + 
#   scale_fill_manual(values = palette_variance) + xlab(label = "Term Model") + 
#   theme(legend.position = "none")+
#   theme(axis.ticks.x = element_blank(),axis.title.x  = element_blank(),
#         axis.ticks.y = element_line(size = 0.5),axis.text.y = element_text(vjust = 0.5))


### Agar ###
Dat_sub <- Dat %>% subset.Dataset(subset = experiment == "T"& 
                                    typebyTissue == "AgarPlant",drop = T,clean = T)
mcap <- oh.cap(Tab = Dat_sub$Tab,Map = Dat_sub$Map,
               formula = "condition + Condition(Rep)",
               distfun = distfun,perms = 5000,sqrt = T)
p_t_cap_agar <- chibi.cap(list_ohpco = mcap,col_val = "condition",
                          comp_a = "CAP1",comp_b = "CAP2" ,size = size_points,alpha=1,
                          size_legend_text = size_legend_text,size_axis_title = size_axis_title,
                          size_axis_text =  size_axis_text,size_axis_line = 2,
                          size_title_text = size_legend_title,
                          font_family = "Arial",legend_proportion_size = 4,lines_zero = T) +
  theme(legend.position = "none") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_fill_manual(values = colores_temp)  + scale_y_continuous(limits = c(ylim_min,ylim_max)) +
  scale_x_continuous(limits = c(xlim_min,xlim_max))+
  theme(axis.ticks = element_blank(),axis.text.x = element_blank(),axis.text.y = element_blank())

Tab_bray <- distfun(t(Dat_sub$Tab))
mypermanova <- adonis(Tab_bray ~  condition  ,
                      data = Dat_sub$Map,strata = Dat_sub$Map$Rep,permutations = 5000)

mypermanova
capture.output(file = outfile_permanova ,
               append = T,print(mypermanova))

# p_t_perm_agar <- chibi.permanova(mypermanova = mypermanova,size_legend_text = size_legend_text,size_axis_title = 15,
#                                  size_axis_text = size_axis_text,size_title_text = size_legend_text,legend_proportion_size = 4) + 
#   scale_fill_manual(values = palette_variance) + xlab(label = "Term Model") + 
#   theme(legend.position = "none")+
#   theme(axis.ticks.x = element_blank(),axis.title.x  = element_blank(),
#         axis.ticks.y = element_line(size = 0.5),axis.text.y = element_text(vjust = 0.5))


### Root ###
Dat_sub <- Dat %>% subset.Dataset(subset = experiment == "T"& 
                                    typebyTissue == "Root",drop = T,clean = T)
mcap <- oh.cap(Tab = Dat_sub$Tab,Map = Dat_sub$Map,
               formula = "condition + Condition(Rep)",
               distfun = distfun,perms = 5000,sqrt = T)
p_t_cap_root <- chibi.cap(list_ohpco = mcap,col_val = "condition",
                          comp_a = "CAP1",comp_b = "CAP2" ,size = size_points,alpha=1,
                          size_legend_text = size_legend_text,size_axis_title = size_axis_title,
                          size_axis_text =  size_axis_text,size_axis_line = 2,
                          size_title_text = size_legend_title,
                          font_family = "Arial",legend_proportion_size = 4,lines_zero = T) +
  theme(legend.position = "none") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_fill_manual(values = colores_temp) + scale_y_continuous(limits = c(ylim_min,ylim_max)) +
  scale_x_continuous(limits = c(xlim_min,xlim_max))+
  theme(axis.ticks = element_blank(),axis.text.x = element_blank(),axis.text.y = element_blank())

Tab_bray <- distfun(t(Dat_sub$Tab))
mypermanova <- adonis(Tab_bray ~  condition  ,
                      data = Dat_sub$Map,strata = Dat_sub$Map$Rep,permutations = 5000)

mypermanova
capture.output(file = outfile_permanova ,
               append = T,print(mypermanova))

# p_t_perm_root <- chibi.permanova(mypermanova = mypermanova,size_legend_text = size_legend_text,size_axis_title = 15,
#                                  size_axis_text = size_axis_text,size_title_text = size_legend_text,legend_proportion_size = 4) + 
#   scale_fill_manual(values = palette_variance) + xlab(label = "Term Model") + 
#   theme(legend.position = "none")+
#   theme(axis.ticks.x = element_blank(),axis.title.x  = element_blank(),
#         axis.ticks.y = element_line(size = 0.5),axis.text.y = element_text(vjust = 0.5))


### Shoot ###
Dat_sub <- Dat %>% subset.Dataset(subset = experiment == "T"& 
                                    typebyTissue == "Shoot",drop = T,clean = T)
mcap <- oh.cap(Tab = Dat_sub$Tab,Map = Dat_sub$Map,
               formula = "condition + Condition(Rep)",
               distfun = distfun,perms = 5000,sqrt = T)
p_t_cap_shoot <- chibi.cap(list_ohpco = mcap,col_val = "condition",
                           comp_a = "CAP1",comp_b = "CAP2" ,size = size_points,alpha=1,
                           size_legend_text = size_legend_text,size_axis_title = size_axis_title,
                           size_axis_text =  size_axis_text,size_axis_line = 2,
                           size_title_text = size_legend_title,
                           font_family = "Arial",legend_proportion_size = 4,lines_zero = T) +
  theme(legend.position = "none") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_fill_manual(values = colores_temp)  + scale_y_continuous(limits = c(ylim_min,ylim_max)) +
  scale_x_continuous(limits = c(xlim_min,xlim_max))+
  theme(axis.ticks = element_blank(),axis.text.x = element_blank(),axis.text.y = element_blank())

Tab_bray <- distfun(t(Dat_sub$Tab))
mypermanova <- adonis(Tab_bray ~  condition  ,
                      data = Dat_sub$Map,strata = Dat_sub$Map$Rep,permutations = 5000)

mypermanova
capture.output(file = outfile_permanova,
               append = T,print(mypermanova))


composition <- egg::ggarrange(ncol = 3,nrow = 4,
                               p_phosphate_cap_agar,p_phosphate_cap_root,p_phosphate_cap_shoot,
                               p_salinity_cap_agar,p_salinity_cap_root,p_salinity_cap_shoot,
                               p_ph_cap_agar,p_ph_cap_root,p_ph_cap_shoot,
                               p_t_cap_agar,p_t_cap_root,p_t_cap_shoot,
                              debug = F,byrow = T
 )


oh.save.pdf(p = composition,outname = "amplicon_4stresses_composition_beta_diversity_vertical.pdf",outdir = "../figures/",
             height = 20,width = 15)
 

