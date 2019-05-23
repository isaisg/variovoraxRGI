library(ohchibi)
library(emmeans)
library(egg)

set.seed(130816)
setwd('/home/isai/Documents/results/github/variovoraxRGI/scripts/')

p_combi <- readRDS("../figures/primaryroot_elongation_modules_plots.RDS")
p_mono <- readRDS("../figures/primaryroot_elongation_monoassociation_plots.RDS")

p_checker <- p_combi$p_checker
p_combi_boxplot <- p_combi$p_boxplot
p_combi_boxplot_quant <- p_combi$p_quant

composition <- egg::ggarrange(nrow = 2,ncol = 2,
                              p_combi_boxplot,p_mono,
                              p_checker,
                              heights = c(1,0.15))
oh.save.pdf(p = composition,
            outname = "primaryroot_elongation_combinatorial_monoassociation_composition.pdf",
            outdir = "../figures/",width = 40,height = 20)

composition <- egg::ggarrange(nrow = 2,ncol = 3,
                              p_combi_boxplot,p_combi_boxplot_quant,p_mono,
                              p_checker,
                              heights = c(1,0.15))
oh.save.pdf(p = composition,
            outname = "primaryroot_elongation_combinatorial_monoassociation_composition_quant.pdf",
            outdir = "../figures/",width = 40,height = 20)
