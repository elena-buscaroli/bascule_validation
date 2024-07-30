

#-------------------------------------------------------------------------------
# LOADING...
#-------------------------------------------------------------------------------
organ_type <- "Breast" # Breast, Lung, Colorectal
path <- "/Users/azadsadr/"
devtools::load_all( paste(path, "Documents/packages/basilica/", sep = "") )
source(paste(path, "Nextcloud/basilica/scripts/utils_plot.R", sep = ""))
source(paste(path, "Nextcloud/basilica/scripts/utils/plot.signatures.R", sep = ""))

data <- readRDS( paste(path, "Nextcloud/basilica/compiled_fit/", tolower(organ_type), "_data.Rds", sep = "") )

x <- basilica:::merge_clusters(data$x_after, cutoff = 0.8)

x_mapped <- basilica:::convert_dn_names(x, reference_cat = list(SBS=COSMIC_sbs_filt, DBS=COSMIC_dbs), cutoff = 0.8)


#-------------------------------------------------------------------------------
# SET UP THEME
#-------------------------------------------------------------------------------
mytheme <- theme(
  text=element_text(size=20), # change font size of all text
  #axis.text=element_text(size=20), # change font size of axis text
  #axis.title=element_text(size=20), # change font size of axis titles
  #plot.title=element_text(size=20), # change font size of plot title
  #legend.text=element_text(size=20), # change font size of legend text
  #legend.title=element_text(size=20) # change font size of legend title
)


#-------------------------------------------------------------------------------
# MODEL SELECTION
#-------------------------------------------------------------------------------
a <- plot_scores_nmf(x, types=get_types(x), remove_outliers=FALSE) + mytheme


#-------------------------------------------------------------------------------
# DENOVO SIGNATURES [SBS + DBS]
#-------------------------------------------------------------------------------

#dn_list <- c(2, 7, 8, 10, 11) # breast
#dn_list <- c(1, 2, 3, 5, 6, 7) # lung
#dn_list <- c(2, 3, 6, 7, 8, 9, 12, 13, 15, 16) # colorectal

denovo_1 <- basilica::plot_signatures(
  x, types = get_types(x), signames = "SBSD4"
  )[[1]] + 
  theme(
    axis.text.x = element_blank(), 
    axis.ticks.x = element_blank(), 
    axis.text.y = element_blank(), 
    axis.ticks.y = element_blank()
  ) + 
  labs(
    #x = "Context", 
    y = "Density"
  ) + mytheme

cosmic_1 <- plot.signatures(SBS = wide_to_long(COSMIC_sbs_filt["SBS98", ], what = "beta"), DBS = NULL, types = c("SBS"), context = T, cls = NULL) + 
  labs(
    x = "Context", 
    y = "Density"
  ) + 
  theme(
    axis.text.x = element_blank(), 
    axis.ticks.x = element_blank(), 
    axis.text.y = element_blank(), 
    axis.ticks.y = element_blank()
  ) + 
  mytheme

denovo_2 <- basilica::plot_signatures(
  x, types = get_types(x_mapped), signames = "SBSD8" 
)[[1]] + 
  theme(
    axis.text.x = element_blank(), 
    axis.ticks.x = element_blank(), 
    axis.text.y = element_blank(), 
    axis.ticks.y = element_blank()
  ) + 
  labs(
    #x = "Context", 
    #y = "Density"
  ) + mytheme

serena_2 <- plot.signatures(
  SBS = wide_to_long(data$sbs_obj$rare["SBS97", ], what = "beta"), DBS = NULL, types = c("SBS"), context = T, cls = NULL) + 
  labs(
    x = "Context", 
    #y = "Density"
  ) + 
  theme(
    axis.text.x = element_blank(), 
    axis.ticks.x = element_blank(), 
    axis.text.y = element_blank(), 
    axis.ticks.y = element_blank()
  ) + 
  mytheme


b <- patchwork::wrap_plots(denovo_1, denovo_2, cosmic_1, serena_2)


#-------------------------------------------------------------------------------
# HEATMAP
#-------------------------------------------------------------------------------
c <- plot_heatmp_comp(data$sbs_map) + mytheme


#-------------------------------------------------------------------------------
# CENTROIDS + AETIOLOGY
#-------------------------------------------------------------------------------

if (organ_type == "Breast") {
  sig_cls_organ <- list(
    G0 = c("SBS2", "SBS13", "DBS2", "DBS13", "DBS11"), 
    G1 = c("SBS3", "DBS2", "DBS13"), 
    G10 = c("SBS1", "SBS3", "SBS2", "SBS13", "SBSD11", "DBS11"), 
    G11 = c("SBS1", "SBS3", "SBS2", "SBS13", "SBSD11", "DBS14", "DBS13"), 
    G13 = c("SBS1", "SBS3", "SBS2", "SBS13", "SBSD11", "DBS13")
  )
} else if (organ_type == "Lung") {
  sig_cls_organ <- list(
    G0 = c("SBS1", "SBS2", "SBS3", "SBS5", "SBS13", "DBS6", "DBS13"), 
    G1 = c("SBS3", "DBS2", "DBS13"), 
    G10 = c("SBS1", "SBS3", "SBS2", "SBS13", "SBSD11", "DBS11"), 
    G11 = c("SBS1", "SBS3", "SBS2", "SBS13", "SBSD11", "DBS14", "DBS13"), 
    G13 = c("SBS1", "SBS3", "SBS2", "SBS13", "SBSD11", "DBS13")
  )
} else if (organ_type == "Colorectal") {
  sig_cls_organ <- list(
    G0 = c("SBS2", "SBS13", "DBS2", "DBS13", "DBS11"), 
    G1 = c("SBS3", "DBS2", "DBS13"), 
    G10 = c("SBS1", "SBS3", "SBS2", "SBS13", "SBSD11", "DBS11"), 
    G11 = c("SBS1", "SBS3", "SBS2", "SBS13", "SBSD11", "DBS14", "DBS13"), 
    G13 = c("SBS1", "SBS3", "SBS2", "SBS13", "SBSD11", "DBS13")
  )
} else {
  warning("wring tumor type!")
}

d <- custom_centroid_plot(
  x=x_mapped, 
  clusters=get_cluster_labels(x_mapped), 
  sig_cls=sig_cls_organ, 
  types=get_types(x_mapped), 
  sbs_aetiology_path=paste(path, "Nextcloud/basilica/sbs_aetiology.csv", sep = ""), 
  dbs_aetiology_path=paste(path, "Nextcloud/basilica/dbs_aetiology.csv", sep = ""), 
  custom_theme = mytheme
)


#d <- plot_centroids_etiology(
#  x = x_mapped, 
#  sbs_aetiology_path = paste(path, "Nextcloud/basilica/sbs_aetiology.csv", sep = ""), 
#  dbs_aetiology_path = paste(path, "Nextcloud/basilica/dbs_aetiology.csv", sep = ""), 
#  exposure_thr = 0.05, 
#  quantile_thr = 0.9, 
#  custom_theme = mytheme
#)


#-------------------------------------------------------------------------------
# CLUSTER SCORE - [one cluster + final score] - [SBS + DBS]
#-------------------------------------------------------------------------------

#plot_cluster_scores(x = x_mapped, type = "SBS", cluster_label = "G1", final_score = T, exposure_thr = 0.05, quantile_thr = 0.9)

#e <- plot_cluster_scores(x_mapped, "SBS", cluster_label = "G1", final_score = T, exposure_thr = 0, quantile_thr = 0.9) + 
#  mytheme

#f <- plot_cluster_scores(x_mapped, "DBS", cluster_label = "G1", final_score = T, exposure_thr = 0,  quantile_thr = 0.9) + 
#  mytheme


breast_G0 <- c("SBS2", "SBS13", "DBS11", "DBS13", "DBS2", "Others")
breast_G1 <- c("SBS3", "DBS13", "DBS2", "Others")

#lung_G0 <- c("SBS2", "SBS13", "DBS11", "DBS13", "DBS2")
#lung_G1 <- c("SBS3", "DBS13", "DBS2")

#colorectal_G0 <- c("SBS2", "SBS13", "DBS11", "DBS13", "DBS2")
#colorectal_G1 <- c("SBS3", "DBS13", "DBS2")

g <- plot_custom_exposures(x=x_mapped, clusters="G0", signature_list=breast_G0) + labs(x = "Samples", y = "Density", title = "Exposure Plot") + mytheme
f <- plot_custom_exposures(x=x_mapped, clusters="G1", signature_list=breast_G1) + labs(x = "Samples", y = "Density", title = "Exposure Plot") + mytheme


#-------------------------------------------------------------------------------
# FINAL PLOT
#-------------------------------------------------------------------------------

# a : model selection
# b : denovo signatures
# c : heatmap
# d : custom centroids
# g : exposure 1
# f : exposure 2


design <- "aaabbb
           ccccdd
           ccccdd
           gggggg
           ffffff"

p <- patchwork::wrap_plots(
  a, b, c, d, g, f, 
  design = design
) + patchwork::plot_annotation(tag_levels = "A") + theme(family = "serif")
# widths = 1, heights = 2)

ggsave(
  paste(path, "Nextcloud/basilica/results/", tolower(organ_type), ".main.png", sep = ""), 
  dpi=300, width = 23, height = 33, plot = p
)


#===============================================================================





