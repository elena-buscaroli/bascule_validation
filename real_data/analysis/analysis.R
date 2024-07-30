
organ_type <- "Breast"
data <- readRDS( paste0("/Users/azadsadr/Nextcloud/basilica/compiled_fit/", tolower(organ_type), "_data.Rds") )
path <- "/Users/azadsadr/"


devtools::load_all( paste(path, "Documents/packages/basilica/", sep = "") )
source(paste(path, "Nextcloud/basilica/scripts/utils_plot.R", sep = ""))

x_before <- data$x_before
x_after <- data$x_after

x_merged <- basilica:::merge_clusters(x_after, cutoff = 0.8)
x_mapped <- basilica:::convert_dn_names(x_merged, reference_cat = list(SBS=COSMIC_sbs_filt, DBS=COSMIC_dbs), cutoff = 0.8)


#===============================================================================
#                                      SBS                                     #
#===============================================================================
context_type <- "SBS"


#                          REFINEMENT ANALYSIS  / SBS                          
#-------------------------------------------------------------------------------
p00 <- basilica::plot_signatures(x_before, types = c(context_type) )
ggsave( paste(path, "Nextcloud/basilica/results/", tolower(context_type), ".", tolower(organ_type), ".sig.before.png", sep = "") , dpi=300, width = 16, height = 13, plot = p00)

p01 <- basilica::plot_signatures(x_after, types = c(context_type) )
ggsave( paste(path, "Nextcloud/basilica/results/", tolower(context_type), ".", tolower(organ_type), ".sig.png", sep = "") , dpi=300, width = 16, height = 12, plot = p01)


#                              NMF ANALYSIS / SBS
#-------------------------------------------------------------------------------
p10 <- basilica::plot_signatures(x_after, types = c(context_type), signames = get_denovo_signames(x_after))
ggsave( paste(path, "Nextcloud/basilica/results/", tolower(context_type), ".", tolower(organ_type), ".denovo.png", sep = "") , dpi=300, width = 12, height = 8, plot = p10)

p20 <- plot_unmapped(
  basilica_exp = data$sbs_obj$BasilicaExposure, 
  serena_exp = data$sbs_obj$SerenaExposure, 
  basilica_singles = data$sbs_map$basilica_singles, 
  serena_singles = data$sbs_map$serena_singles, 
  exposure_thr=0.01
)
ggsave( paste(path, "Nextcloud/basilica/results/", tolower(context_type), ".", tolower(organ_type), ".unmapped.sig.png", sep = "") , dpi=300, width = 10, height = 6, plot = p20)

p30 <- plot_exposure_comp(data$sbs_obj, data$sbs_map)
ggsave( paste(path, "Nextcloud/basilica/results/", tolower(context_type), ".", tolower(organ_type), ".exp.alignment.png", sep = "") , dpi=300, width = 16, height = 9, plot = p30)

p40 <- plot_heatmp_comp(data$sbs_map)
ggsave( paste(path, "Nextcloud/basilica/results/", tolower(context_type), ".", tolower(organ_type), ".sig.heatmap.png", sep = "") , dpi=300, width = 10, height = 8, plot = p40)



#===============================================================================
#                                      DBS                                     #
#===============================================================================
context_type <- "DBS"


#                          REFINEMENT ANALYSIS  / DBS                          
#-------------------------------------------------------------------------------
t00 <- basilica::plot_signatures(x_before, types = c(context_type) )
ggsave( paste(path, "Nextcloud/basilica/results/", tolower(context_type), ".", tolower(organ_type), ".sig.before.png", sep = "") , dpi=300, width = 12, height = 8, plot = t00)

t01 <- basilica::plot_signatures(x_after, types = c(context_type) )
ggsave( paste(path, "Nextcloud/basilica/results/", tolower(context_type), ".", tolower(organ_type), ".sig.png", sep = "") , dpi=300, width = 12, height = 8, plot = t01)


#                              NMF ANALYSIS / DBS
#-------------------------------------------------------------------------------
t10 <- basilica::plot_signatures(x_after, types = c(context_type), signames = get_denovo_signames(x_after))
ggsave( paste(path, "Nextcloud/basilica/results/", tolower(context_type), ".", tolower(organ_type), ".denovo.png", sep = "") , dpi=300, width = 10, height = 6, plot = t10)


t20 <- plot_unmapped(
  basilica_exp = data$dbs_obj$BasilicaExposure, 
  serena_exp = data$dbs_obj$SerenaExposure, 
  basilica_singles = data$dbs_map$basilica_singles, 
  serena_singles = data$dbs_map$serena_singles, 
  exposure_thr=0.01
)
ggsave( paste(path, "Nextcloud/basilica/results/", tolower(context_type), ".", tolower(organ_type), ".unmapped.sig.png", sep = "") , dpi=300, width = 14, height = 10, plot = t20)


t30 <- plot_exposure_comp(
  data$dbs_obj, data$dbs_map
  #basilica_exp = data$dbs_obj$BasilicaExposure, 
  #serena_exp = data$dbs_obj$SerenaExposure, 
  #basilica_serena = data$dbs_map$basilica_serena
)
ggsave( paste(path, "Nextcloud/basilica/results/", tolower(context_type), ".", tolower(organ_type), ".exp.alignment.png", sep = "") , dpi=300, width = 10, height = 6, plot = t30)


t40 <- plot_heatmp_comp(data$dbs_map)
ggsave( paste(path, "Nextcloud/basilica/results/", tolower(context_type), ".", tolower(organ_type), ".sig.heatmap.png", sep = "") , dpi=300, width = 10, height = 8, plot = t40)



#===============================================================================
#                             CLUSTERING ANALYSIS
#===============================================================================

c10 <- plot_mixture_weights (x_merged, empirical = T)
#c10 <- plot_cluster_freq(x_merged)
ggsave( paste(path, "Nextcloud/basilica/results/cls/", tolower(organ_type), ".freq.png", sep = "") , dpi=300, width = 8, height = 5, plot = c10)




for ( cls in basilica:::get_cluster_labels(x_merged) ) {
  cc <- patchwork::wrap_plots(
    basilica:::plot_cluster_scores(x = x_mapped, type = "SBS", cluster_label = cls, final_score = T), 
    basilica:::plot_cluster_scores(x = x_mapped, type = "DBS", cluster_label = cls, final_score = T), 
    basilica:::plot_cluster_scores(x = x_mapped, type = "SBS", cluster_label = cls, final_score = F), 
    basilica:::plot_cluster_scores(x = x_mapped, type = "DBS", cluster_label = cls, final_score = F), 
    ncol = 2
  )
  ggsave( 
    paste(path, "Nextcloud/basilica/results/cls/", tolower(organ_type), ".", cls, ".scores.png", sep = ""), 
    dpi=300, width = 16, height = 8, plot = cc
  )
}



c30 <- plot_cls_score_heatmap(x_mapped, type="SBS", exposure_thr = 0.05)
ggsave( paste(path, "Nextcloud/basilica/results/cls/sbs.", tolower(organ_type), ".cls.heatmap", ".png", sep = "") , dpi=300, width = 8, height = 8, plot = c30)

c30 <- plot_cls_score_heatmap(x_mapped, type="DBS", exposure_thr = 0.05)
ggsave( paste(path, "Nextcloud/basilica/results/cls/dbs.", tolower(organ_type), ".cls.heatmap", ".png", sep = "") , dpi=300, width = 8, height = 8, plot = c30)



# test for different context types
c40 <- plot_centroids(x_mapped, types = get_types(x_merged), cls = NULL, sort_by = NULL, exposure_thr = 0, quantile_thr = 0)
ggsave( paste(path, "Nextcloud/basilica/results/cls/", tolower(organ_type), ".centroids.png", sep = "") , dpi=300, width = 10, height = 5, plot = c40)

c41 <- plot_centroids_etiology(
  x_mapped, 
  sbs_aetiology_path = paste(path, "Nextcloud/basilica/sbs_aetiology.csv", sep = ""), 
  dbs_aetiology_path = paste(path, "Nextcloud/basilica/dbs_aetiology.csv", sep = ""), 
  custom_theme = theme(text=element_text(size=20))
)
ggsave( paste(path, "Nextcloud/basilica/results/cls/", tolower(organ_type), ".centroids.etiology.png", sep = "") , dpi=300, width = 10, height = 8, plot = c41)



get_cluster_labels(x_mapped)
ppp <- plot_exposures(x_mapped, clusters = "G13")
ggsave( paste(path, "Nextcloud/basilica/results/breast/cls/", tolower(organ_type), ".exposure.G13.png", sep = "") , dpi=300, width = 16, height = 8, plot = ppp)






