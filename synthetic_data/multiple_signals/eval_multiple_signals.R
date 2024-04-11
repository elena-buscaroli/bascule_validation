devtools::load_all("~/GitHub/simbasilica/")
devtools::load_all("~/GitHub/basilica/")
run_id = "matched.2011"
save_path = "~/Dropbox/dropbox_shared/2022. Basilica/simulations/stats_dataframes/"
source("~/GitHub/basilica_validation/aux_fns/eval_aux_fns.R")
source("~/GitHub/basilica_validation/aux_fns/plots_aux_fns.R")

# Generate stats #####
# runids = c("Autoguide", "ManualGuide")
# fitnames = c("fit.0.auto", "fit.0.man")
# 
# path = paste0("~/Dropbox/dropbox_shared/2022. Basilica/simulations/fits/fits_dn.", run_id, "/")
# files = list.files(path, full.names=T, pattern=".Rds")
# 
# all_stats = lapply(files, function(fname) {
#   stats_single_data(fname, names_fits=fitnames %>% setNames(runids))
# }) %>% dplyr::bind_rows()
# saveRDS(all_stats, paste0(save_path, "stats_", run_id, "_KM.Rds"))



# Plots #####
all_stats = readRDS(paste0(save_path, "stats_", run_id, "_KM.Rds")) %>% 
  compute_quantiles(colname="K_true")

plot_list = list("p1"=list(), "p2"=list())
## fixed/dn retrieved ####
plot_list$p1[["K"]] = all_stats %>% 
  compute_quantiles(colname="K_true") %>% 
  dplyr::filter(penalty=="Autoguide") %>% 
  plot_K()

plot_list$p2[["K"]] = all_stats %>% 
  compute_quantiles(colname="K_true") %>% 
  dplyr::filter(penalty=="Autoguide") %>% 
  plot_K(fill="K_true_cat")


## quality metrics #####
plot_list$p1[["performance"]] = all_stats %>% 
  dplyr::filter(penalty=="Autoguide") %>% 
  plot_performance()

plot_list$p2[["performance"]] = all_stats %>%
  dplyr::filter(penalty=="Autoguide") %>%
  plot_performance(fill="K_true_cat")



## clustering validation #####
plot_list$p1[["clustering"]] = all_stats %>% plot_performance_clustering()
plot_list$p2[["clustering"]] = all_stats %>% 
  plot_performance_clustering(fill="K_true_cat", facet="penalty~metric")


# Save plots #####
saveRDS(plot_list, paste0(save_path, "stats_",run_id, "_plots.Rds"))

patchwork::wrap_plots(plot_list$p1, design="BBBBB\nAAACC")
ggsave(paste0(save_path, "stats_", run_id, "_nogrps.pdf"), width=10, height=10)

patchwork::wrap_plots(plot_list$p2, design="BBBBB\nAAACC") &
  theme(legend.position="bottom")
ggsave(paste0(save_path, "stats_", run_id, "_grps.pdf"), width=10, height=10)



### save plots for each fit #####
fitnames = c("fit.0.auto", "fit.0.man")
path = paste0("~/Dropbox/dropbox_shared/2022. Basilica/simulations/fits/fits_dn.", run_id, "/")
files = list.files(path, full.names=T, pattern=".Rds")

lapply(files, function(fname) {
  tmp = strsplit(fname, split="/")[[1]]; fit_id = tmp[length(tmp)]
  simul_fit = readRDS(fname)
  x.simul = simul_fit$dataset
  types = get_types(x.simul)

  plots = lapply(fitnames, function(fitname) {
    x.fit = simul_fit[[fitname]] %>%
      convert_dn_names(x.simul=x.simul) %>% 
      merge_clusters()
    assigned_missing = get_assigned_missing(x.fit=x.fit, x.simul=x.simul)
    added = sapply(assigned_missing, function(i) i[["added_fp"]]) %>% unlist() %>% paste(collapse=",")
    missing = sapply(assigned_missing, function(i) i[["missing_fn"]]) %>% unlist() %>% paste(collapse=",")

    caption = paste("Added signatures (FP):", added, "- Missing signatures (FN):", missing)

    patchwork::wrap_plots(plot_fit(x.fit),
                          plot_fit(x.simul), ncol=1) &
      patchwork::plot_annotation(title=paste0(fit_id, " , fitname: ", fitname),
                                 subtitle="Fit (top) and simulated (bottom)",
                                 caption=caption)
  })

  # QC = lapply(fitnames, function(fitname) {
  #   x.fit = simul_fit[[fitname]] %>%
  #     merge_clusters()
  #   assigned_missing = get_assigned_missing(x.fit=x.fit, x.simul=x.simul)
  #   added = sapply(assigned_missing, function(i) i[["added_fp"]]) %>% unlist() %>% paste(collapse=",")
  #   missing = sapply(assigned_missing, function(i) i[["missing_fn"]]) %>% unlist() %>% paste(collapse=",")
  # 
  #   caption = paste("Added signatures (FP):", added, "- Missing signatures (FN):", missing)
  # 
  #   plot_QC(x.fit) & patchwork::plot_annotation(title=paste0(fit_id, " , fitname: ", fitname), caption=caption)
  # })

  pname = stringr::str_replace_all(fname, ".Rds", ".pdf") %>% stringr::str_replace_all("simul_fit","plots_fit")
  pdf(pname, width=20, height=20)
  print(plots)
  dev.off()

  # pdf(stringr::str_replace_all(fname, ".Rds", ".pdf") %>% stringr::str_replace_all("simul_fit","plots_QC"), width=20, height=16)
  # print(QC)
  # dev.off()

})







