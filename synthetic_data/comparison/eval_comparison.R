devtools::load_all("~/GitHub/simbascule/")
devtools::load_all("~/GitHub/bascule/")
run_id = "matched.2011.compare"
save_path = "~/Dropbox/dropbox_shared/2022. Basilica/simulations/stats_dataframes/"
source("~/GitHub/bascule_validation/eval_aux_fns.R")
source("~/GitHub/bascule_validation/plots_aux_fns.R")

# Generate stats dataframe ##### 
# runids = c("BASCULE", "SigProfiler", "SparseSignatures")
# fitnames = c("fit.0", "sigprofiler", "sparsesignatures")
# 
# path = paste0("~/Dropbox/dropbox_shared/2022. Basilica/simulations/fits/fits_dn.", run_id, "/")
# files = list.files(path, full.names=T, pattern=".Rds")
# 
# all_stats = lapply(files, function(fname) {
#   stats_single_data(fname, names_fits=fitnames %>% setNames(runids))
# }) %>% dplyr::bind_rows()
# saveRDS(all_stats, paste0(save_path, "stats_", run_id, ".Rds"))



# Make plots #####
all_stats = readRDS(paste0(save_path, "stats_", run_id, ".Rds")) %>% 
  compute_quantiles(colname="K_true") %>% 
  dplyr::filter(type=="SBS")

pal = RColorBrewer::brewer.pal(3, name="Dark2")

plot_list = list()
plot_list[["K"]] = all_stats %>% 
  dplyr::select(-K_input_ratio, -K_dn_ratio) %>% 
  plot_K(fill="penalty", facet="~metric + K_true_cat", pal=pal)

plot_list[["performance"]] = all_stats %>% plot_performance(fill="penalty", pal=pal)
plot_list[["performance_grps"]] = all_stats %>% plot_performance(fill="penalty", facet="~variable + K_true_cat", pal=pal) &
  patchwork::plot_layout(ncol=1)


# Save plots #####
saveRDS(plot_list, paste0(save_path, "stats_", run_id, "_plots.Rds"))

patchwork::wrap_plots(plot_list, design="CCCCC\nAABBB") &
  theme(legend.position="bottom")
ggsave(paste0(save_path, "stats_", run_id, ".pdf"), width=10, height=10)




# save plots for each fit #####
lapply(files, function(fname) {
  tmp = strsplit(fname, split="/")[[1]]; fit_id = tmp[length(tmp)]
  simul_fit = readRDS(fname)
  x.simul = simul_fit$dataset
  types = get_types(x.simul)
  
  plots = lapply(fitnames, function(fitname) {
    x.fit = simul_fit[[fitname]] %>%
      # rename_dn_expos() %>%
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
  
  QC = lapply(fitnames, function(fitname) {
    x.fit = simul_fit[[fitname]] %>%
      merge_clusters()
    assigned_missing = get_assigned_missing(x.fit=x.fit, x.simul=x.simul)
    added = sapply(assigned_missing, function(i) i[["added_fp"]]) %>% unlist() %>% paste(collapse=",")
    missing = sapply(assigned_missing, function(i) i[["missing_fn"]]) %>% unlist() %>% paste(collapse=",")
    
    caption = paste("Added signatures (FP):", added, "- Missing signatures (FN):", missing)
    
    plot_QC(x.fit) & patchwork::plot_annotation(title=paste0(fit_id, " , fitname: ", fitname), caption=caption)
  })
  
  pdf(stringr::str_replace_all(fname, ".Rds", ".pdf") %>% stringr::str_replace_all("simul_fit","plots_fit"), width=20, height=16)
  print(plots)
  dev.off()
  
  pdf(stringr::str_replace_all(fname, ".Rds", ".pdf") %>% stringr::str_replace_all("simul_fit","plots_QC"), width=20, height=16)
  print(QC)
  dev.off()
  
})

