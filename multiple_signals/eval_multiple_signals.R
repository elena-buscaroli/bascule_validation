devtools::load_all("~/GitHub/simbasilica/")
devtools::load_all("~/GitHub/basilica/")
run_id = "clustering.matched.2011"
save_path = "~/Dropbox/shared/2022. Basilica/simulations/stats_dataframes/"

# Generate stats #####
source("~/GitHub/basilica_validation/eval_aux_fns.R")

runids = c("Autoguide", "ManualGuide")
fitnames = c("x.fit0.auto", "x.fit0.man")

path = paste0("~/Dropbox/shared/2022. Basilica/simulations/fits/fits_dn.", run_id, "/")
files = list.files(path, full.names=T, pattern=".Rds")

all_stats = lapply(files, function(fname) {
  stats_single_data(fname, names_fits=fitnames %>% setNames(runids))
}) %>% dplyr::bind_rows()
# saveRDS(all_stats, paste0(save_path, "stats_", run_id, ".Rds"))



# Plots #####
# all_stats = readRDS(paste0(save_path, "stats_", run_id, ".Rds"))
id_cols = c("N","G","seed","idd","fname","type","penalty")
list_cols = c("assigned_missing","input_sigs","fixed_sigs","dn_sigs")

plot_list = list()
## fixed/dn retrieved ####
plot_list[["K"]] = all_stats %>% 
  compute_quantiles(colname="K_true") %>% 
  dplyr::filter(penalty=="Autoguide") %>% 
  dplyr::rowwise() %>% 
  dplyr::mutate(K_input_found=length(intersect(assigned_missing$assigned_tp, input_sigs)),
                K_input_true=length(input_sigs),
                K_input_ratio=K_input_found / K_input_true,
                
                K_dn_found=length(intersect(assigned_missing$assigned_tp, dn_sigs)),
                K_dn_true=K_true - length(input_sigs),
                K_dn_ratio=K_dn_found / K_dn_true,
                
                K_found=length(assigned_missing$assigned_tp),
                K_ratio=K_found / K_true) %>% 
  dplyr::select(-dplyr::all_of(list_cols)) %>% 
  reshape2::melt(id=c(id_cols,"K_true_cat"), variable.name="metric") %>% 
  dplyr::filter(grepl("_ratio$", metric)) %>% 
  ggplot() +
  geom_violin(aes(x=factor(N), y=value, fill=K_true_cat, color=K_true_cat)) +
  # geom_jitter(aes(x=factor(N), y=value), height=0, width=0.2, size=.5) +
  geom_hline(yintercept=1, color="red4", linetype="dashed") +
  facet_grid(type ~ metric, scales="free_y") + 
  scale_fill_manual(values=c("tan2","#8FBC8B","thistle2")) +
  scale_color_manual(values=c("tan2","#8FBC8B","thistle2")) +
  theme_bw() + ylim(0,NA)




## quality metrics #####
metrics = c("mse_counts","cosine_expos","cosine_expos_missing","cosine_sigs")
cosine = all_stats %>% 
  dplyr::select(-dplyr::all_of(list_cols)) %>% 
  reshape2::melt(id=id_cols, variable.name="metric") %>% 
  dplyr::filter(metric %in% metrics[-1], penalty=="Autoguide") %>% 
  tidyr::separate("metric", into=c("metric","variable"), extra="merge", sep="_") %>% 
  ggplot() +
  geom_boxplot(aes(x=factor(N), y=filter_lims(value))) +
  # geom_jitter(aes(x=factor(N), y=filter_lims(value)), height=0, width=0.2, size=.5) +
  ggh4x::facet_nested(type ~ variable, scales="free_y") +
  theme_bw() + ylab("Cosine similarity")

mse = all_stats %>% 
  dplyr::select(-dplyr::all_of(list_cols)) %>% 
  reshape2::melt(id=id_cols, variable.name="metric") %>% 
  dplyr::filter(metric %in% metrics[1], penalty=="Autoguide") %>% 
  tidyr::separate("metric", into=c("metric","variable"), extra="merge", sep="_") %>% 
  ggplot() +
  geom_boxplot(aes(x=factor(N), y=filter_lims(value))) +
  # geom_jitter(aes(x=factor(N), y=filter_lims(value)), height=0, width=0.2, size=.5) +
  ggh4x::facet_nested(type ~ variable, scales="free_y") +
  theme_bw() + ylab("MSE")
plot_list[["performance"]] = patchwork::wrap_plots(mse, cosine, design="ABBB")

  

## clustering validation #####
plot_list[["clustering"]] = all_stats %>% 
  dplyr::select(-dplyr::all_of(list_cols)) %>%
  reshape2::melt(id=id_cols, variable.name="metric") %>% 
  dplyr::select(-type) %>% unique() %>% 
  dplyr::filter(metric %in% c("ari","nmi")) %>% 
  dplyr::rename(guide_t=penalty) %>% 
  ggplot() +
  geom_boxplot(aes(x=factor(N), y=filter_lims(value), fill=guide_t), na.rm=T) +
  geom_hline(yintercept=1, color="grey70", linetype="dashed") +
  ggh4x::facet_nested(~metric, scales="free_y") +
  scale_fill_manual(values=c("tan2","steelblue2")) +
  theme_bw() + theme(legend.position="bottom")


## save plots #####
patchwork::wrap_plots(plot_list$K, plot_list$performance, plot_list$clustering,
                      design="AAABBBBBCC\nAAABBBBB##")
ggsave(paste0(save_path, "stats_", run_id, ".pdf"), height=5, width=14)


# OLD plots #####
make_plots_stats_compare(all_stats, boxplot=T)
# ggsave(paste0("~/Dropbox/shared/2022. Basilica/simulations/plots_", run_id, "_boxp.pdf"),
#        height=8, width=12)

make_plots_stats_compare(all_stats, boxplot=F)
# ggsave(paste0("~/Dropbox/shared/2022. Basilica/simulations/plots_", run_id, "_line.pdf"),
#        height=8, width=12)




### save plots for each fit #####
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







