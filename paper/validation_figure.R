df_path = "~/Dropbox/dropbox_shared/2022. Basilica/simulations/stats_dataframes/"
stats_basilica = readRDS(paste0(df_path, "stats_clustering.matched.2011.Rds"))
stats_compare = readRDS(paste0(df_path, "stats_clustering.matched.2011.Rds"))


# panelA -> basilica validation: K_ratio by complexity, MSE counts, cosine sigs, NMI
