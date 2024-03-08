library(magrittr)
library(ggplot2)
source("~/GitHub/basilica_validation/paper/figure2/eval_aux_fns.R")
source("~/GitHub/basilica_validation/paper/figure2/plots_aux_fns.R")

df_path = "~/Dropbox/dropbox_shared/2022. Basilica/simulations/stats_dataframes/"
stats_basilica = readRDS(paste0(df_path, "stats_clustering.matched.2011.Rds")) %>% 
  compute_quantiles(colname="K_true")
stats_compare = readRDS(paste0(df_path, "stats_matched.2011.compare.Rds")) %>% 
  compute_quantiles(colname="K_true") %>% 
  dplyr::filter(type=="SBS")

plots = list()

# Basilica ####
pal_k_true = wesanderson::wes_palette("Cavalcanti1", 6, type="continuous")[c(1,4,3)]

stats_basilica %>% 
  dplyr::filter(type=="SBS") %>% 
  compute_quantiles(colname="K_true") %>% 
  dplyr::select(idd, N, K_true_cat, nmi, nmi_KM) %>% 
  
  tidyr::pivot_longer(cols=c("nmi","nmi_KM"), values_to="nmi", names_to="Method") %>% 
  dplyr::mutate(Method=dplyr::case_when(Method=="nmi" ~ "Basilica",
                                        Method=="nmi_KM" ~ "KMeans")) %>% 
  
  dplyr::group_by(N, K_true_cat) %>%
  dplyr::filter(! (nmi %in% boxplot.stats(nmi)$out) ) %>%
  
  ggplot() +
  geom_violin(aes(x=factor(N), y=nmi, fill=Method, color=Method), 
              alpha=.4, lwd=.3, position=position_dodge(width=.6)) +
  
  # geom_jitter(aes(x=factor(N), y=nmi_KM), # fill=K_true_cat, color=K_true_cat),
  #             alpha=.4, size=.3, height=0) + #, position=position_jitterdodge(dodge.width=.6, jitter.width=.2)) +
  
  facet_grid(~K_true_cat) +
  
  scale_color_manual(values=pal_k_true) +
  scale_fill_manual(values=pal_k_true) +
  theme_bw()
  


# Comparison ####
pal_methods = c("#7fb3d5", "#FF8C00", "#8FBC8B") %>% setNames(c("Basilica", "SigProfiler", "SparseSignatures"))

## precision and recall ####
plots[["recall"]] = stats_compare %>% dplyr::rowwise() %>%
  dplyr::mutate(FN=length(assigned_missing$missing_fn),
                FP=length(assigned_missing$added_fp),
                TP=length(assigned_missing$assigned_tp)) %>%
  dplyr::mutate(recall=TP / (TP + FN)) %>% 
  tidyr::pivot_longer(cols=c("recall"), names_to="prec_recall") %>% 
  dplyr::mutate(prec_recall=stringr::str_to_title(prec_recall)) %>% 
  
  dplyr::group_by(prec_recall, penalty, N) %>% 
  dplyr::filter(value > boxplot.stats(value)$stats[1]) %>% 
  
  ggplot(aes(x=factor(N), y=value, fill=factor(penalty), color=factor(penalty))) + 
  # geom_violin(alpha=.4, lwd=.3, 
  #             # color="#FF000000", 
  #             position=position_dodge(width=.6)) +
  stat_summary(aes(group=factor(penalty)), position=position_dodge(width=.1), fun.data="mean_cl_boot") +
  stat_summary(aes(group=factor(penalty)), fun.data="mean_cl_boot", position=position_dodge(width=.1), geom="line") +
  
  theme_bw() + 
  scale_y_continuous(n.breaks=5) +
  scale_fill_manual(values=pal_methods) +
  scale_color_manual(values=pal_methods)

## mse ####
plots[["mse_counts"]] = stats_compare %>%
  dplyr::select(N, idd, mse_counts, penalty) %>% dplyr::rename(value=mse_counts) %>% 
  dplyr::mutate(metric="Counts") %>% 
  
  dplyr::group_by(N, penalty) %>% 
  
  dplyr::filter(value < boxplot.stats(value)$stats[5]) %>% 
  dplyr::mutate(mmean=mean(value)) %>% 

  ggplot(aes(x=factor(N), y=value, fill=factor(penalty), color=factor(penalty))) + 
  # geom_violin(alpha=.4, lwd=.3, 
  #             # color="#FF000000", 
  #             position=position_dodge(width=.6)) +
  stat_summary(aes(group=factor(penalty)), position=position_dodge(width=.1), fun.data="mean_cl_boot") +
  stat_summary(aes(group=factor(penalty)), position=position_dodge2(width=.1), fun.data="mean_cl_boot", geom="line") +
  theme_bw() + facet_grid(~metric) +
  scale_y_continuous(n.breaks=5, labels=function(x) round(x, digits=3)) +
  scale_fill_manual(values=pal_methods) +
  scale_color_manual(values=pal_methods)


plots[["cosine_sigs"]] = stats_compare %>%
  dplyr::select(N, idd, cosine_sigs, penalty) %>% dplyr::rename(value=cosine_sigs) %>% 
  dplyr::mutate(metric="Signatures") %>% 
  dplyr::group_by(N, penalty) %>% 
  dplyr::filter(value > boxplot.stats(value)$stats[1]) %>% 
  dplyr::mutate(mmean=mean(value)) %>% 
  
  ggplot(aes(x=factor(N), y=value, fill=factor(penalty))) + 
  # geom_violin(alpha=.4, lwd=.3, 
  #             # color="#FF000000", 
  #             position=position_dodge(width=.6)) +
  stat_summary(aes(group=factor(penalty), color=factor(penalty)), position=position_dodge(width=.1),
               geom="pointrange", fun.data="mean_cl_boot") +
  stat_summary(aes(group=factor(penalty), color=factor(penalty)), geom="line",
               position=position_dodge(width=.1), fun.data="mean_cl_boot") +
  theme_bw() + facet_grid(~metric) +
  scale_y_continuous(n.breaks=5) +
  scale_fill_manual(values=pal_methods) +
  scale_color_manual(values=pal_methods)


plots[["cosine_expos"]] = stats_compare %>%
  dplyr::select(N, idd, cosine_expos_missing, penalty) %>% 
  dplyr::rename(value=cosine_expos_missing) %>% 
  dplyr::mutate(metric="Exposures") %>% 
  
  dplyr::group_by(N, penalty) %>% 
  dplyr::filter(value > boxplot.stats(value)$stats[1]) %>% 
  dplyr::mutate(mmean=mean(value)) %>% 
  
  ggplot(aes(x=factor(N), y=value, fill=factor(penalty), color=factor(penalty))) + 
  # geom_violin(alpha=.4, lwd=.3, 
  #             # color="#FF000000", 
  #             position=position_dodge(width=.6)) +
  
  stat_summary(aes(group=factor(penalty)), fun.data="mean_cl_boot", position=position_dodge(width=.1)) +
  stat_summary(aes(group=factor(penalty)), fun.data="mean_cl_boot", position=position_dodge(width=.1), geom="line") +
  
  theme_bw() + facet_grid(~metric) +
  scale_y_continuous(n.breaks=5) +
  scale_fill_manual(values=pal_methods) +
  scale_color_manual(values=pal_methods)


## runtimes ####
times_sigpr = read.csv("~/Dropbox/dropbox_shared/2022. Basilica/simulations/runtimes/last/sigprofiler_exectimes.csv") %>% 
  dplyr::mutate(tool="SigProfiler") %>% tibble::as_tibble() %>% 
  dplyr::mutate(execution_time=stringr::str_replace_all(execution_time, " 0:", "00:")) %>% 
  dplyr::mutate(execution_time=stringr::str_remove_all(execution_time, " ")) %>% 
  dplyr::mutate(execution_time=lubridate::period_to_seconds(lubridate::hms(execution_time))) %>% 
  dplyr::rename(simulation_name=simulation) %>% 
  dplyr::select(simulation_name, execution_time, tool)
times_sparsesig = read.csv("~/Dropbox/dropbox_shared/2022. Basilica/simulations/runtimes/last/sparsesignatures_exectimes.csv") %>% 
  dplyr::rename(execution_time=total_mins, simulation_name=name) %>% 
  dplyr::mutate(tool="SparseSignatures",
                simulation_name=stringr::str_remove_all(simulation_name, ".Rds")) %>% 
  tibble::as_tibble() %>% 
  dplyr::select(simulation_name, execution_time, tool)
times_basilica = read.csv("~/Dropbox/dropbox_shared/2022. Basilica/simulations/runtimes/basilica_exectimes.csv") %>% 
  dplyr::rename(execution_time_basilica=execution_time_SBS) %>% tibble::as_tibble() %>% 
  dplyr::mutate(tool="Basilica", execution_time=execution_time_basilica) %>%
  dplyr::select(simulation_name, execution_time, execution_time_basilica, tool)

plots[["runtime"]] = dplyr::bind_rows(times_sigpr, 
                                      times_sparsesig, 
                                      times_basilica %>% dplyr::select(-execution_time_basilica)) %>% 
  dplyr::inner_join(times_basilica %>% dplyr::select(-execution_time, -tool)) %>% 

  dplyr::mutate(time_gain=execution_time / execution_time_basilica) %>% 
  
  dplyr::rowwise() %>% 
  dplyr::mutate(N=strsplit(simulation_name, "[.]")[[1]][2] %>% stringr::str_remove_all("N") %>% as.numeric()) %>% 
  
  ggplot(aes(x=factor(N), y=log(time_gain, base=10), fill=factor(tool), color=factor(tool))) +
  stat_summary(aes(group=factor(tool)), position=position_dodge(width=.1), fun.data="mean_cl_boot") +
  stat_summary(aes(group=factor(tool)), fun.data="mean_cl_boot", position=position_dodge(width=.1), geom="line") +
  theme_bw() +
  scale_fill_manual(values=pal_methods, breaks = names(pal_methods), drop=F) +
  scale_color_manual(values=pal_methods, breaks = names(pal_methods), drop=F) + 
  scale_y_continuous(n.breaks=5, limits=c(1, NA))



# Panels #####
panelD = plots[["recall"]] + ylab("Recall") +
  labs(title="Signatures detection accuracy") +
  theme(legend.position="bottom")

panelE = plots[["mse_counts"]] + 
  labs(title="Reconstruction error") +
  ylab("MSE")

panelF = plots[["cosine_sigs"]] + 
  labs(title="Signatures quality") +
  ylab("Cosine similarity")

panelG = plots[["cosine_expos"]] +
  labs(title="Exposures quality") +
  ylab("Cosine similarity")

panelH = plots[["runtime"]] + 
  labs(title="Runtime increment") +
  ylab("Relative time increment (log10)")


comparison_plots = ((patchwork::wrap_plots(panelD, panelE, panelF, 
                        panelG, panelH,
                        guides="collect",
                        design="AABBCCDDEE") & 
    xlab("# samples")) &
    theme(legend.position="bottom") &
    guides(fill=guide_legend(title="Method"),
           color=guide_legend(title="Method")))

  # patchwork::plot_annotation(tag_levels="A")



