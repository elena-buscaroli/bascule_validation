library(magrittr)
library(ggplot2)
source("~/GitHub/bascule_validation/paper/figure2/eval_aux_fns.R")
source("~/GitHub/bascule_validation/paper/figure2/plots_aux_fns.R")

df_path = "~/Dropbox/dropbox_shared/2022. Basilica/simulations/stats_dataframes/"
stats_bascule = readRDS(paste0(df_path, "stats_clustering.matched.2011.Rds")) %>% 
  compute_quantiles(colname="K_true")
stats_compare = readRDS(paste0(df_path, "stats_matched.2011.compare.Rds")) %>% 
  compute_quantiles(colname="K_true") %>% 
  dplyr::filter(type=="SBS")

## Basilica ####

pal_k_true = pal=c("tan2","#8FBC8B","thistle2")

# panelA -> bascule validation: K_ratio by complexity, MSE counts, cosine sigs, NMI
K_ratio = stats_bascule %>% 
  dplyr::filter(type == "SBS") %>% 
  compute_quantiles(colname="K_true") %>% 
  dplyr::select(-K_input_ratio, -K_dn_ratio) %>% 
  dplyr::filter(penalty=="Autoguide") %>% 
  plot_K(fill="K_true_cat", pal=pal_k_true)

## accuracy : (TP + TN) / (TP + TN + FP + FN)
## precision : TP / (TP + FP)
## recall : TP / (TP + FN)

mse_counts = stats_bascule %>% 
  dplyr::select(-mse_expos, -mse_expos_missing, -dplyr::contains("cosine")) %>% 
  dplyr::filter(penalty=="Autoguide") %>% 
  plot_performance(fill="K_true_cat", pal=pal_k_true)

cosine_sigs_expos = stats_bascule %>% 
  dplyr::select(-dplyr::contains("mse"), -cosine_expos_missing) %>% 
  dplyr::filter(penalty=="Autoguide") %>% 
  plot_performance(fill="K_true_cat", pal=pal_k_true)

clustering = stats_bascule %>% 
  dplyr::filter(penalty=="Autoguide", type=="SBS") %>% 
  dplyr::select(-ari) %>% 
  
  plot_performance_clustering(fill="K_true_cat", facet="~metric", pal=pal_k_true)


## Comparison #####

# panelB -> bascule comparison w other methods
pal_methods = RColorBrewer::brewer.pal(3, name="Dark2")

K_ratio_cmp = stats_compare %>% 
  dplyr::select(-K_input_ratio, -K_dn_ratio) %>% 
  plot_K(fill="penalty", facet="~metric", pal=pal_methods)

K_ratio_cmp_fct = stats_compare %>% 
  dplyr::select(-K_input_ratio, -K_dn_ratio) %>% 
  plot_K(fill="penalty", facet="K_true_cat~metric", pal=pal_methods)

precision_recall = stats_compare %>% dplyr::rowwise() %>%
  dplyr::mutate(FN=length(assigned_missing$missing_fn),
                FP=length(assigned_missing$added_fp),
                TP=length(assigned_missing$assigned_tp)) %>%
  dplyr::mutate(precision=TP / (TP + FP),
                recall=TP / (TP + FN)) %>% 
  tidyr::pivot_longer(cols=c("precision","recall"), names_to="false_pos_neg") %>% 
  # dplyr::mutate(value=value/K_found,
  #               false_pos_neg=paste0(false_pos_neg, "_ratio")) %>% 
  ggplot() + 
  geom_violin(aes(x=factor(N), y=value, fill=factor(penalty)), 
              lwd=.3, draw_quantiles=c(.5),
              position=position_dodge(width=.7)) +
  theme_bw() + facet_grid(false_pos_neg~., scales="free") +
  scale_fill_manual(values=pal_methods)

# ggsave(plot=FP_FN_ratio, filename=paste0(df_path, "FN_FP.png"))
# ggsave(plot=K_ratio_cmp_fct, filename=paste0(df_path, "K_ratio_fct.png"))


mse_counts_cmp = stats_compare %>%
  dplyr::select(-mse_expos, -mse_expos_missing, -dplyr::contains("cosine")) %>% 
  plot_performance(fill="penalty", facet="~variable", pal=pal_methods)

cosine_expos_sigs_cmp = stats_compare %>%
  dplyr::select(-dplyr::contains("mse"), -cosine_expos) %>% 
  plot_performance(fill="penalty", facet="~variable", pal=pal_methods)



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
  dplyr::mutate(tool="SparseSignatures") %>% tibble::as_tibble() %>% 
  dplyr::select(simulation_name, execution_time, tool)
times_bascule = read.csv("~/Dropbox/dropbox_shared/2022. Basilica/simulations/runtimes/bascule_exectimes.csv") %>% 
  dplyr::rename(execution_time=execution_time_SBS) %>% 
  dplyr::mutate(tool="Bascule") %>% tibble::as_tibble() %>% 
  dplyr::select(simulation_name, execution_time, tool)

runtimes = dplyr::bind_rows(times_sigpr, times_sparsesig, times_bascule) %>% 
  dplyr::rowwise() %>% 
  dplyr::mutate(processor=dplyr::case_when(
    tool=="Bascule" && grep("N1000", simulation_name) ~ "GPU",
    .default="CPU"
  )) %>% 
  dplyr::mutate(N=strsplit(simulation_name, "[.]")[[1]][2] %>% stringr::str_remove_all("N") %>% as.numeric()) %>% 
  dplyr::select(simulation_name, execution_time, tool, processor, N) %>% 
  
  dplyr::group_by(tool, N) %>% 
  dplyr::filter(execution_time < boxplot.stats(execution_time)$stats[5]) %>% 
  
  ggplot() +
  geom_boxplot(aes(y=execution_time, x=factor(N), fill=tool)) +
  # facet_wrap(~processor, scales="free") +
  scale_fill_manual(values=pal_methods) +
  theme_bw()

# ggsave("~/Dropbox/dropbox_shared/2022. Basilica/simulations/stats_dataframes/runtimes.png")




# Panels #####

panelA = K_ratio + xlab("# samples") + ylab("Recall") + 
  labs(tag="A") +
  guides(fill=guide_legend(title="# signatures")) + 
  theme(legend.position="bottom")

panelB = patchwork::wrap_plots(mse_counts + labs(tag="B"), 
                               cosine_sigs_expos, guides="collect") &
  xlab("# samples") & theme(legend.position="bottom") &
  guides(fill=guide_legend(title="# signatures"))

panelC = clustering + xlab("# samples") + ylab("NMI") +
  labs(tag="C") +
  guides(fill=guide_legend(title="# signatures"))

panelD = precision_recall + labs(tag="D") + xlab("# samples") + ylab("") +
  guides(fill=guide_legend(title="Method")) +
  theme(legend.position="bottom")

panelE = patchwork::wrap_plots(mse_counts_cmp + labs(tag="E"), 
                               cosine_expos_sigs_cmp, guides="collect") &
  xlab("# samples") & theme(legend.position="bottom") &
  guides(fill=guide_legend(title="Method"))

panelF = runtimes + ylab("Runtime (min)") + xlab("# samples") +
  labs(tag="D") +
  guides(fill=guide_legend(title="Method")) +
  theme(legend.position="bottom")


# panelA = patchwork::wrap_plots(K_ratio + labs(tag="A"),
#                                patchwork::wrap_plots(mse_counts, cosine_sigs_expos,
#                                                      guides="collect"),
#                                design="ABBB") &
#   theme(legend.position="bottom")
# 
# panelB = patchwork::wrap_plots(K_ratio_cmp + labs(tag="B"), 
#                                patchwork::wrap_plots(mse_counts_cmp, cosine_expos_sigs_cmp, 
#                                                      guides="collect"),
#                                design="ABBB") &
#   theme(legend.position="bottom")
# 
# 
# panelC = clustering + labs(tag="C")


patchwork::wrap_plots(panelA, panelB, panelC, panelD, panelE, panelF,
                      design="ABBC\nDEEF")
ggsave(paste0(df_path, "draft_fig.pdf"), height=8, width=12)
ggsave(paste0(df_path, "draft_fig.png"), height=8, width=12)




