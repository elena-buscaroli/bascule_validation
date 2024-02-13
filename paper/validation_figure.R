library(magrittr)
library(ggplot2)
source("~/GitHub/basilica_validation/eval_aux_fns.R")
source("~/GitHub/basilica_validation/plots_aux_fns.R")

df_path = "~/Dropbox/dropbox_shared/2022. Basilica/simulations/stats_dataframes/"
stats_basilica = readRDS(paste0(df_path, "stats_clustering.matched.2011.Rds")) %>% 
  compute_quantiles(colname="K_true")
stats_compare = readRDS(paste0(df_path, "stats_matched.2011.compare.Rds")) %>% 
  compute_quantiles(colname="K_true") %>% 
  dplyr::filter(type=="SBS")


# panelA -> basilica validation: K_ratio by complexity, MSE counts, cosine sigs, NMI
K_ratio = stats_basilica %>% 
  compute_quantiles(colname="K_true") %>% 
  dplyr::select(-K_input_ratio, -K_dn_ratio) %>% 
  dplyr::filter(penalty=="Autoguide") %>% 
  plot_K(fill="K_true_cat")

mse_counts = stats_basilica %>% 
  dplyr::select(-mse_expos, -mse_expos_missing, -dplyr::contains("cosine")) %>% 
  dplyr::filter(penalty=="Autoguide") %>% 
  plot_performance(fill="K_true_cat")

cosine_sigs_expos = stats_basilica %>% 
  dplyr::select(-dplyr::contains("mse"), -cosine_expos_missing) %>% 
  dplyr::filter(penalty=="Autoguide") %>% 
  plot_performance(fill="K_true_cat")

clustering = stats_basilica %>% 
  dplyr::filter(penalty=="Autoguide", type=="SBS") %>% 
  dplyr::select(-ari) %>% 
  plot_performance_clustering(fill="K_true_cat", facet="~metric")


# panelB -> basilica comparison w other methods
pal = RColorBrewer::brewer.pal(3, name="Dark2")

K_ratio_cmp = stats_compare %>% 
  dplyr::select(-K_input_ratio, -K_dn_ratio) %>% 
  plot_K(fill="penalty", facet="~metric", pal=pal)

K_ratio_cmp_fct = stats_compare %>% 
  dplyr::select(-K_input_ratio, -K_dn_ratio) %>% 
  plot_K(fill="penalty", facet="K_true_cat~metric", pal=pal)

FP_FN_ratio = stats_compare %>% dplyr::rowwise() %>%
  dplyr::mutate(FN=length(assigned_missing$missing_fn),
                FP=length(assigned_missing$added_fp)) %>%
  tidyr::pivot_longer(cols=c("FN","FP"), names_to="false_pos_neg") %>% 
  dplyr::mutate(value=value/K_found) %>% 
  ggplot() + 
  geom_violin(aes(x=K_true_cat, y=value, fill=factor(penalty)), 
              lwd=.3, draw_quantiles=c(.5),
              position=position_dodge(width=.7)) +
  theme_bw() + facet_grid(false_pos_neg~., scales="free") +
  scale_fill_manual(values=pal)

ggsave(paste0(df_path, "draft_fig2.pdf"))


mse_counts_cmp = stats_compare %>%
  dplyr::select(-mse_expos, -mse_expos_missing, -dplyr::contains("cosine")) %>% 
  plot_performance(fill="penalty", facet="~variable", pal=pal)

cosine_expos_sigs_cmp = stats_compare %>%
  dplyr::select(-dplyr::contains("mse"), -cosine_expos) %>% 
  plot_performance(fill="penalty", facet="~variable", pal=pal)




# Panels #####

panelA = patchwork::wrap_plots(K_ratio + labs(tag="A"),
                               patchwork::wrap_plots(mse_counts, cosine_sigs_expos,
                                                     guides="collect"),
                               design="ABBB") &
  theme(legend.position="bottom")

panelB = patchwork::wrap_plots(K_ratio_cmp + labs(tag="B"), 
                               patchwork::wrap_plots(mse_counts_cmp, cosine_expos_sigs_cmp, 
                                                     guides="collect"),
                               design="ABBB") &
  theme(legend.position="bottom")


panelC = clustering + labs(tag="C")


patchwork::wrap_plots(panelA, panelB, panelC, ggplot() + labs(tag="D"), design="AAAC\nBBBD")
ggsave(paste0(df_path, "draft_fig2.pdf"), height=8, width=12)
ggsave(paste0(df_path, "draft_fig2.png"), height=8, width=12)




