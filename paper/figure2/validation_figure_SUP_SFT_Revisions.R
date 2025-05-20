library(magrittr)
library(ggplot2)
source("~/GitHub/bascule_validation/synthetic_data/aux_fns/eval_aux_fns.R")
source("~/GitHub/bascule_validation/synthetic_data/aux_fns/plots_aux_fns.R")

df_path = "~/Dropbox/dropbox_shared/2022. Basilica/simulations/stats_dataframes/revisions/"
stats_bascule = readRDS(paste0(df_path, "stats_matched.2011.compare_LAST.SigFitTest.nofits.Rds")) %>% 
  compute_quantiles(colname="K_true") %>% dplyr::filter(penalty=="BASCULE") %>% 
  dplyr::mutate(nMuts_label="Number of mutations")
stats_compare = readRDS(paste0(df_path, "stats_matched.2011.compare_LAST.SigFitTest.nofits.Rds")) %>% 
  compute_quantiles(colname="K_true") %>% 
  dplyr::filter(type=="SBS") %>% 
  dplyr::mutate(nMuts_label="Number of mutations")

theme_legend = theme(legend.text=element_text(size=7.4),
                     legend.title=element_text(size=10.6),
                     legend.key.size=unit(0.1,"cm"),
                     legend.key.height=unit(0.5,"cm"),
                     legend.key.width=unit(0.5,"cm"))

theme_text = theme(axis.title=element_text(size=7.4),
                   axis.text=element_text(size=5.35), 
                   plot.title=element_text(size=10.6),
                   plot.subtitle=element_text(size=7.4),
                   strip.text=element_text(size=5.35), 
                   plot.caption=element_text(size=5.35))

theme_nofacet = theme(strip.background=element_blank(), 
                      strip.text=element_blank())

# SF S1 - BASCULE ####

pal_k_true = pal=c("tan2","#8FBC8B","thistle2")

# panelA -> bascule validation: K_ratio by complexity, MSE counts, cosine sigs, NMI
K_ratio = stats_bascule %>% 
  dplyr::rowwise() %>% 
  dplyr::mutate(FN=length(assigned_missing$missing_fn),
                FP=length(assigned_missing$added_fp),
                TP=length(assigned_missing$assigned_tp)) %>%
  dplyr::mutate(Recall=TP / (TP + FN)) %>% 
  
  dplyr::select(-K_input_ratio, -K_dn_ratio) %>% 
  plot_K(fill="K_true_cat", pal=pal_k_true, pattern="Recall$",
         facet="source~nMuts_label+nMuts")
K_ratio

precision = stats_bascule %>% 
  dplyr::rowwise() %>% 
  dplyr::mutate(FN=length(assigned_missing$missing_fn),
                FP=length(assigned_missing$added_fp),
                TP=length(assigned_missing$assigned_tp)) %>%
  dplyr::mutate(Precision=TP / (TP + FP)) %>% 
  
  compute_quantiles(colname="K_true") %>% 
  dplyr::select(-K_input_ratio, -K_dn_ratio) %>% 
  plot_K(fill="K_true_cat", pal=pal_k_true, pattern="Precision$",
         facet="source~nMuts_label+nMuts")
precision

## accuracy : (TP + TN) / (TP + TN + FP + FN)
## precision : TP / (TP + FP)
## recall : TP / (TP + FN)

mse_counts = stats_bascule %>% 
  dplyr::select(-mse_expos, -mse_expos_missing, -dplyr::contains("cosine")) %>% 
  plot_performance(fill="K_true_cat", pal=pal_k_true,
                   facet="source~nMuts_label+nMuts")
mse_counts

cosine_sigs = stats_bascule %>% 
  dplyr::select(-dplyr::contains("mse"), -cosine_expos_missing,
                -cosine_expos) %>% 
  plot_performance(fill="K_true_cat", pal=pal_k_true,
                   facet="source~nMuts_label+nMuts")
cosine_sigs

cosine_expos = stats_bascule %>% 
  dplyr::select(-dplyr::contains("mse"), -cosine_expos_missing,
                -cosine_sigs) %>% 
  plot_performance(fill="K_true_cat", pal=pal_k_true,
                   facet="source~nMuts_label+nMuts")
cosine_expos

clustering = stats_bascule %>% 
  dplyr::filter(type=="SBS") %>% 
  dplyr::rename(ARI=ari, NMI=nmi) %>% 
  dplyr::select(-ARI) %>% 
  plot_performance_clustering(fill="K_true_cat", 
                              facet="source~nMuts_label+nMuts", pal=pal_k_true)
clustering


# SF S2 - comparison #####

# pal_methods = RColorBrewer::brewer.pal(3, name="Dark2")
pal_methods = c("#7fb3d5", "#FF8C00", "#8FBC8B", "#DB7093", RColorBrewer::brewer.pal(4, name="Dark2")) %>% 
  setNames(c("BASCULE", "SigProfiler", "SparseSignatures","SignatureToolsLib_E","KMeans","KL-KMeans","JS-Spectral","SignatureToolsLib"))


stats_compare = stats_compare %>% 
  dplyr::mutate(penalty=replace(penalty,penalty=="Basilica","BASCULE"))


recall_cmp = stats_compare %>% dplyr::rowwise() %>%
  dplyr::mutate(FN=length(assigned_missing$missing_fn),
                FP=length(assigned_missing$added_fp),
                TP=length(assigned_missing$assigned_tp)) %>%
  dplyr::mutate(Recall=TP / (TP + FN),
                label="# signatures") %>% 
  tidyr::pivot_longer(cols=c("Recall"), 
                      names_to="false_pos_neg") %>% 
  ggplot() + 
  geom_violin(aes(x=factor(N), y=value, color=factor(penalty), fill=factor(penalty)), 
              lwd=.5, alpha=0.7, draw_quantiles=c(.5),
              position=position_dodge(width=.7)) +
  ggh4x::facet_nested(label + K_true_cat ~ source + nMuts_label + nMuts, scales="free") +
  scale_fill_manual(values=pal_methods) +
  scale_color_manual(values=pal_methods) +
  theme_bw()
recall_cmp


precision_cmp = stats_compare %>% dplyr::rowwise() %>%
  dplyr::mutate(FN=length(assigned_missing$missing_fn),
                FP=length(assigned_missing$added_fp),
                TP=length(assigned_missing$assigned_tp)) %>%
  dplyr::mutate(Precision=TP / (TP + FP),
                label="# signatures") %>% 
  tidyr::pivot_longer(cols=c("Precision"), 
                      names_to="false_pos_neg") %>% 
  ggplot() + 
  
  geom_violin(aes(x=factor(N), y=value, color=factor(penalty), fill=factor(penalty)),
              lwd=.5, alpha=0.7, draw_quantiles=c(.5),
              position=position_dodge(width=.7)) +
  ggh4x::facet_nested(label + K_true_cat ~ source + nMuts_label + nMuts, scales="free") +
  scale_fill_manual(values=pal_methods) +
  scale_color_manual(values=pal_methods) +
  theme_bw()
precision_cmp


mse_counts_cmp = stats_compare %>%
  dplyr::mutate(label="# signatures") %>% 
  dplyr::select(-mse_expos, -mse_expos_missing, -dplyr::contains("cosine")) %>% 
  plot_performance(fill="penalty", facet="label+K_true_cat ~ source + nMuts_label + nMuts", pal=pal_methods)
mse_counts_cmp

cosine_expos_cmp = stats_compare %>%
  dplyr::mutate(label="# signatures") %>% 
  dplyr::select(-dplyr::contains("mse"), -cosine_expos_missing, -cosine_sigs) %>% 
  plot_performance(fill="penalty", facet="label+K_true_cat ~ source + nMuts_label + nMuts", pal=pal_methods)
cosine_expos_cmp

cosine_expos_missing_cmp = stats_compare %>%
  dplyr::mutate(label="# signatures") %>% 
  dplyr::select(-dplyr::contains("mse"), -cosine_expos, -cosine_sigs) %>% 
  plot_performance(fill="penalty", facet="label+K_true_cat ~ source + nMuts_label + nMuts", pal=pal_methods)
cosine_expos_missing_cmp

cosine_sigs_cmp = stats_compare %>%
  dplyr::mutate(label="# signatures") %>% 
  dplyr::select(-dplyr::contains("mse"), -cosine_expos, -cosine_expos_missing) %>% 
  plot_performance(fill="penalty", facet="label+K_true_cat ~ source + nMuts_label + nMuts", pal=pal_methods)
cosine_sigs_cmp


# SF S3 - runtimes ####
runtime_path = "~/Dropbox/dropbox_shared/2022. Basilica/simulations/runtimes/SigFitTest/"
times_sigpr = read.csv(file.path(runtime_path, "sigprofiler_exectimes.csv")) %>%
  dplyr::mutate(tool="SigProfiler") %>% tibble::as_tibble() %>%
  dplyr::mutate(execution_time=stringr::str_replace_all(execution_time, " 0:", "00:")) %>%
  dplyr::mutate(execution_time=stringr::str_remove_all(execution_time, " ")) %>%
  dplyr::mutate(execution_time=lubridate::period_to_seconds(lubridate::hms(execution_time))) %>%
  dplyr::rename(simulation_name=simulation) %>%
  dplyr::select(simulation_name, execution_time, tool)
times_sparsesig = read.csv(file.path(runtime_path, "sparsesignatures_exectimes.csv")) %>%
  dplyr::rename(execution_time=total_mins, simulation_name=name) %>%
  dplyr::mutate(tool="SparseSignatures",
                simulation_name=stringr::str_remove_all(simulation_name, ".Rds")) %>%
  tibble::as_tibble() %>%
  dplyr::select(simulation_name, execution_time, tool)
times_sigtoolslib = read.csv(file.path(runtime_path, "signaturetoolslib_exectimes.csv")) %>% 
  dplyr::mutate(tool="SignatureToolsLib",
                simulation_name=stringr::str_remove_all(simulation_name, ".Rds")) %>% 
  tibble::as_tibble() %>% 
  dplyr::select(simulation_name, execution_time, tool)
times_sigtoolslib_E = read.csv(file.path(runtime_path, "signaturetoolslib_E_exectimes.csv")) %>% 
  dplyr::mutate(tool="SignatureToolsLib_E",
                simulation_name=stringr::str_remove_all(simulation_name, ".Rds")) %>% 
  tibble::as_tibble() %>% 
  dplyr::select(simulation_name, execution_time, tool)
times_bascule = read.csv(file.path(runtime_path, "bascule_exectimes.csv")) %>% 
  tidyr::separate(simulation_name, into=c("remove","simulation_name"), sep="/") %>% 
  dplyr::select(-remove) %>% 
  dplyr::rename(execution_time_bascule=execution_time_sbs) %>% tibble::as_tibble() %>% 
  dplyr::mutate(tool="BASCULE", execution_time=execution_time_bascule) %>%
  dplyr::select(simulation_name, execution_time, execution_time_bascule, tool)


runtimes_cmp = dplyr::bind_rows(
  # times_sigpr, 
  # times_sparsesig, 
  times_sigtoolslib,
  times_sigtoolslib_E,
  times_bascule
  )  %>% 
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
  geom_boxplot(aes(y=execution_time, x=factor(N), color=tool, fill=tool),
               alpha=0.7, lwd=0.5) +
  # facet_wrap(~processor, scales="free") +
  scale_fill_manual(values=pal_methods) +
  scale_color_manual(values=pal_methods) +
  theme_bw()
runtimes_cmp


# Figure BASCULE ####

panelAa = K_ratio + xlab("# samples") + ylab("Recall") + 
  labs(tag="A", title="Signatures detection accuracy (recall)",
       subtitle="Recall of identified signatures") +
  guides(fill=guide_legend(title="# signatures"),
         color=guide_legend(title="# signatures")) + 
  theme(legend.position="bottom")

panelAb = precision + xlab("# samples") + ylab("Precision") + 
  labs(tag="B", title="Signatures detection accuracy (precision)",
       subtitle="Precision of identified signatures") +
  guides(fill=guide_legend(title="# signatures"),
         color=guide_legend(title="# signatures")) + 
  theme(legend.position="bottom")

panelB = mse_counts + 
  labs(tag="C", title="Reconstruction error",
       subtitle="MSE between inferred and true mutation counts") +
  xlab("# samples") + theme(legend.position="bottom") +
  guides(fill=guide_legend(title="# signatures"),
         color=guide_legend(title="# signatures"))

panelC = cosine_sigs + 
  labs(tag="D", title="Signatures quality",
       subtitle="CS between inferred and true signature profiles") +
  xlab("# samples") + theme(legend.position="bottom") +
  guides(fill=guide_legend(title="# signatures"),
         color=guide_legend(title="# signatures"))

panelD = cosine_expos + 
  labs(tag="E", title="Exposures quality",
       subtitle="CS between inferred and true exposures") +
  xlab("# samples") + theme(legend.position="bottom") +
  guides(fill=guide_legend(title="# signatures"),
         color=guide_legend(title="# signatures"))

panelE = clustering + xlab("# samples") + ylab("NMI") +
  labs(tag="F", title="Clustering accuracy",
       subtitle="NMI between inferred and true assignments") +
  guides(fill=guide_legend(title="# signatures"),
         color=guide_legend(title="# signatures"))

figure1 = patchwork::wrap_plots(panelAa, panelAb, 
                                panelB, panelC, 
                                panelD, panelE, 
                                # guides="collect",
                                ncol=2) &
  theme_text & theme_legend & theme(legend.position="right") &
  patchwork::plot_annotation(tag_levels="A")
figure1

# ggsave(filename="paper/figure2/revisions/figure2_SUP1.pdf", plot=figure1,
#        width=210, height=210, units="mm")
ggsave(filename="paper/figure2/revisions/SigFitTest/figure2_SUP1.png", plot=figure1,
       width=210, height=210, units="mm")

# Figure comparison SBS #####

panelAa = recall_cmp + xlab("# samples") + ylab("Recall") + 
  labs(tag="A", title="Signatures detection accuracy (recall)",
       subtitle="Recall of identified signatures") +
  guides(fill=guide_legend(title="Method"),
         color=guide_legend(title="Method")) + 
  theme(legend.position="bottom")

panelAb = precision_cmp + xlab("# samples") + ylab("Precision") + 
  labs(tag="B", title="Signatures detection accuracy (precision)",
       subtitle="Precision of identified signatures") +
  guides(fill=guide_legend(title="Method"),
         color=guide_legend(title="Method")) + 
  theme(legend.position="bottom")

panelB = mse_counts_cmp + 
  labs(tag="C", title="Reconstruction error",
       subtitle="MSE between inferred and true mutation counts") +
  xlab("# samples") + theme(legend.position="bottom") +
  guides(fill=guide_legend(title="Method"),
         color=guide_legend(title="Method"))

panelC = cosine_sigs_cmp + 
  labs(tag="D", title="Signatures quality",
       subtitle="CS between inferred and true signature profiles") +
  xlab("# samples") + theme(legend.position="bottom") +
  guides(fill=guide_legend(title="Method"),
         color=guide_legend(title="Method"))

panelD = cosine_expos_cmp + 
  labs(tag="E", title="Exposures quality",
       subtitle="CS between inferred and true exposures of matched signatures") +
  xlab("# samples") + theme(legend.position="bottom") +
  guides(fill=guide_legend(title="Method"),
         color=guide_legend(title="Method"))

panelE = cosine_expos_missing_cmp + xlab("# samples") + 
  ylab("Cosine similarity") +
  labs(tag="F", title="Exposures quality (all signatures)",
       subtitle="CS between inferred and true exposures of all signatures") +
  guides(fill=guide_legend(title="Method"),
         color=guide_legend(title="Method"))

figure2 = patchwork::wrap_plots(panelAa, panelAb, 
                                panelB, panelC, 
                                panelD, panelE, 
                                # guides="collect",
                                ncol=2) &
  theme_text & theme_legend & theme(legend.position="bottom") &
  patchwork::plot_annotation(tag_levels="A")
figure2

# ggsave(filename="paper/figure2/revisions/SigFitTest/figure2_SUP2.pdf", plot=figure2,
#        width=210, height=250, units="mm")
# ggsave(filename="paper/figure2/revisions/SigFitTest/figure2_SUP2.png", plot=figure2,
#        width=210, height=250, units="mm")
ggsave(filename="paper/figure2/revisions/SigFitTest/figure2_SUP2.png", plot=figure2,
       width=210*1.5, height=250*1.5, units="mm")


# Figure runtime ####

runtimes = runtimes_cmp + xlab("# samples") + ylab("Time (minutes)") +
  labs(title="Fit runtime",
       subtitle="Execution time (minutes) required to fit the model") +
  guides(fill=guide_legend(title="Method"),
         color=guide_legend(title="Method")) +
  theme_text + theme_legend + theme(legend.position="bottom")

runtimes

# ggsave(filename="paper/figure2/revisions/SigFitTest/figure2_SUP3.pdf", plot=runtimes,
#        width=120, height=100, units="mm")
ggsave(filename="paper/figure2/revisions/SigFitTest/figure2_SUP3.png", plot=runtimes,
       width=120, height=100, units="mm")
