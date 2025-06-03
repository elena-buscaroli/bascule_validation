library(magrittr)
library(dplyr)
library(ggplot2)
source("~/GitHub/bascule_validation/synthetic_data/aux_fns/eval_aux_fns.R")
source("~/GitHub/bascule_validation/synthetic_data/aux_fns/plots_aux_fns.R")

df_path = "~/Dropbox/dropbox_shared/2022. Basilica/simulations/stats_dataframes/revisions/"
# stats_bascule = readRDS(paste0(df_path, "stats_matched.2011.compare_LAST.SigFitTest.nofits.Rds")) %>% 
#   compute_quantiles(colname="K_true") %>% dplyr::filter(penalty=="BASCULE") %>% 
#   dplyr::mutate(nMuts_label="Number of mutations")
stats_compare = readRDS(paste0(df_path, "stats_matched.2011.compare_LAST.SigFitTest.nofits.Rds")) %>% 
  filter(penalty!="SignatureToolsLib") %>%
  mutate(penalty=stringr::str_replace_all(penalty, "SignatureToolsLib", "FitMS")) %>% 
  compute_quantiles(colname="K_true") %>% 
  dplyr::filter(type=="SBS", source=="WGS") %>% 
  dplyr::mutate(nMuts_label="Number of mutations") %>% 
  mutate(penalty=factor(penalty, levels=c("BASCULE", "SigProfiler", "SparseSignatures", "FitMS_E")))

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

# # SF S1 - BASCULE ####
# 
# pal_k_true = pal =c("tan2","#8FBC8B","thistle2")
# 
# # panelA -> bascule validation: K_ratio by complexity, MSE counts, cosine sigs, NMI
# K_ratio = stats_bascule %>% 
#   dplyr::rowwise() %>% 
#   dplyr::mutate(FN=length(assigned_missing$missing_fn),
#                 FP=length(assigned_missing$added_fp),
#                 TP=length(assigned_missing$assigned_tp)) %>%
#   dplyr::mutate(Recall=TP / (TP + FN)) %>% 
#   
#   dplyr::select(-K_input_ratio, -K_dn_ratio) %>% 
#   # plot_K(fill="K_true_cat", pal=pal_k_true, pattern="Recall$", facet="source~nMuts_label+nMuts")
#   ggplot(aes(x=factor(N), y=Recall, fill=K_true_cat, color=K_true_cat)) +
#   geom_violin(draw_quantiles=c(.5), lwd=.5, alpha=0.7,
#               position=position_dodge(width=.7), show.legend=T) +
#   ggh4x::facet_nested(source~nMuts_label+nMuts) +
#   scale_color_manual(values=pal_k_true) +
#   scale_fill_manual(values=pal_k_true) +
#   theme_bw()
# K_ratio
# 
# precision = stats_bascule %>% 
#   dplyr::rowwise() %>% 
#   dplyr::mutate(FN=length(assigned_missing$missing_fn),
#                 FP=length(assigned_missing$added_fp),
#                 TP=length(assigned_missing$assigned_tp)) %>%
#   dplyr::mutate(Precision=TP / (TP + FP)) %>% 
#   
#   compute_quantiles(colname="K_true") %>% 
#   dplyr::select(-K_input_ratio, -K_dn_ratio) %>% 
#   # plot_K(fill="K_true_cat", pal=pal_k_true, pattern="Precision$",
#   #        facet="source~nMuts_label+nMuts")
#   ggplot(aes(x=factor(N), y=Precision, fill=K_true_cat, color=K_true_cat)) +
#   geom_violin(draw_quantiles=c(.5), lwd=.5, alpha=0.7,
#               position=position_dodge(width=.7), show.legend=T) +
#   ggh4x::facet_nested(source~nMuts_label+nMuts) +
#   scale_color_manual(values=pal_k_true) +
#   scale_fill_manual(values=pal_k_true) +
#   theme_bw()
# precision
# 
# ## accuracy : (TP + TN) / (TP + TN + FP + FN)
# ## precision : TP / (TP + FP)
# ## recall : TP / (TP + FN)
# 
# mse_counts = stats_bascule %>% 
#   dplyr::select(-mse_expos, -mse_expos_missing, -dplyr::contains("cosine")) %>% 
#   # plot_performance(fill="K_true_cat", pal=pal_k_true, facet="source~nMuts_label+nMuts")
#   ggplot(aes(x=factor(N), y=mse_counts, fill=K_true_cat, color=K_true_cat)) +
#   geom_boxplot(outlier.shape=NA, lwd=0.5, alpha=0.7, width=0.5) +
#   ggh4x::facet_nested(source~nMuts_label+nMuts) +
#   scale_color_manual(values=pal_k_true) +
#   scale_fill_manual(values=pal_k_true) +
#   theme_bw()
# mse_counts
# 
# cosine_sigs = stats_bascule %>% 
#   dplyr::select(-dplyr::contains("mse"), -cosine_expos_missing,
#                 -cosine_expos) %>% 
#   # plot_performance(fill="K_true_cat", pal=pal_k_true, facet="source~nMuts_label+nMuts")
#   ggplot(aes(x=factor(N), y=cosine_sigs, fill=K_true_cat, color=K_true_cat)) +
#   geom_boxplot(outlier.shape=NA, lwd=0.5, alpha=0.7, width=0.5) +
#   ggh4x::facet_nested(source~nMuts_label+nMuts) +
#   scale_color_manual(values=pal_k_true) +
#   scale_fill_manual(values=pal_k_true) +
#   theme_bw()
# cosine_sigs
# 
# cosine_expos = stats_bascule %>% 
#   dplyr::select(-dplyr::contains("mse"), -cosine_expos_missing,
#                 -cosine_sigs) %>% 
#   # plot_performance(fill="K_true_cat", pal=pal_k_true, facet="source~nMuts_label+nMuts")
#   ggplot(aes(x=factor(N), y=cosine_expos, fill=K_true_cat, color=K_true_cat)) +
#   geom_boxplot(outlier.shape=NA, lwd=0.5, alpha=0.7, width=0.5) +
#   ggh4x::facet_nested(source~nMuts_label+nMuts) +
#   scale_color_manual(values=pal_k_true) +
#   scale_fill_manual(values=pal_k_true) +
#   theme_bw()
# cosine_expos
# 
# # clustering = stats_bascule %>% 
# #   dplyr::filter(type=="SBS") %>% 
# #   dplyr::rename(NMI=nmi) %>% 
# #   dplyr::select(-ARI) %>% 
# #   plot_performance_clustering(fill="K_true_cat", 
# #                               facet="source~nMuts_label+nMuts", pal=pal_k_true)
# 
# clustering = stats_bascule %>% 
#   dplyr::filter(type=="SBS") %>% 
#   compute_quantiles(colname="K_true") %>% 
#   dplyr::select(idd, N, G, source, nMuts_label, nMuts, K_true_cat, starts_with("nmi")) %>% 
#   
#   tidyr::pivot_longer(cols=starts_with("nmi"), values_to="nmi", names_to="Method") %>% 
#   dplyr::mutate(Method=dplyr::case_when(Method=="nmi" ~ "BASCULE",
#                                         Method=="nmi_KM" ~ "KMeans",
#                                         Method=="nmi_KL" ~ "KL-KMeans",
#                                         Method=="nmi_JS" ~ "JS-Spectral")) %>% 
#   
#   dplyr::mutate(Method=reorder(Method, nmi, mean, decreasing=T)) %>% 
#   
#   ggplot(aes(x=factor(N), y=nmi, fill=Method, color=Method)) +
#   # stat_summary(aes(group=Method), position=position_dodge(width=0.2), 
#   #              fun.data="mean_cl_boot", show.legend=T, size=.2) +
#   # stat_summary(aes(group=Method), position=position_dodge(width=0.2), 
#   #              fun.data="mean_cl_boot", show.legend=T,
#   #              geom="line", linewidth=1) +
#   geom_boxplot(outlier.shape=NA, lwd=.5, alpha=.7, width=0.5) +
#   ggh4x::facet_nested(source ~ nMuts_label+nMuts) +
#   scale_fill_manual(values=pal_methods, 
#                     breaks=names(pal_methods),
#                     name="Method") +
#   scale_color_manual(values=pal_methods, 
#                      breaks=names(pal_methods),
#                      name="Method") +
#   theme_bw() + ylim(0, 1)
# 
# clustering
# 

# SF S2 - comparison #####

# pal_methods = RColorBrewer::brewer.pal(3, name="Dark2")
# pal_methods = c("#7fb3d5", "#FF8C00", "#8FBC8B", "#DB7093", RColorBrewer::brewer.pal(4, name="Dark2")) %>% 
#   setNames(c("BASCULE", "SigProfiler", "SparseSignatures","SignatureToolsLib_E","KMeans","KL-KMeans","JS-Spectral","SignatureToolsLib"))
pal_methods = c("#7fb3d5", "#FF8C00", "#8FBC8B", "#DB7093", RColorBrewer::brewer.pal(4, name="Dark2")) %>%
  setNames(c("BASCULE", "SigProfiler", "SparseSignatures","FitMS_E","KMeans","KL-KMeans","JS-Spectral","FitMS"))


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
  ggplot(aes(x=factor(N), y=value, color=factor(penalty), fill=factor(penalty))) + 
  # geom_violin(aes(x=factor(N), y=value, color=factor(penalty), fill=factor(penalty)), 
  #             lwd=.5, alpha=0.7, draw_quantiles=c(.5),
  #             position=position_dodge(width=.7)) +
  stat_summary(aes(group=penalty), fun.data="mean_cl_boot", position=position_dodge(width=.15), 
               show.legend=TRUE, geom="line", linewidth=1) +
  stat_summary(aes(group=penalty), fun.data="mean_cl_boot", position=position_dodge(width=.15), 
               show.legend=TRUE, size=0.2) +
  ggh4x::facet_nested(label + K_true_cat ~ nMuts_label + nMuts, scales="fixed") +
  scale_fill_manual(values=pal_methods) +
  scale_color_manual(values=pal_methods) + ylim(0,1) +
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
  ggplot(aes(x=factor(N), y=value, color=factor(penalty), fill=factor(penalty))) + 
  
  # geom_violin(aes(x=factor(N), y=value, color=factor(penalty), fill=factor(penalty)),
  #             lwd=.5, alpha=0.7, draw_quantiles=c(.5),
  #             position=position_dodge(width=.7)) +
  stat_summary(aes(group=penalty), fun.data="mean_cl_boot", position=position_dodge(width=.15), 
               show.legend=TRUE, geom="line", linewidth=1) +
  stat_summary(aes(group=penalty), fun.data="mean_cl_boot", position=position_dodge(width=.15), 
               show.legend=TRUE, size=0.2) +
  ggh4x::facet_nested(label + K_true_cat ~ nMuts_label + nMuts, scales="free") +
  scale_fill_manual(values=pal_methods) +
  scale_color_manual(values=pal_methods) + ylim(0,1) +
  theme_bw()
precision_cmp


mse_counts_cmp = stats_compare %>%
  dplyr::mutate(label="# signatures") %>% 
  dplyr::select(-mse_expos, -mse_expos_missing, -dplyr::contains("cosine")) %>% 
  # plot_performance(fill="penalty", facet="label+K_true_cat ~ nMuts_label + nMuts", pal=pal_methods)

  ggplot(aes(x=factor(N), y=mse_counts, fill=penalty, color=penalty)) +
  stat_summary(aes(group=penalty), fun.data="mean_cl_boot", position=position_dodge(width=.15), 
               show.legend=TRUE, geom="line", linewidth=1) +
  stat_summary(aes(group=penalty), fun.data="mean_cl_boot", position=position_dodge(width=.15), 
               show.legend=TRUE, size=0.2) +
  ggh4x::facet_nested(label+K_true_cat ~ nMuts_label + nMuts) +
  scale_fill_manual(values=pal_methods) +
  scale_color_manual(values=pal_methods) +
  theme_bw()
mse_counts_cmp

cosine_expos_cmp = stats_compare %>%
  dplyr::mutate(label="# signatures") %>% 
  dplyr::select(-dplyr::contains("mse"), -cosine_expos_missing, -cosine_sigs) %>% 
  # plot_performance(fill="penalty", facet="label+K_true_cat ~ source + nMuts_label + nMuts", pal=pal_methods)
  
  ggplot(aes(x=factor(N), y=cosine_expos, fill=penalty, color=penalty)) +
  stat_summary(aes(group=penalty), fun.data="mean_cl_boot", position=position_dodge(width=.15), 
               show.legend=TRUE, geom="line", linewidth=1) +
  stat_summary(aes(group=penalty), fun.data="mean_cl_boot", position=position_dodge(width=.15), 
               show.legend=TRUE, size=0.2) +
  ggh4x::facet_nested(label+K_true_cat ~ nMuts_label + nMuts) +
  scale_fill_manual(values=pal_methods) +
  scale_color_manual(values=pal_methods) +
  theme_bw()
cosine_expos_cmp

cosine_expos_missing_cmp = stats_compare %>%
  dplyr::mutate(label="# signatures") %>% 
  dplyr::select(-dplyr::contains("mse"), -cosine_expos, -cosine_sigs) %>% 
  # plot_performance(fill="penalty", facet="label+K_true_cat ~ source + nMuts_label + nMuts", pal=pal_methods)
  
  ggplot(aes(x=factor(N), y=cosine_expos_missing, fill=penalty, color=penalty)) +
  stat_summary(aes(group=penalty), fun.data="mean_cl_boot", position=position_dodge(width=.15), 
               show.legend=TRUE, geom="line", linewidth=1) +
  stat_summary(aes(group=penalty), fun.data="mean_cl_boot", position=position_dodge(width=.15), 
               show.legend=TRUE, size=0.2) +
  ggh4x::facet_nested(label+K_true_cat ~ nMuts_label + nMuts) +
  scale_fill_manual(values=pal_methods) +
  scale_color_manual(values=pal_methods) +
  theme_bw()
cosine_expos_missing_cmp

cosine_sigs_cmp = stats_compare %>%
  dplyr::mutate(label="# signatures") %>% 
  dplyr::select(-dplyr::contains("mse"), -cosine_expos, -cosine_expos_missing) %>% 
  # plot_performance(fill="penalty", facet="label+K_true_cat ~ source + nMuts_label + nMuts", pal=pal_methods)
  
  ggplot(aes(x=factor(N), y=cosine_sigs, fill=penalty, color=penalty)) +
  stat_summary(aes(group=penalty), fun.data="mean_cl_boot", position=position_dodge(width=.15), 
               show.legend=TRUE, geom="line", linewidth=1) +
  stat_summary(aes(group=penalty), fun.data="mean_cl_boot", position=position_dodge(width=.15), 
               show.legend=TRUE, size=0.2) +
  ggh4x::facet_nested(label+K_true_cat ~ nMuts_label + nMuts) +
  scale_fill_manual(values=pal_methods) +
  scale_color_manual(values=pal_methods) +
  theme_bw()
cosine_sigs_cmp


# # SF S3 - runtimes ####
# runtime_df = readRDS(file.path(df_path, "runtime_SigFitTest.Rds")) %>% 
#   mutate(tool=stringr::str_replace_all(tool, "SignatureToolsLib", "FitMS"))
# 
# runtimes_cmp = runtime_df %>% 
#   dplyr::filter(tool != "FitMS") %>%
#   
#   rowwise() %>% 
#   dplyr::mutate(N=strsplit(simulation_name, "[.]")[[1]][2] %>% 
#                   stringr::str_remove_all("N") %>% as.numeric()) %>% 
#   dplyr::select(simulation_name, execution_time, tool, N) %>% 
#   
#   # dplyr::group_by(tool, N) %>% 
#   # dplyr::filter(execution_time < boxplot.stats(execution_time)$stats[5]) %>% 
#   
#   ggplot() +
#   geom_boxplot(aes(y=execution_time, x=factor(N), color=tool, fill=tool),
#                alpha=0.7, lwd=0.5) +
#   # facet_wrap(~processor, scales="free") +
#   scale_fill_manual(values=pal_methods, breaks=names(pal_methods)) +
#   scale_color_manual(values=pal_methods, breaks=names(pal_methods)) +
#   theme_bw()
# runtimes_cmp
# 

# # Figure BASCULE ####
# 
# panelAa = K_ratio + xlab("# samples") + ylab("Recall") + 
#   labs(tag="A", title="Signatures detection accuracy (recall)",
#        subtitle="Recall of identified signatures") +
#   guides(fill=guide_legend(title="# signatures"),
#          color=guide_legend(title="# signatures")) + 
#   theme(legend.position="bottom")
# 
# panelAb = precision + xlab("# samples") + ylab("Precision") + 
#   labs(tag="B", title="Signatures detection accuracy (precision)",
#        subtitle="Precision of identified signatures") +
#   guides(fill=guide_legend(title="# signatures"),
#          color=guide_legend(title="# signatures")) + 
#   theme(legend.position="bottom")
# 
# panelB = mse_counts + 
#   labs(tag="C", title="Reconstruction error",
#        subtitle="MSE between inferred and true mutation counts") +
#   xlab("# samples") + theme(legend.position="bottom") +
#   guides(fill=guide_legend(title="# signatures"),
#          color=guide_legend(title="# signatures"))
# 
# panelC = cosine_sigs + 
#   labs(tag="D", title="Signatures quality",
#        subtitle="CS between inferred and true signature profiles") +
#   xlab("# samples") + theme(legend.position="bottom") +
#   guides(fill=guide_legend(title="# signatures"),
#          color=guide_legend(title="# signatures"))
# 
# panelD = cosine_expos + 
#   labs(tag="E", title="Exposures quality",
#        subtitle="CS between inferred and true matched exposures") +
#   xlab("# samples") + theme(legend.position="bottom") +
#   guides(fill=guide_legend(title="# signatures"),
#          color=guide_legend(title="# signatures"))
# 
# panelE = clustering + xlab("# samples") + ylab("NMI") +
#   labs(tag="F", title="Clustering accuracy",
#        subtitle="NMI between inferred and true assignments") +
#   guides(fill=guide_legend(title="# signatures"),
#          color=guide_legend(title="# signatures"))
# 
# figure1 = patchwork::wrap_plots(panelAa, panelAb, 
#                                 panelB, panelC, 
#                                 panelD, panelE, 
#                                 # guides="collect",
#                                 ncol=2) &
#   theme_text & theme_legend & theme(legend.position="right") &
#   patchwork::plot_annotation(tag_levels="A")
# figure1
# 
# ggsave(filename="paper/figure2/revisions/SigFitTest/figure2_SUPX.png", plot=figure1,
#        width=210, height=210, units="mm", device=png, family="Helvetica")

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
  xlab("# samples") + ylab("Mean squared error") + theme(legend.position="bottom") +
  guides(fill=guide_legend(title="Method"),
         color=guide_legend(title="Method"))

panelC = cosine_sigs_cmp + 
  labs(tag="D", title="Signatures quality",
       subtitle="CS between inferred and true signature profiles") +
  xlab("# samples") + ylab("Cosine similarity") + theme(legend.position="bottom") +
  guides(fill=guide_legend(title="Method"),
         color=guide_legend(title="Method"))

panelD = cosine_expos_cmp + 
  labs(tag="E", title="Exposures quality",
       subtitle="CS between inferred and true exposures of matched signatures") +
  xlab("# samples") + ylab("Cosine similarity") + theme(legend.position="bottom") +
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


ggsave(filename="paper/figure2/revisions/SigFitTest/figure2_SUP4.png", plot=figure2,
       width=210*1.2, height=260*1.2, units="mm", device=png, family="Helvetica")


# Figure runtime ####

runtimes = runtimes_cmp + xlab("# samples") + ylab("Time (minutes)") +
  labs(title="Fit runtime",
       subtitle="Execution time (minutes) required to fit the model") +
  guides(fill=guide_legend(title="Method"),
         color=guide_legend(title="Method")) +
  theme_text + theme_legend + theme(legend.position="right")

runtimes

ggsave(filename="paper/figure2/revisions/SigFitTest/figure2_SUPXX.png", plot=runtimes,
       width=150, height=100, units="mm", device=png, family="Helvetica")
