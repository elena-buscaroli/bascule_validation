library(magrittr)
library(ggplot2)
devtools::load_all("~/GitHub/bascule/")
source("~/GitHub/bascule_validation/synthetic_data/aux_fns/eval_aux_fns.R")
source("~/GitHub/bascule_validation/synthetic_data/aux_fns/plots_aux_fns.R")

df_path = "~/Dropbox/dropbox_shared/2022. Basilica/simulations/stats_dataframes/"
stats_bascule = readRDS(paste0(df_path, "stats_matched.2011_KM.Rds")) %>% 
  compute_quantiles(colname="K_true")
stats_compare = readRDS(paste0(df_path, "stats_matched.2011.compare.Rds")) %>% 
  compute_quantiles(colname="K_true") %>% 
  dplyr::filter(type=="SBS") %>% 
  dplyr::mutate(penalty=replace(penalty, penalty=="Basilica", "BASCULE"))

plots = list()

pal_k_true = wesanderson::wes_palette("Cavalcanti1", 6, type="continuous")[c(1,4,3)]
pal_methods = c("#7fb3d5", "#FF8C00", "#8FBC8B", "#DB7093") %>% setNames(c("BASCULE", "SigProfiler", "SparseSignatures","KMeans"))


# Basilica ####

plots[["nmi"]] = stats_bascule %>% 
  dplyr::filter(type=="SBS") %>% 
  compute_quantiles(colname="K_true") %>% 
  dplyr::select(idd, N, K_true_cat, nmi, nmi_KM) %>% 
  
  tidyr::pivot_longer(cols=c("nmi","nmi_KM"), values_to="nmi", names_to="Method") %>% 
  dplyr::mutate(Method=dplyr::case_when(Method=="nmi" ~ "BASCULE",
                                        Method=="nmi_KM" ~ "KMeans")) %>% 
  
  dplyr::mutate(Method=reorder(Method, nmi, mean, decreasing=T)) %>% 
  
  ggplot(aes(x=factor(N), y=nmi, fill=Method, color=Method)) +
  stat_summary(aes(group=Method), position=position_dodge(width=0.2), fun.data="mean_cl_boot") +
  stat_summary(aes(group=Method), position=position_dodge(width=0.2), fun.data="mean_cl_boot",
               geom="line", linewidth=1) +
  
  scale_fill_manual(values=pal_methods, breaks=names(pal_methods)) +
  scale_color_manual(values=pal_methods, breaks=names(pal_methods)) +
  theme_bw() + ylim(NA, 1)
  


# Comparison ####

## precision and recall ####
plots[["recall"]] = stats_compare %>% dplyr::rowwise() %>%
  dplyr::mutate(FN=length(assigned_missing$missing_fn),
                FP=length(assigned_missing$added_fp),
                TP=length(assigned_missing$assigned_tp)) %>%
  dplyr::mutate(recall=TP / (TP + FN)) %>% 
  tidyr::pivot_longer(cols=c("recall"), names_to="prec_recall") %>% 
  dplyr::mutate(prec_recall=stringr::str_to_title(prec_recall)) %>% 
  
  dplyr::mutate(penalty=reorder(penalty, value, mean, decreasing=T)) %>% 
  
  ggplot(aes(x=factor(N), y=value, fill=penalty, color=penalty)) + 
  stat_summary(aes(group=penalty), position=position_dodge(width=.15), fun.data="mean_cl_boot") +
  stat_summary(aes(group=penalty), fun.data="mean_cl_boot", position=position_dodge(width=.15), 
               geom="line", linewidth=1) +
  
  theme_bw() + 
  scale_y_continuous(breaks=scales::pretty_breaks(n=3), limits=c(NA, 1)) +
  scale_fill_manual(values=pal_methods, breaks=names(pal_methods)) +
  scale_color_manual(values=pal_methods, breaks=names(pal_methods)) 


## mse ####
plots[["mse_counts"]] = stats_compare %>%
  dplyr::select(N, idd, mse_counts, penalty) %>% dplyr::rename(value=mse_counts) %>% 
  dplyr::mutate(metric="Counts") %>% 
  
  dplyr::mutate(penalty=reorder(penalty, value, mean, decreasing=T)) %>% 
  
  ggplot(aes(x=factor(N), y=value, fill=penalty, color=penalty)) + 
  stat_summary(aes(group=penalty), position=position_dodge(width=.15), fun.data="mean_cl_boot") +
  stat_summary(aes(group=penalty), position=position_dodge(width=.15), fun.data="mean_cl_boot", 
               geom="line", linewidth=1) +
  theme_bw() + # facet_grid(~metric) +
  scale_y_continuous(breaks=scales::pretty_breaks(n=3), limits=c(0, NA),
                     labels=function(x) scales::scientific(x)) +
  scale_fill_manual(values=pal_methods, breaks=names(pal_methods)) +
  scale_color_manual(values=pal_methods, breaks=names(pal_methods))


plots[["cosine_sigs"]] = stats_compare %>%
  dplyr::select(N, idd, cosine_sigs, penalty) %>% dplyr::rename(value=cosine_sigs) %>% 
  dplyr::mutate(metric="Signatures") %>% 
  
  dplyr::mutate(penalty=reorder(penalty, value, mean, decreasing=T)) %>% 
  
  ggplot(aes(x=factor(N), y=value, fill=penalty)) + 
  stat_summary(aes(group=penalty, color=penalty), position=position_dodge(width=.15),
               geom="pointrange", fun.data="mean_cl_boot") +
  stat_summary(aes(group=penalty, color=penalty), geom="line",
               position=position_dodge(width=.15), fun.data="mean_cl_boot", linewidth=1) +
  theme_bw() + # facet_grid(~metric) +
  # scale_y_continuous(n.breaks=5) +
  scale_y_continuous(breaks=scales::pretty_breaks(n=3), limits=c(NA, 1)) +
  scale_fill_manual(values=pal_methods, breaks=names(pal_methods)) +
  scale_color_manual(values=pal_methods, breaks=names(pal_methods)) 


plots[["cosine_expos"]] = stats_compare %>%
  dplyr::select(N, idd, cosine_expos_missing, penalty) %>% 
  dplyr::rename(value=cosine_expos_missing) %>% 
  dplyr::mutate(metric="Exposures") %>% 
  
  dplyr::mutate(penalty=reorder(penalty, value, mean, decreasing=T)) %>% 
  
  ggplot(aes(x=factor(N), y=value, fill=penalty, color=penalty)) + 
  stat_summary(aes(group=penalty), fun.data="mean_cl_boot", position=position_dodge(width=.2)) +
  stat_summary(aes(group=penalty), fun.data="mean_cl_boot", position=position_dodge(width=.2), 
               geom="line", linewidth=1) +
  
  theme_bw() + # facet_grid(~metric) +
  scale_y_continuous(breaks=scales::pretty_breaks(n=3), limits=c(NA, 1)) +
  scale_fill_manual(values=pal_methods, breaks = names(pal_methods)) +
  scale_color_manual(values=pal_methods, breaks = names(pal_methods)) 


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
times_bascule = read.csv("~/Dropbox/dropbox_shared/2022. Basilica/simulations/runtimes/bascule_exectimes.csv") %>% 
  dplyr::rename(execution_time_bascule=execution_time_SBS) %>% tibble::as_tibble() %>% 
  dplyr::mutate(tool="BASCULE", execution_time=execution_time_bascule) %>%
  dplyr::select(simulation_name, execution_time, execution_time_bascule, tool)

plots[["runtime"]] = dplyr::bind_rows(times_sigpr, 
                                      times_sparsesig, 
                                      times_bascule %>% dplyr::select(-execution_time_bascule)) %>% 
  dplyr::inner_join(times_bascule %>% dplyr::select(-execution_time, -tool)) %>% 

  dplyr::mutate(time_gain=execution_time / execution_time_bascule) %>% 
  
  dplyr::filter(tool!="BASCULE") %>% 
  
  dplyr::rowwise() %>% 
  dplyr::mutate(N=strsplit(simulation_name, "[.]")[[1]][2] %>% stringr::str_remove_all("N") %>% as.numeric()) %>% 
  
  dplyr::mutate(tool=reorder(tool, time_gain, mean, decreasing=T)) %>% 
  
  ggplot(aes(x=factor(N), y=time_gain, fill=tool, color=tool)) +
  stat_summary(aes(group=tool), position=position_dodge(width=.15), fun.data="mean_cl_boot") +
  stat_summary(aes(group=tool), fun.data="mean_cl_boot", position=position_dodge(width=.15), 
               geom="line", linewidth=1) +
  theme_bw() +
  scale_fill_manual(values=pal_methods, breaks = names(pal_methods)) +
  scale_color_manual(values=pal_methods, breaks = names(pal_methods)) + 
  scale_y_continuous(limits=c(1, NA)
                     # labels=function(x) round(10^x, digits = 0),
                     # breaks=scales::pretty_breaks(n=4), 
                     )


## example fit #####

# fit_simul = readRDS("~/Dropbox/dropbox_shared/2022. Basilica/simulations/fits/fits_dn.matched.2011/simul_fit.N500.G3.s11.matched.2011.Rds")
# fit_simul = readRDS("~/Dropbox/dropbox_shared/2022. Basilica/simulations/fits/fits_dn.matched.2011/simul_fit.N150.G3.s22.matched.2011.Rds")
fit_simul = readRDS("~/Dropbox/dropbox_shared/2022. Basilica/simulations/fits/fits_dn.matched.2011/simul_fit.N150.G3.s12.matched.2011.Rds")
fit_simul = readRDS("~/Dropbox/dropbox_shared/2022. Basilica/simulations/fits/fits_dn.matched.2011/simul_fit.N500.G3.s14.matched.2011.Rds")


bas_mapped = fit_simul$fit.0.auto %>% 
  convert_dn_names(reference_cat=get_signatures(fit_simul$dataset, matrix=T), cutoff=0.7) %>% 
  merge_clusters()

plot_exposures(bas_mapped)
plot_exposures(fit_simul$dataset)

# clusters_new = c("G3","G2","G1") %>% setNames(c("G0","G1","G3"))
# clusters_new = c("G1","G3","G2") %>% setNames(c("G3","G1","G2"))
# clusters_new = c("G1","G3","G2","Unmatched") %>% setNames(c("G0","G1","G2","G3"))
clusters_new = c("G3","G1","G2","Unmatched") %>% setNames(c("G0","G1","G2","G3"))

plots[["example"]] = get_exposure(bas_mapped, add_groups=T)[["SBS"]] %>% 
  dplyr::rowwise() %>% 
  dplyr::mutate(clusters=clusters_new[clusters], method="Predicted") %>% 
  dplyr::ungroup() %>% 
  
  dplyr::bind_rows(
    get_exposure(fit_simul$dataset, add_groups=T)[["SBS"]] %>% 
      dplyr::mutate(method="Ground truth") 
  ) %>% 
  
  dplyr::group_by(samples) %>%
  dplyr::mutate(clusters=replace(clusters, length(unique(clusters))>1, "UM")) %>%
  
  ggplot() +
  geom_bar(aes(x=samples, y=value, fill=sigs), stat="identity") +
  facet_grid(factor(method, levels=c("Ground truth","Predicted")) ~ clusters, scales="free_x", space="free_x") +
  scale_fill_manual(values=gen_palette(n=9)) +
  
  theme_bw() +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), 
        panel.grid.major.x=element_blank(), panel.grid.major.y=element_blank()) +
  labs(fill="Signatures")



# Panels #####
# panelsAB = ggplot()
panelsAB = plots[["example"]] + ylab("") +
  labs(title="Inference on a simulated dataset",
       subtitle="N=500, K=5, G=3") +
  theme(legend.position="bottom") + xlab("Samples") +
  labs(fill="Signatures")

panelC = plots[["recall"]] + ylab("Recall") +
  labs(title="Signatures detection accuracy") +
  theme(legend.position="none") + xlab("# samples") +
  guides(fill=guide_legend(title="Method"),
         color=guide_legend(title="Method"))

panelD = plots[["mse_counts"]] + 
  labs(title="Reconstruction error") +
  ylab("MSE") + xlab("# samples") +
  theme(legend.position="none") +
  guides(fill=guide_legend(title="Method"),
         color=guide_legend(title="Method"))

panelE = plots[["cosine_sigs"]] + 
  labs(title="Signatures quality") +
  ylab("Cosine similarity") + xlab("# samples") +
  theme(legend.position="none") +
  guides(fill=guide_legend(title="Method"),
         color=guide_legend(title="Method"))

panelF = plots[["cosine_expos"]] +
  labs(title="Exposures quality") +
  ylab("Cosine similarity") + xlab("# samples") +
  theme(legend.position="none") +
  guides(fill=guide_legend(title="Method"),
         color=guide_legend(title="Method"))

panelG = plots[["nmi"]] +
  labs(title="Clustering accuracy") +
  ylab("Normalized mutual information") + 
  xlab("# samples") + theme(legend.position="none") +
  guides(fill=guide_legend(title="Method"),
         color=guide_legend(title="Method"))

panelH = plots[["runtime"]] + 
  labs(title="Runtime increment") +
  ylab("Relative time increment") + 
  xlab("# samples") + theme(legend.position="none") +
  guides(fill=guide_legend(title="Method"),
         color=guide_legend(title="Method"))


(patchwork::wrap_plots(panelsAB, panelC, panelG,
                        panelD, panelE, 
                        panelF, panelH,
                        # guides="collect",
                        design="AAAABBCC\nDDEEFFGG")
    # theme(legend.position="bottom") &
  ) & patchwork::plot_annotation(tag_levels="A")

ggsave("~/Dropbox/dropbox_shared/2022. Basilica/paper/figure2/draft_fig2BIS.png", height=8, width=14)
ggsave("~/Dropbox/dropbox_shared/2022. Basilica/paper/figure2/draft_fig2BIS.pdf", height=8, width=14)






