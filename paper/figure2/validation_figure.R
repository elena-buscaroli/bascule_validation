library(magrittr)
library(ggplot2)
source("~/GitHub/basilica_validation/aux_fns/eval_aux_fns.R")
source("~/GitHub/basilica_validation/aux_fns/plots_aux_fns.R")

df_path = "~/Dropbox/dropbox_shared/2022. Basilica/simulations/stats_dataframes/"
stats_basilica = readRDS(paste0(df_path, "stats_matched.2011_KM.Rds")) %>% 
  compute_quantiles(colname="K_true")
stats_compare = readRDS(paste0(df_path, "stats_matched.2011.compare.Rds")) %>% 
  compute_quantiles(colname="K_true") %>% 
  dplyr::filter(type=="SBS")

plots = list()

pal_k_true = wesanderson::wes_palette("Cavalcanti1", 6, type="continuous")[c(1,4,3)]
pal_methods = c("#7fb3d5", "#FF8C00", "#8FBC8B", "#DB7093") %>% setNames(c("Basilica", "SigProfiler", "SparseSignatures","KMeans"))


# Basilica ####

plots[["nmi"]] = stats_basilica %>% 
  dplyr::filter(type=="SBS") %>% 
  compute_quantiles(colname="K_true") %>% 
  dplyr::select(idd, N, K_true_cat, nmi, nmi_KM) %>% 
  
  tidyr::pivot_longer(cols=c("nmi","nmi_KM"), values_to="nmi", names_to="Method") %>% 
  dplyr::mutate(Method=dplyr::case_when(Method=="nmi" ~ "Basilica",
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



# Panels #####
panelsAB = ggplot()

panelC = plots[["recall"]] + ylab("Recall") +
  labs(title="Signatures detection accuracy") +
  theme(legend.position="bottom") + xlab("# samples") 

panelD = plots[["mse_counts"]] + 
  labs(title="Reconstruction error") +
  ylab("MSE") + xlab("# samples") 

panelE = plots[["cosine_sigs"]] + 
  labs(title="Signatures quality") +
  ylab("Cosine similarity") + xlab("# samples") 

panelF = plots[["cosine_expos"]] +
  labs(title="Exposures quality") +
  ylab("Cosine similarity") + xlab("# samples") 

panelG = plots[["nmi"]] +
  labs(title="Clustering accuracy") +
  ylab("Normalized mutual information") + 
  xlab("# samples") 

panelH = plots[["runtime"]] + 
  labs(title="Runtime increment") +
  ylab("Relative time increment") + 
  xlab("# samples") 


(patchwork::wrap_plots(panelsAB, panelC, panelG,
                        panelD, panelE, 
                        panelF, panelH,
                        guides="collect",
                        design="AAAABBCC\nDDEEFFGG") & 
    theme(legend.position="bottom") &
    guides(fill=guide_legend(title="Method"),
           color=guide_legend(title="Method"))) +
  patchwork::plot_annotation(tag_levels="A")

ggsave("~/Dropbox/dropbox_shared/2022. Basilica/paper/figure2/draft_fig2.png", height=8, width=14)
ggsave("~/Dropbox/dropbox_shared/2022. Basilica/paper/figure2/draft_fig2.pdf", height=8, width=14)






