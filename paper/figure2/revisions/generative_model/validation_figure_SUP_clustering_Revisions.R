library(magrittr)
library(dplyr)
library(ggplot2)
source("~/GitHub/bascule_validation/synthetic_data/aux_fns/eval_aux_fns.R")
source("~/GitHub/bascule_validation/synthetic_data/aux_fns/plots_aux_fns.R")

df_path = "~/Dropbox/dropbox_shared/2022. Basilica/simulations/stats_dataframes/revisions/"
stats_bascule = readRDS(paste0(df_path, "stats_matched.2011.compare_LAST.generative_model.nofits.Rds")) %>% 
  compute_quantiles(colname="K_true") %>% dplyr::filter(penalty=="BASCULE")

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

pal_methods = c("#7fb3d5", "#FF8C00", "#8FBC8B", "#DB7093", RColorBrewer::brewer.pal(4, name="Dark2")) %>%
  setNames(c("BASCULE", "SigProfiler", "SparseSignatures","FitMS_E","KMeans","KL-KMeans","JS-Spectral","FitMS"))


# Clustering ####

input_df = stats_bascule %>% 
  dplyr::filter(type=="SBS") %>% 
  
  dplyr::select(idd, N, G, source, K_true_cat, starts_with("nmi")) %>%
  
  tidyr::pivot_longer(cols=starts_with("nmi"), values_to="nmi", names_to="Method") %>%
  dplyr::mutate(Method=dplyr::case_when(Method=="nmi" ~ "BASCULE",
                                        Method=="nmi_KM" ~ "KMeans",
                                        Method=="nmi_KL" ~ "KL-KMeans",
                                        Method=="nmi_JS" ~ "JS-Spectral")) %>%
  
  dplyr::mutate(Method=reorder(Method, nmi, mean, decreasing=T),
                N=factor(N, levels=unique(sort(N))))

pvalues = ggpubr::compare_means(nmi ~ Method, data=input_df, method="wilcox.test", 
                                exact=FALSE, group.by=c("N","K_true_cat")) %>%
  filter(group1=="BASCULE" | group2=="BASCULE") %>% 
  group_by(N, K_true_cat) %>% 
  mutate(y.position=0.95 - 0.05 * row_number()) %>% 
  ungroup() %>% 
  mutate(p=scales::scientific(p, digits=2))

clustering = input_df %>% 
  
  ggplot(aes(x=factor(N), y=nmi)) +
  geom_boxplot(aes(fill=Method, color=Method), outlier.shape=NA, lwd=.5, alpha=.7, width=0.5) +
  ggpubr::stat_pvalue_manual(
    pvalues %>% ungroup(), label="p.signif", x="N", color="group2",
    tip.length=0.01, remove.bracket=F,
    bracket.size=0.4, position=position_nudge(y=0.25)
  ) +
  ggh4x::facet_nested(~ "# signatures" + K_true_cat) +
  scale_fill_manual(values=pal_methods,
                    breaks=names(pal_methods),
                    name="Method") +
  scale_color_manual(values=pal_methods,
                     breaks=names(pal_methods),
                     name="Method") +
  theme_bw() + ylim(0, NA)
  
clustering


figure = clustering + xlab("# samples") + ylab("NMI") +
  labs(title="Clustering accuracy",
       subtitle="NMI between inferred and true assignments") +
  theme_text + theme_legend + theme(legend.position="bottom")

figure

ggsave(filename="paper/figure2/revisions/generative_model/figure2_SUP3.png", plot=figure,
       width=210, height=120, units="mm", device=png, family="Helvetica")

