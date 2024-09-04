# organ_type = "Breast" # Breast, Lung, Colorectal
# path = "/Users/azadsadr/"
# devtools::load_all( paste(path, "Documents/packages/basilica/", sep="") )
# source(paste(path, "Nextcloud/basilica/scripts/utils_plot.R", sep=""))
# source(paste(path, "Nextcloud/basilica/scripts/utils/plot.signatures.R", sep=""))
# 
# data = readRDS( paste(path, "Nextcloud/basilica/compiled_fit/", tolower(organ_type), "_data.Rds", sep="") )

library(tidyverse)

organ_type = "Breast" # Breast, Lung, Colorectal
scripts_path = "~/Dropbox/dropbox_shared/2022. Basilica/"
devtools::load_all("~/GitHub/bascule/")
source("real_data/analysis/utils_plot.R")
# source(paste(path, "~/Google Drive/My Drive/work/basilica_shared//basilica/scripts/utils/plot.signatures.R", sep=""))

data = readRDS(paste0("~/Google Drive/My Drive/work/basilica_shared/compiled_fits/", tolower(organ_type), "_data.Rds", sep="") )
x_orig = merge_clusters(data$x_after, cutoff=0.8)
x = convert_dn_names(x_orig, reference_cat=list(SBS=COSMIC_sbs_filt, DBS=COSMIC_dbs), cutoff=0.8)
class(x_orig) = "bascule_obj"
class(x) = "bascule_obj"


sig_cls_organ = sig_cls_organ_all[[tolower(organ_type)]]
clusters_to_remove = c()
cluster_names = setdiff(gtools::mixedsort(names(sig_cls_organ)),
                        clusters_to_remove)
sig_cls_organ = sig_cls_organ[cluster_names]


theme_legend = theme(legend.text=element_text(size=7.4),
                     legend.title=element_text(size=10.6),
                     legend.key.size=unit(0.2,"cm"),
                     legend.key.height=unit(0.2,"cm"),
                     legend.key.width=unit(0.2,"cm"))

theme_text = theme(axis.title=element_text(size=7.4),
                   axis.text=element_text(size=5.35), 
                   plot.title=element_text(size=10.6),
                   plot.subtitle=element_text(size=7.4),
                   strip.text=element_text(size=5.35), 
                   plot.caption=element_text(size=5.35))

theme_noY = theme(axis.ticks.y=element_blank(),
                  axis.text.y=element_blank(),
                  axis.line.y=element_blank(),
                  axis.ticks.length.y=unit(0,"pt"))

theme_noX = theme(axis.ticks.x=element_blank(),
                  axis.text.x=element_blank(),
                  axis.line.x=element_blank(),
                  axis.ticks.length.x=unit(0,"pt"))

plots = list()

# MODEL SELECTION #####
plots[["scores"]] = plot_scores_nmf(x, types=get_types(x), remove_outliers=FALSE) +
  theme(panel.grid.minor=element_blank()) +
  theme_text


# MAPPING ####
cls_catalogues = c("chocolate", "#5499c7", "darkseagreen") %>% 
  setNames(c("De novo","COSMIC v3.4","Degasperi et al."))
plots[["mapped_sigs"]] = plot_mapping_barplot(x_orig, cls_catalogues) +
  theme_text + theme_legend


# CENTROIDS + AETIOLOGY #####

aetiology = read.csv("~/Google Drive/My Drive/work/basilica_shared/codes/sbs_aetiology.csv") %>% 
  dplyr::bind_rows(read.csv("~/Google Drive/My Drive/work/basilica_shared/codes/dbs_aetiology.csv"))

cls = get_color_palette(cosmic, degasperi, sig_cls_organ_all)[[tolower(organ_type)]]

plots_tmp = custom_centroid_plot(
  x=x, 
  clusters=cluster_names,
  sig_cls=sig_cls_organ, 
  col_palette=cls,
  sbs_aetiology_path="~/Google Drive/My Drive/work/basilica_shared/codes/sbs_aetiology.csv", 
  dbs_aetiology_path="~/Google Drive/My Drive/work/basilica_shared/codes/dbs_aetiology.csv"
) %>% lapply(function(i) i + theme_text)

plots_tmp$heatmap = plots_tmp$heatmap +
  theme(axis.text.y.left=element_text(colour="white"),
        axis.title.y.left=element_text(colour="white"), 
        strip.background=element_rect(color="white", fill="white"),
        strip.text=element_text(color="white"),
        strip.placement="inside",
        legend.position="bottom") + theme_legend
# guides(fill=guide_legend(title="Signatures", position="none")) +
# theme_legend

plots[["aetiology"]] = patchwork::wrap_plots(
  plots_tmp$centroids + theme_legend + theme_noX + theme(axis.title.x=element_blank()), 
  # plots_tmp$heatmap,
  patchwork::wrap_elements(full=plots_tmp$heatmap, ignore_tag=TRUE),
  ncol=1, heights=c(4,1), guides="collect") + 
  patchwork::plot_layout(axes="collect_x")

plots[["aetiology"]]


# DENOVO SIGNATURES [SBS + DBS] #####

aetiology %>% dplyr::filter(signature %in% (dplyr::bind_rows(assigned) %>% 
                                         dplyr::filter(signature_type=="De novo") %>%
                                         dplyr::pull(reference))) %>% 
  dplyr::filter(aetiology!="UNKNOWN")
denovo_catalogue = get_denovo_signatures(x_orig)

dn1 = "SBSD10"
dn2 = "SBSD11"
ref1 = assigned$SBS %>% dplyr::filter(bascule==dn1) %>% dplyr::pull(reference)
ref2 = assigned$SBS %>% dplyr::filter(bascule==dn2) %>% dplyr::pull(reference)
cat1 = assigned$SBS %>% dplyr::filter(bascule==dn1) %>% dplyr::pull(catalogue)
cat2 = assigned$SBS %>% dplyr::filter(bascule==dn2) %>% dplyr::pull(catalogue)

fn_tmp = function(cat_name) if (cat_name=="COSMIC v3.4") return(cosmic) else return(degasperi)
sigs1 = plot_mirrored_sigs(denovo_catalogue[["SBS"]], fn_tmp(cat1), dn1, ref1, cat1, type="SBS", thr=0)
sigs2 = plot_mirrored_sigs(denovo_catalogue[["SBS"]], fn_tmp(cat2), dn2, ref2, cat2, type="SBS", thr=0)

plots[["denovo_sig1"]] = sigs1 + theme(legend.position="bottom") + theme_text + theme_legend
plots[["denovo_sig2"]] = sigs2 + theme(legend.position="bottom") + theme_text + theme_legend


# plots[["denovo_sigs"]] = patchwork::wrap_plots(sigs1, sigs2, 
#                                                guides="collect", 
#                                                nrow=1) & 
#   theme(legend.position="bottom") & theme_text & theme_legend &
#   theme_noX & theme_noY 


# EXPOSURES #####

cls["Other"] = "gainsboro"
cls = cls[!is.na(names(cls))]
sigs_order = c(gtools::mixedsort(unique(unlist(sig_cls_organ))), "Other")

n_samples = sapply(cluster_names, function(cl_id) get_cluster_assignments(x, clusters=cl_id) %>% nrow())
plots[["exposures"]] = lapply(cluster_names, function(cl_id) {
  pl = plot_exposures(x, signatures_list=sigs_order, 
                      clusters=cl_id, color_palette=cls) +
    scale_fill_manual(values=cls, breaks=sigs_order, limits=sigs_order) +
    scale_y_continuous(breaks=c(0,1)) +
    xlab("Samples") + ylab("Relative exposures") +
    theme(panel.grid=element_blank()) +
    xlab(paste0("n=", n_samples[cl_id]))
    # labs(caption=paste0("N=", n_samples[cl_id]))

  if (cl_id != cluster_names[1]) 
    pl = pl + theme(plot.margin=margin(t=5.5,r=1,l=1,b=5.5,unit="pt"),
                    axis.title.y=element_blank(),
                    axis.ticks.length.y=unit(0,"pt")) +
      theme_noY
  
  if (cl_id != cluster_names[length(cluster_names)])
    pl = pl + theme(strip.background.y=element_blank(),
                    strip.text.y=element_blank(),
                    plot.margin=margin(t=5.5,r=1,l=1,b=5.5,unit="pt"))
  return(pl)
}) %>% patchwork::wrap_plots(nrow=1, widths=n_samples, guides="collect", tag_level="new")

plots[["exposures"]] = plots[["exposures"]] + 
  patchwork::plot_layout(axis_titles="collect", axes="collect") &
  theme_text & theme_legend & theme(plot.caption=element_text(hjust=.5))

# FINAL PLOT #####

panelA = plots$scores + 
  labs(title="Model selection", subtitle="BIC of the tested bNMF models")

panelB = plots$mapped_sigs + 
  labs(title="Mapping of de novo signatures",
       subtitle="De novo signatures mapped to known catalogues") +
  theme(axis.title.y=element_blank())

plots$aetiology[[1]] = plots$aetiology[[1]] + 
  labs(title="Relevant signatures in clusters",
       subtitle="Prevalence and aetiology of relevant signatures in clusters") +
  guides(fill=guide_legend(ncol=1, title="Signatures", override.aes=list(size=0.5)))
panelC = plots$aetiology

panelD = plots$denovo_sig1 + 
  theme(axis.text.x=element_text(angle=90, size=3), axis.text.y=element_text(),
        legend.position="right") +
  scale_y_continuous(label=function(x) abs(x))
panelE = plots$denovo_sig2 + 
  theme(axis.text.x=element_text(angle=90, size=3), axis.text.y=element_text(),
        legend.position="right") +
  scale_y_continuous(label=function(x) abs(x))

plots$exposures[[1]] = plots$exposures[[1]] + 
  labs(title="Samples exposures", tag="F",
       subtitle="Exposures to relevant signatures across all samples")
plots$exposures[[2]] = plots$exposures[[2]] + theme(plot.title=element_blank(), plot.background=element_blank())
panelF = plots$exposures & 
  theme(legend.position="bottom") & 
  guides(fill=guide_legend(nrow=2, title="Signatures", override.aes=list(size=0.5)))


a = patchwork::wrap_plots(
  panelA,
  panelB,
  panelC,
  panelD,
  panelE,
  panelF,
  design="ABC\nABC\nDDD\nEEE\nFFF\nFFF"
) + patchwork::plot_annotation(tag_levels="A")


ggsave(paste0("paper/figure3/", tolower(organ_type), ".pdf"), plot=a, width=210, height=297, units="mm")

