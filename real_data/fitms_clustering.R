library(ggplot2)
library(dplyr)
devtools::load_all("~/GitHub/bascule/")
source("real_data/fitms_clustering_utils.R")


lapply(c("Breast","Colorectal","Lung"), function(organ_name) {
  
  x_organ = readRDS(paste0("~/Google Drive/My Drive/work/bascule_shared/compiled_fits/", tolower(organ_name), "_data.Rds"))
  x_orig = merge_clusters(x_organ$x_after, cutoff=0.8)
  x = convert_dn_names(x_orig, reference_cat=list(SBS=COSMIC_sbs_filt, DBS=COSMIC_dbs), cutoff=0.8)
  class(x_orig) = "bascule_obj"
  class(x) = "bascule_obj"
  
  sample_names = get_samples(x)
  
  exposures_sbs = read.csv("~/Google Drive/My Drive/work/bascule_shared/zenodo/real_data_validation/data/Degasperi/SBS_v2.03/RefSig_SBS_Exposures_v2.03.tsv", sep="\t") %>% 
    filter(organ==organ_name) %>% select(-cohort, -organ) %>% 
    select(where(~ any(. != 0)))
  exposures_sbs = exposures_sbs[rowSums(exposures_sbs) > 0, ]
  
  exposures_dbs =  read.csv("~/Google Drive/My Drive/work/bascule_shared/zenodo/real_data_validation/data/Degasperi/DBS_v1.01/RefSig_DBS_Exposures_v1.01.tsv", sep="\t") %>% 
    filter(organ==organ_name) %>% select(-cohort, -organ) %>%
    select(where(~ any(. != 0)))
  exposures_dbs = exposures_dbs[rowSums(exposures_dbs) > 0, ]
  
  common_samples = intersect(sample_names, intersect(exposures_sbs %>% rownames(), exposures_dbs %>% rownames()))
  
  
  counts_sbs = readxl::read_xlsx("~/Dropbox/dropbox_shared/2022. Basilica/real_data/processed_data/science_supmat/SupplementaryTables.xlsx", 
                                 sheet="Table S7") %>% 
    filter(sample %in% common_samples) %>% 
    tibble::column_to_rownames(var="sample")
  counts_dbs = readxl::read_xlsx("~/Dropbox/dropbox_shared/2022. Basilica/real_data/processed_data/science_supmat/SupplementaryTables.xlsx", 
                                 sheet="Table S8") %>% 
    filter(sample %in% common_samples) %>% 
    tibble::column_to_rownames(var="sample")
  
  
  exposures_sbs = exposures_sbs[common_samples, ] / rowSums(counts_sbs[common_samples, ])
  exposures_sbs[exposures_sbs <= 0] = 1e-15
  exposures_sbs = exposures_sbs / rowSums(exposures_sbs)
  
  exposures_dbs = exposures_dbs[common_samples, ] / rowSums(counts_dbs[common_samples, ])
  exposures_dbs[exposures_dbs <= 0] = 1e-15
  exposures_dbs = exposures_dbs / rowSums(exposures_dbs)
  
  
  sigs_sbs = read.csv("~/Google Drive/My Drive/work/bascule_shared/zenodo/real_data_validation/data/Degasperi/SBS_v2.03/RefSig_SBS_v2.03.tsv", sep="\t") %>% 
    select(colnames(exposures_sbs))
  
  sigs_dbs = read.csv("~/Google Drive/My Drive/work/bascule_shared/zenodo/real_data_validation/data/Degasperi/DBS_v1.01/RefSig_DBS_v1.01.tsv", sep="\t") %>% 
    select(colnames(exposures_dbs))
  
  bascule_obj = list(
    input=list(
      SBS=list(counts=wide_to_long(counts_sbs[common_samples, ], what="counts")),
      DBS=list(counts=wide_to_long(counts_dbs[common_samples, ], what="counts"))
    ),
    nmf=list(
      SBS=list(
        exposure=wide_to_long(exposures_sbs, what="exposures"),
        beta_fixed=wide_to_long(t(sigs_sbs), what="beta"),
        beta_denovo=NULL
      ),
      DBS=list(
        exposure=wide_to_long(exposures_dbs, what="exposures"),
        beta_fixed=wide_to_long(t(sigs_dbs), what="beta"),
        beta_denovo=NULL
      )
    )
  )
  class(bascule_obj) = "bascule_obj"
  bascule_obj
  
  
  if (file.exists(paste0("real_data/fitms_clustering/fit_", organ_name, ".Rds"))) {
    fitms_clust_merg = readRDS(paste0("real_data/fitms_clustering/fit_", organ_name, ".Rds"))
  } else {
    reticulate::use_condaenv("bascule-env")
    py = reticulate::import_from_path("pybascule", "~/GitHub/pybascule/")
    fitms_clust = fit_clustering(bascule_obj, cluster=10,
                                 seed_list=c(19,255,18321,331),
                                 py=py,
                                 autoguide=TRUE)
    fitms_clust_merg = merge_clusters(fitms_clust)
    saveRDS(fitms_clust_merg, paste0("real_data/fitms_clustering/fit_", organ_name, ".Rds"))
  }
  
  
  col_palette = get_color_palette(get_cosmic(), get_degasperi(), sig_cls_organ_all)
  missing_sigs = setdiff(c(get_signames(x) %>% unlist(use.names=F),
                           get_signames(fitms_clust_merg) %>% unlist(use.names=F)),
                         unlist(lapply(col_palette, names), use.names=F) %>% unique())
  col_palette = c(col_palette[[tolower(organ_name)]], 
                  gen_palette_aux(signames=list("tmp"=missing_sigs)))
  sigs_order = c(gtools::mixedsort(names(col_palette)), "Other")
  
  exposures_f = lapply(get_types(fitms_clust_merg), function(t)
    get_exposure(fitms_clust_merg, types=get_types(fitms_clust_merg), 
                 samples=common_samples,
                 clusters=get_cluster_labels(fitms_clust_merg), add_groups=TRUE)[[t]] %>%
      dplyr::mutate(type=t)) %>%
    do.call(rbind, .) %>% unique()

  pl_exposures_f = exposures_f %>% 
    group_by(sigs) %>%
    mutate(sigs=replace(sigs, sum(value > 0.05) < 10, "Other")) %>%
    ungroup() %>%
    group_by(clusters) %>% 
    filter(length(unique(samples)) >= 20) %>% unique() %>% 
    ggplot() +
    geom_bar(aes(x=samples, y=value, fill=sigs), stat="identity") +
    facet_grid(type ~ clusters, scales="free", space="free") +
    scale_fill_manual(values=c(col_palette, "Other"="gainsboro"), 
                      breaks=sigs_order) +
    theme_bw() +
    theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), legend.position="right")
  
  exposures_b = lapply(get_types(x), function(t)
    get_exposure(x, types=get_types(x), 
                 samples=common_samples,
                 clusters=get_cluster_labels(x), add_groups=TRUE)[[t]] %>%
      dplyr::mutate(type=t)) %>%
    do.call(rbind, .)
  
  pl_exposures_b = exposures_b %>% 
    group_by(sigs) %>%
    mutate(sigs=replace(sigs, sum(value > 0.05) < 10, "Other")) %>%
    ungroup() %>%
    group_by(clusters) %>% 
    filter(length(unique(samples)) >= 20) %>% unique() %>% 
    ggplot() +
    geom_bar(aes(x=samples, y=value, fill=sigs), stat="identity") +
    facet_grid(type ~ clusters, scales="free", space="free") +
    scale_fill_manual(values=c(col_palette, "Other"="gainsboro"), 
                      breaks=sigs_order) +
    theme_bw() +
    theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), legend.position="right")
  
  pl_exposures = patchwork::wrap_plots(pl_exposures_b + labs(title="BASCULE exposures and clustering", tag="A"),
                                       pl_exposures_f + labs(title="FitMS exposures and clustering", tag="B"), 
                                       ncol=1) & 
    xlab("Samples") & ylab("Relative exposures") & theme_legend & theme_text &
    theme(legend.position="bottom") & guides(fill=guide_legend(title="Signatures", nrow=2))
  
  wide_f = exposures_f %>% filter(type=="SBS") %>% select(-type, -clusters) %>%
    tidyr::pivot_wider(names_from="sigs", values_from="value") %>% 
    tibble::column_to_rownames(var="samples")
  wide_b = exposures_b %>% filter(type=="SBS") %>% select(-type, -clusters) %>%
    tidyr::pivot_wider(names_from="sigs", values_from="value") %>% 
    tibble::column_to_rownames(var="samples")
  
  common_sigs = intersect(colnames(wide_f), colnames(wide_b))

  cosine_rowwise = function(mat1, mat2) {
    cs = matrix(0, nrow=nrow(mat1), ncol=nrow(mat2))
    
    for (i in 1:nrow(mat1)) {
      for (j in 1:nrow(mat2)) {
        if (i == j)
          cs[i, j] = lsa::cosine(as.matrix(mat1)[i,], as.matrix(mat2)[j,])
      }
    }
    rownames(cs) = rownames(mat1)
    colnames(cs) = rownames(mat1)
    return(cs)
  }
  
  cs = cosine_rowwise(wide_f[common_samples,common_sigs], wide_b[common_samples,common_sigs]) %>% 
    as.data.frame() %>% 
    tibble::rownames_to_column(var="samples1") %>% 
    tidyr::pivot_longer(cols=-"samples1", names_to="samples2") %>% 
    filter(samples1 == samples2) %>% 
    mutate(organ=organ_name)
  
  fitms_clust_merg$median_cs = cs$value %>% median
  
  saveRDS(cs, paste0("real_data/fitms_clustering/cs_", organ_name, ".Rds"))
  ggsave(paste0("real_data/fitms_clustering/plot_expos_", organ_name, ".png"), pl_exposures,
         height=18*1.5, width=21*1.5, units="cm", device=png)
  
})


# Cosine Similarity #####

organs = c("Breast", "Colorectal", "Lung")
pl_cs = readRDS("real_data/fitms_clustering/cs_Breast.Rds") %>% 
  bind_rows(readRDS("real_data/fitms_clustering/cs_Colorectal.Rds")) %>% 
  bind_rows(readRDS("real_data/fitms_clustering/cs_Lung.Rds")) %>% 
  
  ggplot() +
  geom_violin(aes(x=organ, y=value, fill=organ, color=organ), 
              alpha=0.5, width=0.5, draw_quantiles=0.5) +
  scale_color_manual(name="Organ", values=yarrr::piratepal(palette="nemo") %>% setNames(organs)) +
  scale_fill_manual(name="Organ", values=yarrr::piratepal(palette="nemo") %>% setNames(organs)) +
  ylab("Cosine similarity") +
  theme_bw() + theme_text + theme_legend + theme(axis.title.x=element_blank()) +
  theme(legend.position="bottom")

ggsave("real_data/fitms_clustering/plot_cs.png", pl_cs,
       height=8, width=10, units="cm", device=png)

            