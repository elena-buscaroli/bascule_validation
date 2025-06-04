library(ggplot2)
library(dplyr)
devtools::load_all("~/GitHub/bascule/")

for (organ_name in c("Breast","Colorectal","Lung")) {

  x_organ = readRDS(paste0("~/Google Drive/My Drive/work/bascule_shared/compiled_fits/", tolower(organ_name), "_data.Rds"))
  x_orig = merge_clusters(x_organ$x_after, cutoff=0.8)
  x = convert_dn_names(x_orig, reference_cat=list(SBS=COSMIC_sbs_filt, DBS=COSMIC_dbs), cutoff=0.8)
  class(x_orig) = "bascule_obj"
  class(x) = "bascule_obj"
  
  
  exposures_sbs = read.csv("~/Google Drive/My Drive/work/bascule_shared/zenodo/real_data_validation/data/Degasperi/SBS_v2.03/RefSig_SBS_Exposures_v2.03.tsv", sep="\t") %>% 
    filter(organ==organ_name) %>% select(-cohort, -organ) %>% 
    select(where(~ any(. != 0)))
  exposures_sbs = exposures_sbs[rowSums(exposures_sbs) > 0, ]
  # exposures_sbs = exposures_sbs %>% mutate_if(is.numeric, function(x) replace(x, x==0, 1e-10))
  
  exposures_dbs =  read.csv("~/Google Drive/My Drive/work/bascule_shared/zenodo/real_data_validation/data/Degasperi/DBS_v1.01/RefSig_DBS_Exposures_v1.01.tsv", sep="\t") %>% 
    filter(organ==organ_name) %>% select(-cohort, -organ) %>%
    select(where(~ any(. != 0)))
  exposures_dbs = exposures_dbs[rowSums(exposures_dbs) > 0, ]
  # exposures_dbs = exposures_dbs %>% mutate_if(is.numeric, function(x) replace(x, x==0, 1e-10))
  
  common_samples = intersect(exposures_sbs %>% rownames(), exposures_dbs %>% rownames())
  
  exposures_sbs %>% colSums() %>% min()
  exposures_dbs %>% colSums() %>% min()
  
  
  counts_sbs = readxl::read_xlsx("~/Dropbox/dropbox_shared/2022. Basilica/real_data/processed_data/science_supmat/SupplementaryTables.xlsx", 
                                 sheet="Table S7") %>% 
    filter(sample %in% common_samples) %>% 
    tibble::column_to_rownames(var="sample")
  counts_dbs = readxl::read_xlsx("~/Dropbox/dropbox_shared/2022. Basilica/real_data/processed_data/science_supmat/SupplementaryTables.xlsx", 
                                 sheet="Table S8") %>% 
    filter(sample %in% common_samples) %>% 
    tibble::column_to_rownames(var="sample")
  
  
  exposures_sbs = exposures_sbs[common_samples, ] / rowSums(counts_sbs[common_samples, ])
  exposures_sbs[exposures_sbs <= 0] = 1e-10
  exposures_sbs = exposures_sbs / rowSums(exposures_sbs)
  
  exposures_dbs = exposures_dbs[common_samples, ] / rowSums(counts_dbs[common_samples, ])
  exposures_dbs[exposures_dbs <= 0] = 1e-10
  exposures_dbs = exposures_dbs / rowSums(exposures_dbs)
  
  
  
  # with_unass = read.csv("~/Google Drive/My Drive/work/bascule_shared/zenodo/real_data_validation/data/Degasperi/DBS_v1.01/organSpecificExposures/GEL/GEL-Breast_DBS_exposures_finalT.tsv", sep="\t") %>% 
  #   rowSums()
  # 
  # unass = read.csv("~/Google Drive/My Drive/work/bascule_shared/zenodo/real_data_validation/data/Degasperi/DBS_v1.01/organSpecificExposures/GEL/GEL-Breast_DBS_exposures_finalT.tsv", sep="\t") %>% 
  #   pull(unassigned)
  # 
  # without_unass = read.csv("~/Google Drive/My Drive/work/bascule_shared/zenodo/real_data_validation/data/Degasperi/DBS_v1.01/organSpecificExposures/GEL/GEL-Breast_DBS_exposures_finalT.tsv", sep="\t") %>% 
  #   select(-unassigned) %>% rowSums()
  # 
  # unass[1:10] / with_unass[1:10] * 100
  
  
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
  
  
  reticulate::use_condaenv("bascule-env")
  py = reticulate::import_from_path("pybascule", "~/GitHub/pybascule/")
  fitms_clust = fit_clustering(bascule_obj, cluster=10,
                               seed_list=c(19,255,18321,331),
                               py=py,
                               autoguide=TRUE)
  fitms_clust_merg = merge_clusters(fitms_clust)
  
  col_palette = gen_palette(fitms_clust_merg)
  
  exposures = lapply(get_types(fitms_clust_merg), function(t)
    get_exposure(fitms_clust_merg, types=get_types(fitms_clust_merg), 
                 samples=get_samples(fitms_clust_merg),
                 clusters=get_cluster_labels(fitms_clust_merg), add_groups=TRUE)[[t]] %>%
      dplyr::mutate(type=t)) %>%
    do.call(rbind, .)
  
  pl_exposures = ggplot(exposures) +
    geom_bar(aes(x=samples, y=value, fill=sigs), stat="identity") +
    facet_grid(type ~ clusters, scales="free", space="free") +
    scale_fill_manual(values=col_palette, breaks=names(col_palette)) +
    theme_bw() +
    theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), legend.position="right")
  
  
  saveRDS(fitms_clust_merg, paste0("real_data/fitms_clustering/fit_", organ_name, ".Rds"))
  ggsave(paste0("real_data/fitms_clustering/plot_expos_", organ_name, ".png"), pl_exposures,
         height=10, width=30, units="cm", device=png)
  
}

