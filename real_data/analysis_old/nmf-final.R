
# ==============================================================================
# ==============================================================================
# ==============================================================================

# choosing the main path
linux_path <- "/home/azad/"
mac_path <- "/Users/azadsadr/"
path <- mac_path

devtools::load_all(paste(path, "Documents/packages/basilica", sep = ""))

x0 <- readRDS(file = paste(path, "Nextcloud/basilica/basilica_fit/fit_refined_clustered.Rds", sep = ""))
x <- basilica:::merge_clusters(x = x0, cutoff = 0.9)


# ==============================================================================
# ======================== MAPPING BASILICA WITH SEERENA =======================
# ==============================================================================
serena_reference_path <- paste(path, "Nextcloud/basilica/serena/SBS_v2.03/RefSig_SBS_v2.03.tsv", sep = "") # set serena reference file path
serena_exposure_path <- paste(path, "Nextcloud/basilica/serena/SBS_v2.03/RefSig_SBS_Exposures_v2.03.tsv", sep = "") # set serena exposure file path

dbs_serena_ref_path <- paste(path, "Nextcloud/basilica/serena/DBS_v1.01/RefSig_DBS_v1.01.tsv", sep = "") # set serena reference file path
dbs_serena_exp_path <- paste(path, "Nextcloud/basilica/serena/DBS_v1.01/RefSig_DBS_Exposures_v1.01.tsv", sep = "") # set serena exposure file path


serena_reference_path <- "~/Google Drive/My Drive/work/basilica_shared/processed_data/SBS_v2.03/RefSig_SBS_v2.03.tsv"
serena_exposure_path <- "~/Google Drive/My Drive/work/basilica_shared/processed_data/SBS_v2.03/RefSig_SBS_Exposures_v2.03.tsv"



organ_type <- "Breast" # set oragan type

# loading data
source(paste(path, "Nextcloud/basilica/scripts/nmf-analysis/utils_nmf/utils_load.R", sep = "")) # load.all
obj <- load.all(
  basilica_obj = x,
  serena_ref_path = serena_reference_path,
  serena_exp_path = serena_exposure_path,
  organ = organ_type,
  reference = basilica::COSMIC_sbs_filt,
  type = "SBS"
)

# mapping basilica and serena inference
source(paste(path, "Nextcloud/basilica/scripts/nmf-analysis/utils_nmf/utils_map.R", sep = "")) # map
map <- map.data(
  reference = basilica::COSMIC_sbs_filt,
  fixed = obj$fixed,
  denovo = obj$denovo,
  common = obj$common,
  rare = obj$rare,
  threshold = 0.8
)

rare_long = rare %>% wide_to_long(what="beta") %>% reformat_contexts(what="SBS") %>% dplyr::mutate(type="SBS")
common_long = common %>% wide_to_long(what="beta") %>% reformat_contexts(what="SBS") %>% dplyr::mutate(type="SBS")
reference_long = reference %>% wide_to_long(what="beta") %>% reformat_contexts(what="SBS") %>% dplyr::mutate(type="SBS")

plot_signatures_aux(common_long %>% dplyr::filter(sigs%in%map$denovo_common$second)) %>%
  patchwork::wrap_plots(plot_signatures(x, signames=map$denovo_common$first))

plot_signatures_aux(rare_long %>% dplyr::filter(sigs%in%map$denovo_rare$second)) %>%
  patchwork::wrap_plots(plot_signatures(x, signames=map$denovo_rare$first))

plot_signatures_aux(reference_long %>% dplyr::filter(sigs%in%map$denovo_reference$second)) %>%
  patchwork::wrap_plots(plot_signatures(x, signames=map$denovo_reference$first))

plot_signatures_aux(common_long %>% dplyr::filter(sigs%in%map$fixed_common$second)) %>%
  patchwork::wrap_plots(plot_signatures(x, signames=map$fixed_common$first))

plot_signatures_aux(rbind(common_long, rare_long) %>% dplyr::filter(sigs%in%map$basilica_serena$Serena)) %>%
  patchwork::wrap_plots(plot_signatures(x, signames=map$basilica_serena$Basilica))



# visualization
source(paste(path, "Nextcloud/basilica/scripts/nmf-analysis/utils_nmf/utils_visual.R", sep = ""))
p1 <- plot_unexplained(
  basilica_exposure = obj$BasilicaExposure,
  serena_exposure = obj$SerenaExposure,
  basilica_singles = map$basilica_singles,
  serena_singles = map$serena_singles
)

p2 <- plot_exposure_comp(
  basilica_exposure = obj$BasilicaExposure,
  serena_exposure = obj$SerenaExposure,
  basilica_serena = map$basilica_serena
)

#map$denovo_reference
#map$fixed_reference
#map$common_reference
#map$rare_reference

# serena | basilica | reference
# SBS1   |     -    | SBS1
# SBS3   |     -    |    -
# SBS127 |     -    |    -
#   -    |  SBSD7   |    -
#   -    |  SBSD8   |    -
#   -    |  SBSD17  | SBS90

# ==============================================================================
# ================ VISUALIZE UN-EXPLAINED SIGNATURES (BASILICA) ================
# ==============================================================================
# plot the signatures
source(paste(path, "Nextcloud/basilica/scripts/visualization/utils_plot.R", sep = ""))
source(paste(path, "Nextcloud/basilica/scripts/visualization/plot_signatures.R", sep = "")) # plot_signatures
SBSD7 <- plot_signatures(
  SBS = basilica:::wide_to_long(obj$denovo, what = "beta") %>% subset(sigs == "SBSD7"),
  DBS = NULL, types = c("SBS")
)
SBSD8 <- plot_signatures(
  SBS = basilica:::wide_to_long(obj$denovo, what = "beta") %>% subset(sigs == "SBSD8"),
  DBS = NULL, types = c("SBS")
)
SBSD17 <- plot_signatures(
  SBS = basilica:::wide_to_long(obj$denovo, what = "beta") %>% subset(sigs == "SBSD17"),
  DBS = NULL, types = c("SBS")
)
SBS90 <- plot_signatures(
  SBS = basilica:::wide_to_long(basilica::COSMIC_filt, what = "beta") %>% subset(sigs == "SBS90"),
  DBS = NULL, types = c("SBS")
)

# ==============================================================================
# ================ VISUALIZE UN-EXPLAINED SIGNATURES (SERENA) ==================
# ==============================================================================
SBS1_C <- plot_signatures(
  SBS = basilica:::wide_to_long(obj$common, what = "beta") %>% subset(sigs == "SBS1"),
  DBS = NULL, types = c("SBS")
)
SBS3_C <- plot_signatures(
  SBS = basilica:::wide_to_long(obj$common, what = "beta") %>% subset(sigs == "SBS3"),
  DBS = NULL, types = c("SBS")
)
SBS127_C <- plot_signatures(
  SBS = basilica:::wide_to_long(obj$common, what = "beta") %>% subset(sigs == "SBS127"),
  DBS = NULL,
  types = c("SBS")
)
SBS1 <- plot_signatures(
  SBS = basilica:::wide_to_long(basilica::COSMIC_filt, what = "beta") %>% subset(sigs == "SBS1"),
  DBS = NULL, types = c("SBS")
)

#===============================================================================

p1
p2

SBS1_C
SBS1

SBS3_C
SBS127_C

SBSD7
SBSD8

SBSD17
SBS90

#===============================================================================

plot_signatures(
  SBS = basilica:::wide_to_long(obj$common, what = "beta") %>% subset(sigs == "SBS3"),
  DBS = NULL, types = c("SBS")
)

plot_signatures(
  SBS = basilica:::wide_to_long(obj$rare, what = "beta") %>% subset(sigs == "SBS96"),
  DBS = NULL, types = c("SBS")
)

plot_signatures(
  SBS = basilica:::wide_to_long(obj$fixed, what = "beta") %>% subset(sigs == "SBS3"),
  DBS = NULL, types = c("SBS")
)

plot_signatures(
  SBS = basilica:::wide_to_long(obj$denovo, what = "beta") %>% subset(sigs == "SBS3"),
  DBS = NULL, types = c("SBS")
)

plot_signatures(
  SBS = basilica:::wide_to_long(basilica::COSMIC_filt, what = "beta") %>% subset(sigs == "SBS3"),
  DBS = NULL, types = c("SBS")
)


cosine.matrix(obj$common["SBS127", ], basilica::COSMIC_filt)

plot_exposure_dist(obj$BasilicaExposure, obj$SerenaExposure, map)



