devtools::load_all("~/GitHub/simbasilica/")
devtools::load_all("~/GitHub/basilica/")
library(Rtsne)
library(umap)

# x_old = readRDS("~/Google Drive/My Drive/work/basilica_shared/matched.2011/fits_1102/fit_wcat_penalty0.Breast.Rds")
# x_old2 = x_old %>% refinement()

x = readRDS("~/Google Drive/My Drive/work/basilica_shared/fit_refined_clustered.Rds")

serena_cat = readRDS("~/Google Drive/My Drive/work/basilica_shared/processed_data/catalogues/signatures_sbs.Rds")
columns_order = colnames(get_signatures(x, matrix=T)[["SBS"]])
serena_private = setdiff(rownames(serena_cat), rownames(COSMIC_sbs_filt))

reference_sbs = rbind(COSMIC_sbs_filt[,columns_order], 
                      serena_cat[serena_private,columns_order])

convert_dn_names(x, reference_cat=list("SBS"=reference_sbs), cutoff=0.8)
get_assigned_missing(x, reference_cat=list("SBS"=reference_sbs), cutoff=0.8)


exposures = get_exposure(x, matrix=T) %>% do.call(cbind, .)
counts = get_input(x, matrix=T) %>% do.call(cbind, .)
n_muts = lapply(get_types(x), function(tid) get_exposure(x, matrix=T)[[tid]] * rowSums(get_input(x, matrix=T)[[tid]])) %>% 
  do.call(cbind, .)

pca = prcomp(exposures, center=TRUE, scale.=TRUE)
plot(pca)

corr_counts = scale(counts)
pca = prcomp(corr_counts, center=T, scale.=T)
plot(pca)

set.seed(6543)
tsne = Rtsne(counts, dims=3, initial_dims=100)

set.seed(6543)
umap = umap(as.matrix(counts), n_components=3, metric="euclidean")


pca$x %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column(var="samples") %>% 
  dplyr::full_join(get_cluster_assignments(x %>% merge_clusters(cutoff=0.9))) %>% 
  dplyr::select(-samples) %>% 
  ggpairs(aes(color=clusters), columns=1:2)


tsne$Y %>% 
  as.data.frame() %>% 
  dplyr::mutate(samples=rownames(exposures)) %>% 
  dplyr::full_join(get_cluster_assignments(x %>% merge_clusters(cutoff=0.9))) %>% 
  dplyr::select(-samples) %>% 
  ggplot() +
  geom_point(aes(x=V1, y=V2, color=clusters)) +
  facet_grid(~clusters)
  ggpairs(aes(color=clusters), columns=1:2)


umap$layout %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column(var="samples") %>%
  dplyr::full_join(get_cluster_assignments(x %>% merge_clusters(cutoff=0.9))) %>%
  dplyr::select(-samples) %>% 
  ggpairs(aes(color=clusters), columns=1:2)
