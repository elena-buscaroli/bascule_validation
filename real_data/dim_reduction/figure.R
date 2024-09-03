devtools::load_all("~/GitHub/simbasilica/")
devtools::load_all("~/GitHub/basilica/")
library(Rtsne)
library(umap)
library(GGally)

# x_old = readRDS("~/Google Drive/My Drive/work/basilica_shared/matched.2011/fits_1102/fit_wcat_penalty0.Breast.Rds")
# x_old2 = x_old %>% refinement()

x = readRDS("~/Google Drive/My Drive/work/basilica_shared/fit_refined_clustered.Rds")

serena_cat = readRDS("~/Google Drive/My Drive/work/basilica_shared/processed_data/catalogues/signatures_sbs.Rds")
columns_order = colnames(get_signatures(x, matrix=T)[["SBS"]])
serena_private = setdiff(rownames(serena_cat), rownames(COSMIC_sbs_filt))

reference_sbs = rbind(COSMIC_sbs_filt[,columns_order], 
                      serena_cat[serena_private,columns_order])

convert_dn_names(x, reference_cat=list("SBS"=reference_sbs), cutoff=0.8)
get_assigned_missing(x=x, reference_cat=list("SBS"=reference_sbs), cutoff=0.8)


exposures = get_exposure(x, matrix=T) %>% do.call(cbind, .)
counts = get_input(x, matrix=T) %>% do.call(cbind, .)

make_dimreduction_plots = function(df) {
  pca = prcomp(df, center=TRUE, scale.=TRUE)

  set.seed(032)
  tsne = Rtsne(df, dims=3, initial_dims=100)
  
  set.seed(483)
  umap = umap(as.matrix(df), n_components=3, metric="euclidean")
  
  pca_plt = pca$x %>% 
    as.data.frame() %>% 
    tibble::rownames_to_column(var="samples") %>% 
    dplyr::full_join(get_cluster_assignments(x %>% merge_clusters(cutoff=0.9))) %>% 
    dplyr::select(-samples) %>% 
    ggplot() +
    geom_point(aes(x=PC1, y=PC2, color=clusters)) +
    facet_grid(~clusters) + theme_bw() +
    scale_color_manual(values=ggsci::pal_simpsons()(7)) +
    labs(title="PCA")
  # ggpairs(aes(color=clusters), columns=1:2)
  
  tsne_plt = tsne$Y %>% 
    as.data.frame() %>% 
    dplyr::mutate(samples=rownames(exposures)) %>% 
    dplyr::full_join(get_cluster_assignments(x %>% merge_clusters(cutoff=0.9))) %>% 
    dplyr::select(-samples) %>% 
    dplyr::rename(PC1=V1, PC2=V2) %>% 
    ggplot() +
    geom_point(aes(x=PC1, y=PC2, color=clusters)) +
    facet_grid(~clusters) + theme_bw() +
    scale_color_manual(values=ggsci::pal_simpsons()(7)) +
    labs(title="tSNE")
  # ggpairs(aes(color=clusters), columns=1:2)
  
  umap_plt = umap$layout %>% 
    as.data.frame() %>% 
    tibble::rownames_to_column(var="samples") %>%
    dplyr::full_join(get_cluster_assignments(x %>% merge_clusters(cutoff=0.9))) %>%
    dplyr::select(-samples) %>% 
    dplyr::rename(PC1=V1, PC2=V2) %>% 
    ggplot() +
    geom_point(aes(x=PC1, y=PC2, color=clusters)) +
    facet_grid(~clusters) + theme_bw() +
    scale_color_manual(values=ggsci::pal_simpsons()(7)) +
    labs(title="UMAP")
  # ggpairs(aes(color=clusters), columns=1:2)
  return(list("pca"=pca_plt, "tsne"=tsne_plt, "umap"=umap_plt))
}

plts_counts = make_dimreduction_plots(counts)
plts_expos = make_dimreduction_plots(exposures)

pdf("~/GitHub/basilica_validation/real_data/dimensionality_reduction.pdf", height=10, width=12)
patchwork::wrap_plots(plts_counts$pca, plts_counts$tsne, plts_counts$umap, ncol=1, guides="collect") &
  theme(legend.position="bottom", legend.direction="horizontal") &
  patchwork::plot_annotation(title="Dimensionality reduction on mutation counts")

patchwork::wrap_plots(plts_expos$pca, plts_expos$tsne, plts_expos$umap, ncol=1, guides="collect") &
  theme(legend.position="bottom", legend.direction="horizontal") &
  patchwork::plot_annotation(title="Dimensionality reduction on exposures")
dev.off()


