devtools::load_all("~/GitHub/simbasilica/")
devtools::load_all("~/GitHub/basilica/")

# x_old = readRDS("~/Google Drive/My Drive/work/basilica_shared/matched.2011/fits_1102/fit_wcat_penalty0.Breast.Rds")
# x_old2 = x_old %>% refinement()

x = readRDS("~/Google Drive/My Drive/work/basilica_shared/fit_refined_clustered.Rds")

serena_cat = readRDS("~/Google Drive/My Drive/work/basilica_shared/processed_data/catalogues/signatures_sbs.Rds")
columns_order = colnames(get_signatures(x, matrix=T)[["SBS"]])
serena_private = setdiff(rownames(serena_cat), rownames(COSMIC_sbs_filt))

reference_sbs = rbind(COSMIC_sbs_filt[,columns_order], 
                      serena_cat[serena_private,columns_order])

get_assigned_missing(x, reference_cat=list("SBS"=reference_sbs), cutoff=0.8)


exposures = get_exposure(x, matrix=T) %>% do.call(cbind, .)
scale_exposures = scale(exposures)
pca = prcomp(scale_exposures, center=TRUE, scale.=TRUE)

pca$x[,1:2] %>% 
  ggplot() +
  geom_point(aes(x=PC1, y=PC2))





