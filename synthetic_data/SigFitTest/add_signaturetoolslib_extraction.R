library(tidyverse)
devtools::load_all("~/GitHub/bascule/")

create_bascule_obj = function(counts, expos, sigs) {
  obj = list(); class(obj) = "bascule_obj"
  obj[["input"]][["SBS"]] = list("counts"=counts,
                                 "reference"=NULL)
  obj[["nmf"]][["SBS"]] = list("exposure"=wide_to_long(expos, what="exposures"),
                               "beta_fixed"=NULL,
                               "beta_denovo"=wide_to_long(sigs, what="beta"))
  return(obj)
}


main_path = "~/Dropbox/dropbox_shared/2022. Basilica/simulations/fits_SigFitTest//"
fits_path = file.path(main_path, "signaturetoolslib_E/")

save_path = file.path(main_path, "fits_dn.matched.2011.compare_LAST/")

fitsnames = list.files(fits_path, pattern="simul_fit")

lapply(fitsnames, function(fname) {
  # if (file.exists(paste0(save_path, fname))) return()
  
  print(fname)
  old_compare_fit = readRDS(paste0(save_path, fname))
  
  counts = get_input(old_compare_fit$dataset)[["SBS"]] %>% dplyr::select(-clusters)
  
  stl_res = readRDS(paste0(fits_path, fname))
  stl_alpha = stl_res$exposures
  stl_alpha = stl_alpha / rowSums(stl_alpha)
  stl_sigs = stl_res$signatures %>% t() %>% as.data.frame()
  x.stlE = create_bascule_obj(counts, stl_alpha, stl_sigs)
  
  old_compare_fit$signaturetoolslib_E = x.stlE
  old_compare_fit$fit.0 = NULL
  
  saveRDS(old_compare_fit, paste0(save_path, fname))
})




