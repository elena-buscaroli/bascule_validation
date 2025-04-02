library(tidyverse)
devtools::load_all("~/GitHub/bascule/")


main_path = "~/Dropbox/dropbox_shared/2022. Basilica/simulations/fits_generative_model/"
fits_path = file.path(main_path, "all_fits/fits_dn.clustering.matched.2011/")

save_path = file.path(main_path, "all_fits/fits_dn.matched.2011.compare_LAST/")

fitsnames = list.files(fits_path, pattern="simul_fit")

lapply(fitsnames, function(fname) {
  # if (file.exists(paste0(save_path, fname))) return()
  
  print(fname)
  
  simul_fit = readRDS(paste0(fits_path, fname))
  old_compare_fit = readRDS(paste0(save_path, fname))
  
  old_compare_fit$x.fit0.auto = simul_fit$x.fit0.auto
  old_compare_fit$fit.0 = NULL
  old_compare_fit$fit.N = NULL
  
  saveRDS(old_compare_fit, paste0(save_path, fname))
})










