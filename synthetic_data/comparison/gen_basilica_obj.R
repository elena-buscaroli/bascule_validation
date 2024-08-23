library(magrittr)
devtools::load_all("~/GitHub/simbascule/")
load_bascule()

create_bascule_obj = function(counts, expos, sigs) {
  obj = list(); class(obj) = "bascule_obj"
  obj[["input"]][["SBS"]] = list("counts"=counts,
                                 "reference"=NULL)
  obj[["nmf"]][["SBS"]] = list("exposure"=wide_to_long(expos, what="exposures"),
                               "beta_fixed"=NULL,
                               "beta_denovo"=wide_to_long(sigs, what="beta"))
  
  return(obj)
}

fits_path = "~/Dropbox/shared/2022. Basilica/simulations/fits/fits_dn.matched.2011/"
sp_path = "~/Dropbox/shared/2022. Basilica/simulations/fits/sigprofiler/NMF_100/"
ss_path = "~/Dropbox/shared/2022. Basilica/simulations/fits/sparsesignatures/"
save_path = "~/Dropbox/shared/2022. Basilica/simulations/fits/fits_dn.matched.2011.compare/"
dir.create(save_path)

fitsnames = list.files(fits_path, pattern="simul_fit")

lapply(fitsnames, function(fname) {
  if (file.exists(paste0(save_path, fname))) return()

  simu_fit = readRDS(paste0(fits_path, fname))
  
  counts = get_input(simu_fit$dataset)[["SBS"]] %>% dplyr::select(-clusters)
  x.sp = x.ss = NULL
  ## sigprofiler ####
  try({
    tmp = stringr::str_replace_all(fname, ".Rds", "")
    sp_res = paste0(sp_path, tmp, "/SBS96/Suggested_Solution/SBS96_De-Novo_Solution/")
    sp_alpha = read.csv(paste0(sp_res, "Activities/SBS96_De-Novo_Activities_refit.txt"), sep="\t", row.names=1)
    sp_alpha = sp_alpha / rowSums(sp_alpha)
    sp_sigs = read.csv(paste0(sp_res, "Signatures/SBS96_De-Novo_Signatures.txt"), sep="\t", row.names=1) %>% t()
    x.sp = create_bascule_obj(counts, sp_alpha, sp_sigs)
  })
  
  ## sparsesignatures ####
  try({
    ss_res = readRDS(paste0(ss_path, fname))
    rownames(ss_res$alpha) = rownames(long_to_wide(counts, what="counts"))
    ss_res$alpha = ss_res$alpha / rowSums(ss_res$alpha)
    x.ss = create_bascule_obj(counts, ss_res$alpha, ss_res$beta)
  })
  
  simu_fit[["sigprofiler"]] = x.sp
  simu_fit[["sparsesignatures"]] = x.ss
  
  saveRDS(simu_fit, paste0(save_path, fname))
})









