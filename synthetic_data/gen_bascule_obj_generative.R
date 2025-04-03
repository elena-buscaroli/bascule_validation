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

main_path = "~/Dropbox/dropbox_shared/2022. Basilica/simulations/fits_generative_model/"
fits_path = file.path(main_path, "all_fits/fits_dn.clustering.matched.2011/")
sp_path = file.path(main_path, "sigprofiler_BestSolution/NMF_100/")
ss_path = file.path(main_path, "sparsesignatures/last/")
ss_signature_selection = file.path(ss_path, "signature_selection.csv") %>% 
  read.csv(header=FALSE, col.names=c("filename", "lambd_a", "lambd_b", "K")) %>% 
  mutate(file_id=paste(stringr::str_remove_all(K, "_Signatures"), lambd_a, lambd_b, sep="_")) %>% 
  as_tibble()
stl_path = file.path(main_path, "signaturetoolslib_E/")

save_path = file.path(main_path, "all_fits/fits_dn.matched.2011.compare_LAST/")
dir.create(save_path)

stopifnot(dir.exists(c(main_path, fits_path, sp_path, ss_path, stl_path)))

cli::cli_text("BASCULE fits: {fits_path}\n
              SigProfiler fits: {sp_path}\n
              SparseSignatures fits: {ss_path}\n
              SignatureToolsLib fits: {stl_path}")

fitsnames = list.files(fits_path, pattern="simul_fit")

lapply(fitsnames, function(fname) {
  if (file.exists(paste0(save_path, fname))) return()
  
  print(fname)

  simu_fit = readRDS(paste0(fits_path, fname))

  counts = get_input(simu_fit$dataset)[["SBS"]] %>% dplyr::select(-clusters)
  x.sp = x.ss = x.stl = NULL
  
  ## signature tools lib ####
  try({
    stl_res = readRDS(paste0(stl_path, fname))
    stl_alpha = stl_res$exposures
    stl_alpha = stl_alpha / rowSums(stl_alpha)
    stl_sigs = stl_res$signatures %>% t() %>% as.data.frame()
    x.stl = create_bascule_obj(counts, stl_alpha, stl_sigs)
  })
  
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
    tmp = stringr::str_replace_all(fname, ".Rds", "")
    ss_fname = ss_signature_selection %>% filter(filename==fname) %>% pull(file_id)
    ss_res = readRDS(paste0(ss_path, "signatures/", tmp, "/", ss_fname, ".Rds"))
    rownames(ss_res$alpha) = rownames(long_to_wide(counts, what="counts"))
    ss_res$alpha = ss_res$alpha / rowSums(ss_res$alpha)
    x.ss = create_bascule_obj(counts, ss_res$alpha, ss_res$beta)
  })
  
  simu_fit[["sigprofiler"]] = x.sp
  simu_fit[["sparsesignatures"]] = x.ss
  simu_fit[["signaturetoolslib"]] = x.stl
  
  saveRDS(simu_fit, paste0(save_path, fname))
})









