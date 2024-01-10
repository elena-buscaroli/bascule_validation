library(magrittr)
devtools::load_all("~/GitHub/simbasilica/")
load_basilica()

main_path = "~/GitHub/basilica_validation/comparison_sigprofiler_sparsesig/"
fits_path = "~/Dropbox/shared/2022. Basilica/simulations/fits/fits_dn.matched.2011/"

fitsnames = list.files(fits_path, pattern="simul_fit")

x.simul = readRDS("~/Dropbox/shared/2022. Basilica/simulations/fits/fits_dn.clustering.matched.2011/simul_fit.N150.G3.s5.matched.2011.Rds")$dataset
x = readRDS("~/Dropbox/shared/2022. Basilica/simulations/fits/fits_dn.clustering.matched.2011/simul_fit.N150.G3.s5.matched.2011.Rds")$fit.0

counts = get_input(x)[["SBS"]] %>% dplyr::select(-clusters)
## sigprofiler ####
sp_res = paste0(main_path, "simul_fit.N150.G3.s5.matched.2011/SBS96/Suggested_Solution/SBS96_De-Novo_Solution/")
sp_alpha = read.csv(paste0(sp_res, "Activities/SBS96_De-Novo_Activities_refit.txt"), sep="\t", row.names=1)
sp_sigs = read.csv(paste0(sp_res, "Signatures/SBS96_De-Novo_Signatures.txt"), sep="\t", row.names=1) %>% t()
x.sp = create_basilica_obj(counts, sp_alpha, sp_sigs)



## sparsesignatures ####
ss_res = readRDS(paste0(main_path, "simul_fit.N150.G3.s5.matched.2011.Rds"))
rownames(ss_res$alpha) = rownames(long_to_wide(counts, what="counts"))
ss_res$alpha = ss_res$alpha / rowSums(ss_res$alpha)
x.ss = create_basilica_obj(counts, ss_res$alpha, ss_res$beta)


create_basilica_obj = function(counts, expos, sigs) {
  obj = list(); class(obj) = "basilica_obj"
  obj[["input"]][["SBS"]] = list("counts"=counts,
                                 "reference"=NULL)
  obj[["nmf"]][["SBS"]] = list("exposure"=wide_to_long(expos, what="exposures"),
                               "beta_fixed"=NULL,
                               "beta_denovo"=wide_to_long(sigs, what="beta"))
  
  return(obj)
}






