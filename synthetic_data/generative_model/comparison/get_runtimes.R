library(dplyr)

main_path = "~/Dropbox/dropbox_shared/2022. Basilica/simulations/fits_generative_model/"
save_path = "~/Dropbox/dropbox_shared/2022. Basilica/simulations/"

# BASCULE #####

fits_path = file.path(main_path, "all_fits/fits_dn.matched.2011.compare_LAST")
files = list.files(fits_path, full.names=T, pattern=".Rds")

runtimes = lapply(files, function(fname) {
  cat(fname); cat("\n")
  simul_fit = readRDS(fname)
  t_sbs = simul_fit$x.fit0.auto$nmf$SBS$pyro$time %>% as.numeric(units="mins")
  t_dbs = simul_fit$x.fit0.auto$nmf$DBS$pyro$time %>% as.numeric(units="mins")
  data.frame(simulation_name=strsplit(fname,"//")[[1]][2] %>% stringr::str_remove_all(".Rds") %>% 
               stringr::str_remove_all("all_fits/fits_dn.matched.2011.compare_LAST/"),
             execution_time_SBS=t_sbs,
             execution_time_DBS=t_dbs)
}) %>% bind_rows()

write.csv(runtimes, file.path(save_path, "runtimes/generative_model/bascule_exectimes.csv"), row.names=F)


# SignatureToolsLib #####

fits_path = file.path(main_path, "signaturetoolslib")
files = list.files(fits_path, full.names=T, pattern=glob2rx("simul_fit*.Rds"))

runtimes = lapply(files, function(fname) {
  cat(fname); cat("\n")
  simul_fit = readRDS(fname)
  t1 = simul_fit$time %>% as.numeric(units="mins")
  # t2 = simul_fit$extraction_time %>% stringr::str_remove_all("Time difference of | mins") %>% as.numeric()
  data.frame(simulation_name=strsplit(fname,"//")[[1]][2] %>% stringr::str_remove_all(".Rds") %>% 
               stringr::str_remove_all("signaturetoolslib/"),
             execution_time=t1
             # execution_time_extraction=t2
  )
}) %>% bind_rows()

write.csv(runtimes, file.path(save_path, "runtimes/generative_model/signaturetoolslib_exectimes.csv"), row.names=F)


# SignatureToolsLib_E #####

fits_path = file.path(main_path, "signaturetoolslib_E")
files = list.files(fits_path, full.names=T, pattern=glob2rx("simul_fit*.Rds"))

runtimes = lapply(files, function(fname) {
  cat(fname); cat("\n")
  simul_fit = readRDS(fname)
  t1 = simul_fit$time %>% as.numeric(units="mins")
  t2 = simul_fit$extraction_time %>% as.numeric(units="mins")
  data.frame(simulation_name=strsplit(fname,"//")[[1]][2] %>% stringr::str_remove_all(".Rds") %>% 
               stringr::str_remove_all("signaturetoolslib_E/"),
             execution_time=t1,
             execution_time_extraction=t2
  )
}) %>% bind_rows()

write.csv(runtimes, file.path(save_path, "runtimes/generative_model/signaturetoolslib_E_exectimes.csv"), row.names=F)

