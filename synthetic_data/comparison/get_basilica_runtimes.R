library(dplyr)

fits_path = "~/Dropbox/dropbox_shared/2022. Basilica/simulations/fits/fits_dn.matched.2011/"
files = list.files(fits_path, full.names=T, pattern=".Rds")

runtimes = lapply(files, function(fname) {
  cat(fname); cat("\n")
  simul_fit = readRDS(fname)
  t_sbs = simul_fit$fit.0.auto$nmf$SBS$pyro$time %>% stringr::str_remove_all("Time difference of | mins") %>% as.numeric()
  t_dbs = simul_fit$fit.0.auto$nmf$DBS$pyro$time %>% stringr::str_remove_all("Time difference of | mins") %>% as.numeric()
  data.frame(simulation_name=strsplit(fname,"//")[[1]][2] %>% stringr::str_remove_all(".Rds"),
             execution_time_SBS=t_sbs,
             execution_time_DBS=t_dbs)
}) %>% bind_rows()

write.csv(runtimes, "~/Dropbox/dropbox_shared/2022. Basilica/simulations/runtimes/bascule_exectimes.csv", row.names=F)
