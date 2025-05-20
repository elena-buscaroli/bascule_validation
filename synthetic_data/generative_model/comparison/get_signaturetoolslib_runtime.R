library(dplyr)

tool_name = "signaturetoolslib"

fits_path = file.path("~/Dropbox/dropbox_shared/2022. Basilica/simulations/fits_generative_model/",
                      tool_name)
files = list.files(fits_path, full.names=T, pattern=glob2rx("simul_fit*.Rds"))

runtimes = lapply(files, function(fname) {
  cat(fname); cat("\n")
  simul_fit = readRDS(fname)
  t1 = simul_fit$time %>% stringr::str_remove_all("Time difference of | mins") %>% as.numeric()
  # t2 = simul_fit$extraction_time %>% stringr::str_remove_all("Time difference of | mins") %>% as.numeric()
  data.frame(simulation_name=strsplit(fname,"//")[[1]][2] %>% stringr::str_remove_all(".Rds") %>% 
               stringr::str_remove_all("signaturetoolslib/"),
             execution_time=t1
             # execution_time_extraction=t2
             )
}) %>% bind_rows()

write.csv(runtimes, paste0("~/Dropbox/dropbox_shared/2022. Basilica/simulations/runtimes/generative_model/",
                           tool_name, "_exectimes.csv"), row.names=F)
