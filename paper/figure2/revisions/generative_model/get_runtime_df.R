library(dplyr)

runtime_path = "~/Dropbox/dropbox_shared/2022. Basilica/simulations/runtimes/generative_model/"
save_path = "~/Dropbox/dropbox_shared/2022. Basilica/simulations/stats_dataframes/revisions/"


# times_sigpr = read.csv(file.path(runtime_path, "last/sigprofiler_exectimes.csv")) %>% 
#   mutate(tool="SigProfiler") %>% tibble::as_tibble() %>% 
#   mutate(execution_time=stringr::str_replace_all(execution_time, " 0:", "00:")) %>% 
#   mutate(execution_time=stringr::str_remove_all(execution_time, " ")) %>% 
#   mutate(execution_time=lubridate::period_to_seconds(lubridate::hms(execution_time))) %>%
#   mutate(execution_time=execution_time/60) %>% 
#   rename(simulation_name=simulation) %>% 
#   select(simulation_name, execution_time, tool)

# this file is with N runs NMF=100
times_sigpr = read.csv(file.path(runtime_path, "old/sigprofiler_exectimes.csv")) %>% 
  mutate(tool="SigProfiler") %>% tibble::as_tibble()

times_sparsesig = read.csv(file.path(runtime_path, "last/sparsesignatures_exectimes.csv")) %>% 
  rename(execution_time=total_mins, simulation_name=name) %>% 
  mutate(tool="SparseSignatures",
                simulation_name=stringr::str_remove_all(simulation_name, ".Rds")) %>% 
  tibble::as_tibble() %>% 
  select(simulation_name, execution_time, tool)

times_sigtoolslib = read.csv(file.path(runtime_path, "signaturetoolslib_exectimes.csv")) %>% 
  mutate(tool="SignatureToolsLib",
                simulation_name=stringr::str_remove_all(simulation_name, ".Rds")) %>% 
  tibble::as_tibble() %>% 
  select(simulation_name, execution_time, tool)

times_sigtoolslib_E = read.csv(file.path(runtime_path, "signaturetoolslib_E_exectimes.csv")) %>% 
  mutate(tool="SignatureToolsLib_E",
                simulation_name=stringr::str_remove_all(simulation_name, ".Rds")) %>% 
  tibble::as_tibble() %>% 
  select(simulation_name, execution_time, tool)

times_bascule = read.csv(file.path(runtime_path, "bascule_exectimes.csv")) %>% 
  rename(execution_time_bascule=execution_time_SBS) %>% tibble::as_tibble() %>% 
  mutate(tool="BASCULE", execution_time=execution_time_bascule) %>%
  select(simulation_name, execution_time, execution_time_bascule, tool)


runtime_df = bind_rows(times_sigpr, 
                       times_sparsesig, 
                       times_sigtoolslib,
                       times_sigtoolslib_E,
                       times_bascule %>% select(-execution_time_bascule)) %>% 
  inner_join(times_bascule %>% select(-execution_time, -tool))

saveRDS(runtime_df, file.path(save_path, "runtime_generative_model.Rds"))

