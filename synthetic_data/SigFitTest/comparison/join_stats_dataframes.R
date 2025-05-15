library(tidyverse)

# main_path = "~/Dropbox/dropbox_shared/2022. Basilica/simulations/"
main_path = "/orfeo/cephfs/scratch/cdslab/ebusca00/signatures/"

save_path = file.path(main_path, "stats_dataframes/")
files_path = file.path(main_path, "stats_dataframes/SigFitTest/")

files = list.files(files_path, full.names=T)

dataset_id = "SigFitTest"
run_id = "matched.2011.compare_LAST"
out_id = paste0(run_id, ".", dataset_id %>% stringr::str_remove_all("/all_fits/"))
out_file = paste0(save_path, "stats_", out_id, ".Rds")

final_stats_dataframe = tibble::tibble()
for (i in 1:length(files)) {
  fname = files[i]
  cli::cli_text("File: {fname}")
  stats_dataframe_i = readRDS(fname)
  final_stats_dataframe = rbind(
    final_stats_dataframe,
    stats_dataframe_i
  )
  gc()
}

saveRDS(final_stats_dataframe, out_file)

# for (i in 1:length(files)) {
#   fname = files[i]
#   stats_dataframe_i = readRDS(fname)
#   write.table(stats_dataframe_i, out_file, sep=",", col.names=!file.exists(out_file), append=TRUE)
#   # write.table(myDF, "myDF.csv", sep = ",", col.names = !file.exists("myDF.csv"), append = T)
# }

