args = commandArgs(trailingOnly=TRUE)
dataset_id = args[1]  # either "generative_model" or "SigFitTest"
i = as.integer(args[2])

cli::cli_text("Value of i = {i}\n\n")

devtools::load_all("~/GitHub/simbascule/")
devtools::load_all("~/GitHub/bascule/")

run_id = "matched.2011.compare_LAST"

# out_id = paste0(run_id, ".", dataset_id %>% stringr::str_remove_all("/all_fits/"))

main_path_local = "~/Dropbox/dropbox_shared/2022. Basilica/simulations/"
main_path_orfeo = "/orfeo/cephfs/scratch/cdslab/ebusca00/signatures/"

main_path = dplyr::case_when(
  dir.exists(main_path_local) ~ main_path_local,
  dir.exists(main_path_orfeo) ~ main_path_orfeo,
  .default=NULL
)
stopifnot(!is.null(main_path))
rm(main_path_local); rm(main_path_orfeo)

save_path = file.path(main_path, "stats_dataframes/generative_model/")

source("~/GitHub/bascule_validation/synthetic_data/aux_fns/eval_aux_fns.R")
source("~/GitHub/bascule_validation/synthetic_data/aux_fns/plots_aux_fns.R")

# Generate stats dataframe ##### 
if ( grepl("generative_model", dataset_id) ) {
  runids = c("BASCULE", "SigProfiler", "SparseSignatures", "SignatureToolsLib_E", "SignatureToolsLib")
  fitnames = c("x.fit0.auto", "sigprofiler", "sparsesignatures", "signaturetoolslib_E", "signaturetoolslib")
}

if ( grepl("SigFitTest", dataset_id) ) {
  runids = c("BASCULE", "SigProfiler", "SparseSignatures", "SignatureToolsLib_E", "SignatureToolsLib")
  fitnames = c("fit_refined.0", "sigprofiler", "sparsesignatures", "signaturetoolslib_E", "signaturetoolslib")
}

path_local = file.path(main_path, paste0("fits_generative_model/all_fits/fits_dn.", run_id, "/"))
path_orfeo = file.path(main_path, paste0("fits_", dataset_id, "/fits_dn.", run_id, "/"))
path = dplyr::case_when(
  dir.exists(path_local) ~ path_local,
  dir.exists(path_orfeo) ~ path_orfeo,
  .default=NULL
)
stopifnot(!is.null(path))
rm(path_local); rm(path_orfeo)

cli::cli_text("Files path: {path}\n
              Output path: {save_path}\n\n")

if (!dir.exists(save_path)) dir.create(save_path)
stopifnot(all(file.exists(c(path, save_path))))

files = list.files(path, full.names=F, pattern=glob2rx("simul_fit*.Rds"))

fname = file.path(path, files[i])
out_fname = file.path(save_path, files[i] %>% stringr::str_replace_all("simul_fit", "stats_dataframe"))
if (!dir.exists(file.path(path, "clustering"))) dir.create(file.path(path, "clustering"))

if (file.exists(out_fname)) {
  cli::cli_text("File {out_fname} already present. Not saving new file.")
} else {
  stats_dataframe_i = stats_single_data(fname, names_fits=fitnames %>% setNames(runids))
  
  cli::cli_text("Saving object `stats_dataframe_i` in file {out_fname}\n\n")
  saveRDS(stats_dataframe_i, out_fname)
  
}


