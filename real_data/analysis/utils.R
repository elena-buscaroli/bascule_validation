

# load basilica fit
#-----------------------------
# input:
#   - x        : basilica object
#   - type     : context e.g., ("SBS", "DBS")
# output:
#   - fixed    : fixed signatures  (wide format)
#   - denovo   : denovo signatures (wide format)
#   - exposure : exposure matrix   (wide format)
load.basilica <- function(x, type) {
  fixed <- basilica:::get_fixed_signatures(x, matrix = TRUE)[[type]]
  denovo <- basilica:::get_denovo_signatures(x, types = type, matrix = TRUE)[[type]]
  exposure <- basilica:::get_exposure(x, matrix = TRUE)[[type]]
  return(list(fixed=fixed, denovo=denovo, exposure=exposure))
}


#===============================================================================


# load Degasperi fit (based on organ type)
#-----------------------------
# input:
#   - x        : basilica object   -  
#   - reference_path : path to the reference catalogue from Degasperi study (character)
#   - exposure_path  : path to the exposure matrix from Degasperi study (character)
#   - common_names   : list of common signatures present in corresponding organ type in Degasperi study (character vector)
#   - rare_names     : list of rare signatures present in corresponding organ type in Degasperi study (character vector)
#   - organ_type     : tumor type (character)
# output:
#   - common   : common signatures (dataframe, wide format)
#   - rare     : rare signatures   (dataframe, wide format)
#   - exposure : exposure matrix   (dataframe, wide format)

load.serena <- function(reference_path, exposure_path, common_names, rare_names, organ_type) {
  
  # ------------------------------ REFERENCE -----------------------------------
  # load serena reference data.frame
  reference.T <- read.table(file = reference_path, sep = '\t', header = TRUE)
  reference <- as.data.frame(t(reference.T)) # data.frame (wide)
  
  # QUALITY CONTROL - if all common and rare signatures are present in reference catalogue
  if ( !( all(common_names %in% (reference %>% rownames)) & all(rare_names %in% (reference %>% rownames)) ) ) {
    warning("not all common or rare signatures are present in reference")
  }
  
  # ------------------------------- EXPOSURE -----------------------------------
  # load exposure data.frame (all - unnormalized)
  exposure0 <- read.table(file = exposure_path, sep = '\t', header = TRUE)
  
  # filter exposure based on organ type (unnormalized)
  exposure0 <- subset(exposure0, organ==organ_type)
  
  # removing the cohort and organ columns  (unnormalized)
  exposure0 <- exposure0[!(names(exposure0) %in% c("cohort", "organ"))]
  
  # removing columns (signatures) where all are zeros (unnormalized)
  exposure0 <- exposure0 %>% dplyr::select(where(~ any(. != 0))) # dplyr required
  
  # delete rows with NaN values
  #exposure1 <- na.omit(exposure)
  
  # normalizing
  exposure <- sweep(exposure0, 1, rowSums(exposure0), "/")
  
  
  # QUALITY CONTROL - if all common and rare signatures are present in exposures matrix
  if ( !( all(common_names %in% (exposure %>% colnames)) & all(rare_names %in% (exposure %>% colnames)) ) ) {
    warning("not all common or rare signatures are present in exposure")
  }
  
  
  # selecting exposures involved in specified organ
  exposure <- exposure[, c(common_names, rare_names)]
  
  # delete rows with NaN values
  exposure <- na.omit(exposure)
  
  common <- reference[common_names, ]
  rare <- reference[rare_names, ]
  
  return(list(common=common, rare=rare, exposure=exposure))
}


#===============================================================================


# quality control of basilica and Degasperi fit data
#-----------------------------
# input:
#   - fixed            : fixed signatures (dataframe, wide format)
#   - denovo           : denovo signatures (dataframe, wide format)
#   - BasilicaExposure : exposure matrixx from basilica (dataframe, wide format)
#   - common           : common signatures (dataframe, wide format)
#   - rare             : rare signatures (dataframe, wide format)
#   - SerenaExposure   : exposure matrix from Degasperi (dataframe, wide format)
#   - reference        : reference catalogue (dataframe, wide format)
# output:
#   - fixed            : fixed signatures (dataframe, wide format)
#   - denovo           : denovo signatures (dataframe, wide format)
#   - BasilicaExposure : exposure matrixx from basilica (dataframe, wide format)
#   - common           : common signatures (dataframe, wide format)
#   - rare             : rare signatures (dataframe, wide format)
#   - SerenaExposure   : exposure matrix from Degasperi (dataframe, wide format)
#   - reference        : reference catalogue (dataframe, wide format)

input.qc <- function(fixed, denovo, BasilicaExposure, common, rare, SerenaExposure, reference) {
#run.QC <- function(basilica, serena, reference) {
  
  a <- all((colnames(fixed) %in% colnames(reference))==TRUE)
  b <- all((colnames(denovo) %in% colnames(reference))==TRUE)
  c <- all((colnames(common) %in% colnames(reference))==TRUE)
  d <- all((colnames(rare) %in% colnames(reference))==TRUE)
  
  if (!(a & b & c & d)) {
    warning("there are mismatches between contexts!")
    #return(NULL)
  }
  
  fixed <- fixed[, colnames(reference)]
  denovo <- denovo[, colnames(reference)]
  common <- common[, colnames(reference)]
  rare <- rare[, colnames(reference)]
  
  # extracting common samples between basilica and serena exposure matrices
  intersected_samples <- intersect(rownames(BasilicaExposure), rownames(SerenaExposure))
  BasilicaExposure <- BasilicaExposure[intersected_samples, ]
  SerenaExposure <- SerenaExposure[intersected_samples, ]
  
  #return(list(basilica=basilica, serena=serena))
  return(
    list(
    fixed=fixed, denovo=denovo, BasilicaExposure=BasilicaExposure, 
    common=common, rare=rare, SerenaExposure=SerenaExposure
    )
  )
}


#===============================================================================


# loading all fit data from both bascule and Degasperi
#-----------------------------
# input
#   - basilica_obj    : basilica object
#   - serena_ref_path : path to the reference catalogue from Degasperi study (character)
#   - serena_exp_path : path to the exposure matrix from Degasperi study (character)
#   - common_names    : list of common signatures present in corresponding organ type in Degasperi study (character vector)
#   - rare_names      : list of rare signatures present in corresponding organ type in Degasperi study (character vector)
#   - organ           : tumor type (character)
#   - reference       : reference catalogue (dataframe, wide format)
#   - type            : context e.g., ("SBS", "DBS")
# output
#   - fixed            : fixed signatures from bascule (dataframe, wide format)
#   - denovo           : denovo signatures from bascule (dataframe, wide format)
#   - BasilicaExposure : exposure matrix from bascule (dataframe, wide format)
#   - common           : common signatures from Degasperi (dataframe, wide format)
#   - rare             : rare signatures from Degasperi (dataframe, wide format)
#   - SerenaExposure   : exposure matrix from Degasperi (dataframe, wide format)

load.all <- function(basilica_obj, serena_ref_path, serena_exp_path, common_names, rare_names, organ, reference, type) {
  
  basilica <- load.basilica(
    basilica_obj, 
    type
  )
  
  serena <- load.serena(
    reference_path=serena_ref_path, 
    exposure_path=serena_exp_path, 
    common_names, 
    rare_names, 
    organ_type=organ
  )
  
  #qc <- run.QC(basilica=basilica, serena=serena, reference=basilica::COSMIC_sbs)
  qc <- input.qc(
    fixed = basilica$fixed, 
    denovo = basilica$denovo, 
    BasilicaExposure = basilica$exposure, 
    common = serena$common, 
    rare = serena$rare, 
    SerenaExposure = serena$exposure, 
    reference = reference
    )
  
  return(list(
    fixed=qc$fixed, denovo=qc$denovo, BasilicaExposure=qc$BasilicaExposure, 
    common=qc$common, rare=qc$rare, SerenaExposure=qc$SerenaExposure
    )
  )
}


#===============================================================================


# return the ordered (cosine similarity) paired signature (rows) datafrmae 
#-----------------------------
# input
#   - cosine_matrix   : ?
#   - threshold       : ?
# output
#   - filtered_data_B : ?

mapper2 <- function(cosine_matrix, threshold=0.8) {
  # Convert to long format
  long_df <- cosine_matrix %>%
    tibble::rownames_to_column(var = "first") %>%
    tidyr::pivot_longer(-first, names_to = "second", values_to = "values")
  
  # Order by cosine similarity
  ordered_pairs <- long_df %>%
    dplyr::arrange(desc(values)) %>% subset(values > threshold)
  
  filtered_data_A <- ordered_pairs %>% dplyr::distinct(.data[["first"]], .keep_all = TRUE)
  filtered_data_B <- filtered_data_A %>% dplyr::distinct(.data[["second"]], .keep_all = TRUE)
  
  return(filtered_data_B)
}


#===============================================================================


# return the mapped signatures from Bascule and Degasperi study
#-----------------------------
# input
#   - reference : reference catalogue (dataframe, wide format)
#   - fixed     : fixed signatures from bascule (dataframe, wide format)
#   - denovo    : denovo signatures from bascule (dataframe, wide format)
#   - common    : common signatures from Degasperi study (dataframe, wide format)
#   - rare      : rare signatures from Degasperi study (dataframe, wide format)
#   - threshold : map two signatures with cosine similarity higher than threshold (numeric)
# output
#   - denovo_common    : paired mapped signatures - denovo vs. common (dataframe, wide format)
#   - denovo_rare      : paired mapped signatures - denovo vs. rare (dataframe, wide format)
#   - denovo_reference : paired mapped signatures - denovo vs. reference (dataframe, wide format)
#   - fixed_common     : paired mapped signatures - fixed vs. common (dataframe, wide format)
#   - fixed_rare       : paired mapped signatures - fixed vs. rare (dataframe, wide format)
#   - fixed_reference  : paired mapped signatures - fixed vs. reference (dataframe, wide format)
#   - common_reference : paired mapped signatures - common vs. reference (dataframe, wide format)
#   - rare_reference   : paired mapped signatures - rare vs. reference (dataframe, wide format)
#   - basilica_serena  : paired mapped signatures - fixed + denovo vs. common + rare (dataframe, wide format)
#   - basilica_singles : signatures from Bascule not mapped to any signatures
#   - serena_singles   : signatures from Degasperi not mapped to any signatures

map.data2 <- function(
    reference, 
    fixed, 
    denovo, 
    common, 
    rare, 
    threshold
) {
  
  cmatrix1 <- basilica:::cosine.matrix(denovo, common)
  denovo_common <- mapper2(cmatrix1, threshold)
  
  cmatrix2 <- basilica:::cosine.matrix(denovo, rare)
  denovo_rare <- mapper2(cmatrix2, threshold)
  
  cmatrix3 <- basilica:::cosine.matrix(denovo, reference)
  denovo_reference <- mapper2(cmatrix3, threshold)
  
  cmatrix4 <- basilica:::cosine.matrix(fixed, common)
  fixed_common <- mapper2(cmatrix4, threshold)
  
  cmatrix5 <- basilica:::cosine.matrix(fixed, rare)
  fixed_rare <- mapper2(cmatrix5, threshold)
  
  cmatrix6 <- basilica:::cosine.matrix(fixed, reference)
  fixed_reference <- mapper2(cmatrix6, threshold)
  
  cmatrix7 <- basilica:::cosine.matrix(common, reference)
  common_reference <- mapper2(cmatrix7, threshold)
  
  cmatrix8 <- basilica:::cosine.matrix(rare, reference)
  rare_reference <- mapper2(cmatrix8, threshold)
  
  cmatrix9 <- basilica:::cosine.matrix(rbind(fixed, denovo), rbind(common, rare))
  basilica_serena <- mapper2(cmatrix9, threshold)
  colnames(basilica_serena) <- c("Basilica", "Serena", "value")
  
  # aggregated fixed and denovo signatures (basilica)
  all_basilica <- rbind(fixed, denovo)
  # list of signatures name detected by basilica but un-explained by serena
  ne_basilica <- rownames(all_basilica[!(rownames(all_basilica) %in% basilica_serena$Basilica), ])
  
  # aggregated common and rare signatures (serena)
  all_serena <- rbind(common, rare)
  # list of signatures name detected by serena but un-explained by basilica
  ne_serena <- rownames(all_serena[!(rownames(all_serena) %in% basilica_serena$Serena), ])
  
  data <- list(
    denovo_common = denovo_common, 
    denovo_rare = denovo_rare, 
    denovo_reference = denovo_reference, 
    fixed_common = fixed_common, 
    fixed_rare = fixed_rare, 
    fixed_reference = fixed_reference, 
    common_reference = common_reference, 
    rare_reference = rare_reference, 
    basilica_serena = basilica_serena, 
    basilica_singles = ne_basilica, 
    serena_singles = ne_serena
  )
  return(data)
}


#===============================================================================


get_data <- function(system_type, organ_type, context_type) {
  
  basilica_package_path_r <- "Documents/packages/basilica"
  
  
  basilica_obj_before_r <- paste("Nextcloud/basilica/fit/fit_", organ_type, ".Rds", sep = "")
  basilica_obj_after_r <- paste("Nextcloud/basilica/fit/fit_", organ_type, "_RC.Rds", sep = "")
  
  serena_ref_path_sbs_r <- "Nextcloud/basilica/serena/SBS_v2.03/RefSig_SBS_v2.03.tsv"
  serena_ref_path_dbs_r <- "Nextcloud/basilica/serena/DBS_v1.01/RefSig_DBS_v1.01.tsv"
  
  serena_exp_path_sbs_r <- "Nextcloud/basilica/serena/SBS_v2.03/RefSig_SBS_Exposures_v2.03.tsv"
  serena_exp_path_dbs_r <- "Nextcloud/basilica/serena/DBS_v1.01/RefSig_DBS_Exposures_v1.01.tsv"
  
  tumors_sigs_r <- "Nextcloud/basilica/3T_sigs.rds"
  
  #-------------------------------------------------------------------------------
  
  if ( system_type == "mac" ) {
    path <- "/Users/azadsadr/"
  } else if ( system_type == "linux" ) {
    path <- "/home/azad/"
  } else warning("invalid system type! choose between mac and linux")
  
  if ( !(organ_type %in% c("Breast", "Lung", "Colorectal")) ) warning("invalid context type! choose between Breast, Lung and Colorectal")
  
  if ( !(context_type %in% c("SBS", "DBS")) ) warning("invalid context type! choose between SBS and DBS")
  
  devtools::load_all( paste(path, basilica_package_path_r, sep = "") )
  tumor_sigs <- readRDS( paste(path, tumors_sigs_r, sep = "") )
  tumor_sigs$sbs_colorectal_common[tumor_sigs$sbs_colorectal_common == "SBS17a"] <- "SBS17" # CORRECTION SBS17a removed - it has zero exposure for all samples in colorectal
  
  x_before <- readRDS( paste(path, basilica_obj_before_r, sep = "") )
  x_after <- readRDS( paste(path, basilica_obj_after_r, sep = "") )
  
  obj <- load.all(
    basilica_obj = x_after, 
    serena_ref_path = paste(path, ifelse( (context_type=="SBS") , serena_ref_path_sbs_r, serena_ref_path_dbs_r), sep = ""), 
    serena_exp_path = paste(path, ifelse( (context_type=="SBS") , serena_exp_path_sbs_r, serena_exp_path_dbs_r), sep = ""), 
    common_names = tumor_sigs[ paste( tolower(context_type), tolower(organ_type), "common", sep = "_") ][[1]], 
    rare_names = tumor_sigs[ paste( tolower(context_type), tolower(organ_type), "rare", sep = "_") ][[1]], 
    organ = organ_type, 
    reference = ( if(context_type=="SBS") basilica::COSMIC_sbs_filt else basilica::COSMIC_dbs ), 
    type = context_type
  )
  
  map <- map.data2(
    reference = ( if(context_type=="SBS") basilica::COSMIC_sbs_filt else basilica::COSMIC_dbs ), 
    fixed = obj$fixed, 
    denovo = obj$denovo, 
    common = obj$common, 
    rare = obj$rare, 
    threshold = 0.8
  )
  
  return( list(x_before=x_before, x_after=x_after, obj=obj,  map=map) )
}








