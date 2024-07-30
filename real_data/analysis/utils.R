

# load basilica fit
#-----------------------------
# input:
#   - [basilica object, context]
# output:
#   - fixed    (wide format)
#   - denovo   (wide format)
#   - exposure (wide format)
load.basilica <- function(x, type) {
  fixed <- basilica:::get_fixed_signatures(x, matrix = TRUE)[[type]]
  denovo <- basilica:::get_denovo_signatures(x, types = type, matrix = TRUE)[[type]]
  exposure <- basilica:::get_exposure(x, matrix = TRUE)[[type]]
  return(list(fixed=fixed, denovo=denovo, exposure=exposure))
}


#===============================================================================


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


# input
#   x -----> basilica object
#   type --> ("SBS", "DBS")
# output
#   list of:
#   - fixed
#   - denovo
#   - BasilicaExposure
#   - common
#   - rare
#   - SerenaExposure
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
#===============================================================================
#===============================================================================


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
#===============================================================================
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


'
source("/Users/azadsadr/Nextcloud/basilica/scripts/utils.R") # get_data

saveRDS(
  list(
    breast_sbs = get_data(system_type = "mac", organ_type = "Breast", context_type = "SBS"), 
    breast_dbs = get_data(system_type = "mac", organ_type = "Breast", context_type = "DBS"), 
    lung_sbs = get_data(system_type = "mac", organ_type = "Lung", context_type = "SBS"), 
    lung_dbs = get_data(system_type = "mac", organ_type = "Lung", context_type = "DBS"), 
    colorectal_sbs = get_data(system_type = "mac", organ_type = "Colorectal", context_type = "SBS"), 
    colorectal_dbs = get_data(system_type = "mac", organ_type = "Colorectal", context_type = "DBS")
  ), file = "/Users/azadsadr/Nextcloud/basilica/viz_data.Rds"
)
'


#===============================================================================
#===============================================================================
#===============================================================================

# input:
#   - obj       : load.all function output
#   - map       : map.data2 function output
#   - threshold : consider exposure values below this in basilica
exposure_comp_zero <- function(obj, map, threshold) {
  
  total <- obj$SerenaExposure %>% nrow
  
  a <- map$basilica_serena
  
  svec <- c()
  bvec <- c()
  nvec <- c()
  
  for (i in 1:nrow(a)) {
    b <- a[i,] %>% dplyr::pull(Basilica)
    s <- a[i,] %>% dplyr::pull(Serena)
    n <- paste0(a[i,] %>% dplyr::pull(1), "-", a[i,] %>% dplyr::pull(2))
    
    samples <- obj$SerenaExposure %>% basilica:::wide_to_long(what = "exposures") %>% subset(sigs == s & value == 0) %>% dplyr::pull(samples)
    serena_zero <- samples %>% length
    svec[length(svec) + 1] <- serena_zero
    
    vec <- obj$BasilicaExposure[samples, b]
    basilica_low <- vec[vec < threshold] %>% length
    bvec[length(bvec) + 1] <- basilica_low
    
    nvec[length(nvec) + 1] <- n
    
    #print(
    #  paste0(
    #    "total samples: ", total, 
    #    " | serena-zero: ", serena_zero, 
    #    " | basilica-low: ", basilica_low, 
    #    " | percentage: ", round((basilica_low / serena_zero), 2)
    #  )
    #)
  }
  
  df <- tibble::tibble(name = nvec, serena = svec, basilica = bvec, total = rep(total, nrow(a)))
  
  
  return(df)
}


#===============================================================================
#===============================================================================
#===============================================================================

# input  : 
#  - denovo signatures (long format)
#  - data object
# output : 
# - denovo signatures with mapped labels from reference catalogue (long format)
map_denovo2reference <- function(x, data, type) {
  
  dn <- get_denovo_signatures(x)[[type]]
  
  map_obj <- data[[paste(tolower(type), "_map", sep = "")]][["denovo_reference"]]
  
  if ( (map_obj %>% nrow) > 0 ) {
    denovo <- dn %>% left_join(map_obj, by = c("sigs" = "first"))
    
    denovo <- denovo %>%
      mutate(sigs = ifelse(!is.na(second), paste0("R-", second), sigs)) %>%
      select(-c(second, values))
  }
  
  return(denovo)
}

#-------------------------------------------------------------------------------

name_mapping <- function(x, data, type) {
  dn <- map_denovo2reference(x, data, type)
  if (type == "SBS") {
    x$nmf$SBS$beta_denovo <- dn
  } else if (type == "DBS") {
    x$nmf$DBS$beta_denovo <- dn
  } else {
    warning("Context is unknown!")
  }
  return(x)
}

#-------------------------------------------------------------------------------


map_basilica2serena <- function(signatures, map) {
  
  map_obj <- map$basilica_serena %>% dplyr::rename("value2" = "value")
  
  signatures <- signatures %>%
    left_join(map_obj, by = c("sigs" = "Basilica"))
  
  signatures <- signatures %>%
    mutate(sigs = ifelse(!is.na(Serena), paste0("S-", Serena), sigs)) %>%
    select(-c(Serena, value2))
  
  return(signatures)
}


#===============================================================================
#=================================== STORAGE ===================================
#===============================================================================

'
# input: data.frame of cosine similarity matrix
# output: return list of:
#   - maximum value of cosine similarity
#   - respective signatures
#   - cosine similarity matrix where its max value is set to zero
maxfinder <- function(x) {
  
  # check if data.frame has both dimnesion > 0
  if ( (nrow(x)==0) | (ncol(x)==0) ) {
    return(list(matrix=NULL, first=NULL, second=NULL, value=NULL, end=FALSE))
  }
  
  # check if all values of data.frame are zero
  if (all(x == 0)) {
    return(list(matrix=NULL, first=NULL, second=NULL, value=NULL, end=TRUE))
  }
  
  if ( (nrow(x)==1) & (ncol(x)==1) ) {
    a <- rownames(x)
    b <- colnames(x)
    c <- x[1,1]
    x[1,1] <- 0
    return(list(matrix=x, first=a, second=b, value=c, end=FALSE))
  }
  else {
    row <- which(x == max(x), arr.ind = T)[1]
    col <- which(x == max(x), arr.ind = T)[2]
    a <- rownames(x[row, ])
    b <- colnames(x[col])
    c <- max(x)
    x[row, ] <- 0
    x[, col] <- 0
    return(list(matrix=x, first=a, second=b, value=c, end=FALSE))
  }
}


#-------------------------------------------------------------------------------

mapper <- function(x) {
  ind <- max(2, min(nrow(x), ncol(x)))
  
  first <- c()
  second <- c()
  values <- c()
  
  for (i in 1:(ind-1)) {
    
    res <- maxfinder(x)
    
    if ( is.null(res[[1]]) ) {
      print("invalid input!")
      return(0)
    }
    else if ( res[[5]] ) {
      break
    }
    else {
      x <- res[[1]]
      first[i] <- res[[2]]
      second[i] <- res[[3]]
      values[i] <- res[[4]]
    }
  }
  
  df <- data.frame(first, second, values)
  return(df)
}

#-------------------------------------------------------------------------------

# basilica_exposure # basilica brca exposure
# serena_exposure   # serena brca exposure


# reference --> reference signature (e.g., COSMIC)
# fixed ------> basilica fixed signatures (wide)
# denovo -----> basilica denovo signatures (wide)
# common -----> serena common signatures (wide)
# rare -------> serena rare signatures (wide)
# threshold --> real number 
## (if cosine similarity between a basilica signature and serena signature is
## higher than threshold, we map them to eachother)
map.data <- function(
    reference, 
    fixed, 
    denovo, 
    common, 
    rare, 
    threshold
) {
  
  cmatrix1 <- basilica:::cosine.matrix(denovo, common)
  denovo_common <- mapper(cmatrix1)
  denovo_common <- subset(denovo_common, values >= threshold )
  
  cmatrix2 <- basilica:::cosine.matrix(denovo, rare)
  denovo_rare <- mapper(cmatrix2)
  denovo_rare <- subset(denovo_rare, values >= threshold )
  
  cmatrix3 <- basilica:::cosine.matrix(denovo, reference)
  denovo_reference <- mapper(cmatrix3)
  denovo_reference <- subset(denovo_reference, values >= threshold )
  
  cmatrix4 <- basilica:::cosine.matrix(fixed, common)
  fixed_common <- mapper(cmatrix4)
  fixed_common <- subset(fixed_common, values >= threshold )
  
  cmatrix5 <- basilica:::cosine.matrix(fixed, rare)
  fixed_rare <- mapper(cmatrix5)
  fixed_rare <- subset(fixed_rare, values >= threshold )
  
  cmatrix6 <- basilica:::cosine.matrix(fixed, reference)
  fixed_reference <- mapper(cmatrix6)
  fixed_reference <- subset(fixed_reference, values >= threshold )
  
  cmatrix7 <- basilica:::cosine.matrix(common, reference)
  common_reference <- mapper(cmatrix7)
  common_reference <- subset(common_reference, values >= threshold )
  
  cmatrix8 <- basilica:::cosine.matrix(rare, reference)
  rare_reference <- mapper(cmatrix8)
  rare_reference <- subset(rare_reference, values >= threshold )
  
  cmatrix9 <- basilica:::cosine.matrix(rbind(fixed, denovo), rbind(common, rare))
  basilica_serena <- mapper(cmatrix9)
  basilica_serena <- subset(basilica_serena, values >= threshold )
  colnames(basilica_serena) <- c("Basilica", "Serena", "value")
  
  # aggregated fixed and denovo signatures (basilica)
  all_basilica <- rbind(fixed, denovo)
  # list of signatures name detected by basilica but un-explained by serena
  ne_basilica <- rownames(all_basilica[!(rownames(all_basilica) %in% basilica_serena$Basilica), ])
  
  # aggregated common and rare signatures (serena)
  all_serena <- rbind(common, rare)
  # list of signatures name detected by serena but un-explained by basilica
  ne_serena <- rownames(all_serena[!(rownames(all_serena) %in% basilica_serena$Serena), ])
  
  # -------------------------------- [ OUTPUT ] --------------------------------
  
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

#-------------------------------------------------------------------------------

'


