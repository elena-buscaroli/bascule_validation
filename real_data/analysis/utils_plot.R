library(ggplot2)
library(ggrepel)
library(ggpubr)
library(reshape2)


'
remove.sparsity <- function(x, type, exposure_thr=0.01) {
  x[x < exposure_thr] <- 0
  df <- data.frame(type = c(), name = c(), mean = c(), freq = c())
  
  if ( (x %>% ncol) < 1 ) { return(df) }
  
  for (i in 1:ncol(x)) {
    ls <- list(
      type = type, 
      name = colnames(x)[i], 
      mean = sum(x[,i]) / sum(x[,i] != 0), 
      freq = sum(x[,i] != 0)
    )
    df <- rbind(df, ls)
    df$type <- factor(df$type)
  }
  return(df)
}
'

# input:
#  - x: exposure (wide format)
#  - type: basilica | serena
#  - exposure_thr: threshold to eliminate exposure values
# output:
#  - dataframe of signatures, their mean exposure, frequency and type
remove_sparsity <- function(x, type, exposure_thr=0.01) {
  aaa <- basilica:::wide_to_long(x, what = "exposures") %>% subset( value > exposure_thr )
  bbb <- data.frame(
    mean = tapply(aaa$value, aaa$sigs, mean), 
    freq = tapply(aaa$value, aaa$sigs, length), 
    type = type
  ) %>% tibble::rownames_to_column("name")
  return(bbb)
}

# ==============================================================================


library(dplyr)
# input
#   basilica_exposure --> wide dataframe
#   serena_exposure ----> wide dataframe
#   basilica_singles ---> unmapped basilica signatures (character vector)
#   serena_singles -----> unmapped serena signatures (character vector)
# output
#   ggplot object
plot_unmapped <- function(basilica_exp, serena_exp, basilica_singles, serena_singles, exposure_thr=0.01) {
  
  if ( (basilica_singles %>% length) > 0 ) {
    a <- remove_sparsity(
      ( basilica_exp %>% select( all_of(basilica_singles) ) ), 
      "Basilica", 
      exposure_thr
    )
  } else {
    a <- NULL
  }
  
  if ( (serena_singles %>% length) > 0 ) {
    b <- remove_sparsity(
      (serena_exp %>% select( all_of(serena_singles[serena_singles %in% colnames(serena_exp)]) )), 
      "Serena", 
      exposure_thr
    )
  } else {
    b <- NULL
  }
  
  df <- rbind(a, b)
  
  #df <- df %>% replace(is.na(.), 0)
  
  p <- ggplot(
    data = df, 
    aes(x = freq, y = mean, label=name)) + 
    geom_point(aes(color = factor(type)), size=2) + 
    #ggtitle(paste0("Unexplained Signatures (exposure > ", exposure_thr*100, "%)")) + 
    #xlab("Number of Samples") + 
    #ylab("Mean Exposure") + 
    labs(
      title = "Unexplained Signatures Mean Exposures", 
      subtitle = paste0("(exposures > ", exposure_thr*100, "%)"), 
      caption = "samples with exposure below the threshold are removed", 
      x="Number of Samples", 
      y="Mean Exposure", 
      color='Tools'
    ) + 
    geom_label_repel(
      aes(label = name),
      box.padding = 0.35, 
      point.padding = 0.5,
      segment.color = 'grey50', 
      max.overlaps = 25
    ) +
    expand_limits(y=c(0,1)) + 
    theme_bw()
  
  return(p)
}

# ==============================================================================


plot_exposure_comp <- function(obj, map) {
  
  basilica_exp <- obj$BasilicaExposure
  serena_exp <- obj$SerenaExposure
  basilica_serena <- map$basilica_serena
  
  p_list <- list()
  for (i in 1:nrow(basilica_serena)) {
    df <- data.frame(
      b_pair = basilica_exp[, basilica_serena[i,]$Basilica], 
      s_pair = serena_exp[, basilica_serena[i,]$Serena]
    )
    p_list[[length(p_list)+1]] <- ggplot(df, aes(x=b_pair, y=s_pair)) + 
      geom_point() + 
      labs(
        title = "", 
        x = paste0(basilica_serena[i,]$Basilica, " exposure (basilica)"), 
        y = paste0(basilica_serena[i,]$Serena, " exposure (serena)")
      ) + 
      #xlab(paste0(basilica_serena[i,]$Basilica, " exposure (basilica)")) +
      #ylab(paste0(basilica_serena[i,]$Serena, " exposure (serena)")) + 
      ylim(c(0,1)) + 
      theme_bw()
  }
  
  p <- patchwork::wrap_plots(p_list) & patchwork::plot_annotation(
    title = "Exposure Alignment of Mapped Signatures", 
    subtitle = "(Basilica vs. Serena)"
    )
  
  return(p)
}


#===============================================================================

# input : map function output
# output: ggplot object - heatmap of similarity between basilica, serena and refrence signatures
plot_heatmp_comp <- function(map) {
  
  fixed_common <- map$fixed_common %>% dplyr::mutate(rowType="Fixed", colType="Common (serena)")
  fixed_rare <- map$fixed_rare %>% dplyr::mutate(rowType="Fixed", colType="Rare (serena)")
  fixed_reference <- map$fixed_reference %>% dplyr::mutate(rowType="Fixed", colType="COSMIC")
  
  denovo_common <- map$denovo_common %>% dplyr::mutate(rowType="Denovo", colType="Common (serena)")
  denovo_rare <- map$denovo_rare %>% dplyr::mutate(rowType="Denovo", colType="Rare (serena)")
  denovo_reference <- map$denovo_reference %>% dplyr::mutate(rowType="Denovo", colType="COSMIC")
  
  all <- rbind(
    fixed_common, 
    fixed_rare, 
    fixed_reference, 
    denovo_common, 
    denovo_rare, 
    denovo_reference
  )
  
  all$colType <- factor(all$colType, levels = c("Common (serena)", "Rare (serena)", "COSMIC"))
  
  p <- ggplot(data = all, aes(x = first, y = second, fill = values)) + 
    geom_tile(color="black") + 
    facet_grid(colType ~ rowType, scales="free", space = "free") + 
    geom_text(aes(label=round(values, digits = 2))) + 
    scale_fill_gradient(low = "gray", high = "dodgerblue1") + 
    labs(
      title = "Signatures Mapping", 
      x = "Basilica Signatures", 
      y = "[COSMIC + Serena] Signatures", 
      subtitle = "Basilica vs. [COSMIC + Serena]", 
      caption = "COSMIC v3.4 - October 2023", 
      fill = "Cosine Similarity") + 
    theme_bw() + 
    theme(
      legend.title = element_text(angle = 90, vjust = 1), 
      axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)
      )
  
  return(p)
}


#===============================================================================

plot_scores_nmf <- function(x, types=get_types(x), remove_outliers=FALSE) {
  
  a <- get_scores(x) %>% subset(score_id == "bic" & type %in% types & parname == "K") %>%
    dplyr::group_by(type) %>%
    dplyr::mutate(is.min=score==min(score),
                  label=replace(is.min, is.min==T, "Best fit"), 
                  label=replace(label, label==F, NA))
  
  seeds <- a %>% subset( (label == "Best fit") & (type %in% types) ) %>% dplyr::pull(seed) %>% unique
  
  scores <- a %>% subset( (seed %in% seeds) & (type %in% types) )
  
  if (remove_outliers)
    scores = scores %>%
    dplyr::group_by(score_id, parname, type) %>%
    dplyr::mutate(is.outlier=score %in% boxplot.stats(score)$out) %>%
    dplyr::filter(!is.outlier)
  
  scores_nmf = scores %>%
    #dplyr::filter(parname == "K") %>%
    #dplyr::filter(type %in% types) %>%
    ggplot(aes(x=as.integer(value), y=score)) +
    geom_point(aes(shape=label), size=5) +
    geom_line(linetype="solid") +
    
    #ggrepel::geom_label_repel(aes(label=label), box.padding=0.05, size=3,
    #                          na.rm=T, show.legend=F) +
    
    ggh4x::facet_nested_wrap(type ~ ., scales="free", nrow=length(types)) +
    theme_bw() + labs(title="Scores", x="No. of Denovo Signatures", y="BIC") +
    guides(color=guide_legend(title="Seed")) +
    scale_x_continuous(breaks=scores %>% dplyr::filter(parname == "K") %>% dplyr::pull(value) %>% unique()) + 
    theme(legend.title = element_blank(), legend.position = "none")
  
  return(scores_nmf)
  
}


#-------------------------------------------------------------------------------


library(tidyverse)


plot_centroids_etiology <- function(
    x, 
    sbs_aetiology_path, 
    dbs_aetiology_path, 
    exposure_thr = 0.05, 
    quantile_thr = 0.9, 
    custom_theme = NULL
    #manual=F, 
    #clusters=get_cluster_labels(x), 
    #sig_cls=NULL, 
    #types=get_types(x)
  ) {
  
  aetiology <- rbind(
    read.table(file = sbs_aetiology_path, sep = ',', header = TRUE), 
    read.table(file = dbs_aetiology_path, sep = ',', header = TRUE)
  )
  
    
  scores <- get_clusters_score(x, get_types(x), exposure_thr, quantile_thr) %>% subset(significance == T)
  
  # add denovo signatures to aetiology dataframe as UNKNOWN aetiology
  #(scores %>% dplyr::pull(signature) %>% unique) %in% (aetiology %>% dplyr::pull(signature))
  dn_sigs <- setdiff(
    (scores %>% dplyr::pull(signature) %>% unique), 
    (aetiology %>% dplyr::pull(signature))
  )
  aetiology <- rbind( aetiology, data.frame(signature = dn_sigs, aetiology = rep("UNKNOWN", length(dn_sigs))) )
  
  # quality control
  if ( !(all( (scores %>% dplyr::pull(signature) %>% unique) %in% (aetiology %>% dplyr::pull(signature)) ) == TRUE) ) {
    warning("not all signatures in aetiology")
  }
  
  # need tidyverse package to laod
  df <- list(aetiology, scores) %>% reduce(full_join, by='signature') %>% select(signature, aetiology, cluster) %>% na.omit()
  
  # Create a dataframe with counts of cluster and aetiology combinations
  heatmap_data <- df %>%
    group_by(cluster, aetiology) %>%
    summarise(count = n(), .groups = 'drop') %>%
    mutate(present = ifelse(count > 0, 1, 0))
  
  #a <- plot_centroids(x_mapped, get_types(x_mapped), cls = NULL, sort_by = NULL, exposure_thr = 0, quantile_thr = 0.9) + 
  #  theme(
  #    text=element_text(size=20), 
  #    legend.title=element_text(size=18)
  #  )
  
  a <- plot_centroids(
    x, types = get_types(x), cls = NULL, sort_by = NULL, exposure_thr = exposure_thr, quantile_thr = quantile_thr
  ) + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), plot.margin = unit(c(5.5, 5.5, 0, 5.5), "pt")) + labs(y = "Density") + custom_theme
  
  
  # Plot the heatmap
  b <- ggplot(heatmap_data, aes(x = cluster, y = aetiology)) +
    geom_tile(aes(fill = factor(present)), color = "white", width = 0.9, height = 0.9) +
    scale_fill_manual(values = c("0" = "white", "1" = "gray"), guide = "none") +
    theme_bw() +
    labs(x = "Clusters", y = "Aetiology") + # title = "Cluster vs. Aetiology"
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position="none", plot.margin = unit(c(5.5, 5.5, 0, 5.5), "pt")) + custom_theme
  
  design <- "aaa
             aaa
             bbb"
  
  p <- patchwork::wrap_plots(a, b, ncol = 1, design = design)
  
  return(p)
}


#-------------------------------------------------------------------------------


# plot centroids by giving the significant signatures in each cluster
custom_centroid_plot <- function(x, clusters, sig_cls, types, sbs_aetiology_path, dbs_aetiology_path, custom_theme = NULL) {
  
  # CENTROIDS
  #-----------------------------------------------------------------------------
  df <- setNames(data.frame(matrix(ncol = 2, nrow = 0)), c("clusters", "sigs"))
  for (cluster in clusters) {
    d <- data.frame( samples=rep(cluster, length(sig_cls[[cluster]])) , sigs=sig_cls[[cluster]] )
    df <- rbind(df, d)
  }
  df <- df %>% dplyr::mutate(type=basilica:::match_type(types, sigs)) %>% as_tibble
  
  
  centr = basilica:::get_centroids(x, matrix = F)
  a_pr00 = centr %>%
    dplyr::mutate(type=basilica:::match_type(types, sigs)) %>%
    dplyr::filter(!is.na(type)) %>%
    dplyr::rename(samples=clusters)
  a_pr0 <- df %>% dplyr::left_join(a_pr00, by = c('samples', 'sigs', 'type'))
  
  
  others <- aggregate(value ~ samples + type, a_pr0, sum) %>% 
    mutate(other=1-value) %>% 
    select(samples, type, other) %>% 
    dplyr::rename("value"="other") %>% 
    mutate(sigs="Others")
  
  a_pr <- rbind(a_pr0, others)
  
  # signatures colors
  signames <- lapply(types, function(tid) { (a_pr %>% subset(type==tid) %>% dplyr::pull(sigs) %>% unique) })
  
  unique_signames <- signames %>% unlist %>% unique
  unique_signames <- unique_signames[! unique_signames %in% c('Others')]
  
  names(signames) <- types
  cls <- basilica:::gen_palette_aux(signames)
  cls["Others"] <- "gainsboro"
  
  
  a <- plot_exposures_aux(exposures=a_pr,
                          cls=cls,
                          titlee="Centroids",
                          sample_name=TRUE,
                          sample_levels=NULL,
                          sigs_levels=NULL) +
    scale_fill_manual(values=cls) + theme(axis.text.x=element_text(angle=0)) + 
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), plot.margin = unit(c(5.5, 5.5, 0, 5.5), "pt")) + labs(y = "Density") + custom_theme
  
  
  # ETIOLOGY
  #-----------------------------------------------------------------------------
  aetiology <- rbind(
    read.table(file = sbs_aetiology_path, sep = ',', header = TRUE), 
    read.table(file = dbs_aetiology_path, sep = ',', header = TRUE)
  )
  
  # check if all signatures are present in etiology signatures list
  dn_sigs <- setdiff(
    unique_signames, 
    (aetiology %>% dplyr::pull(signature))
  )
  # if not add them to etiology list as UNKNOWN
  aetiology <- rbind( aetiology, data.frame(signature = dn_sigs, aetiology = rep("UNKNOWN", length(dn_sigs))) )
  
  
  # need tidyverse package to laod
  a_pr <- a_pr %>% dplyr::rename("signature"="sigs", "cluster"="samples")
  df <- list(aetiology, a_pr) %>% reduce(full_join, by='signature') %>% select(signature, aetiology, cluster) %>% na.omit()
  
  # Create a dataframe with counts of cluster and aetiology combinations
  heatmap_data <- df %>%
    group_by(cluster, aetiology) %>%
    summarise(count = n(), .groups = 'drop') %>%
    mutate(present = ifelse(count > 0, 1, 0))
  
  
  #a_pr <- a_pr %>% dplyr::rename("clusters"="cluster", "sigs"="signature")
  
  # Plot the heatmap
  b <- ggplot(heatmap_data, aes(x = cluster, y = aetiology)) +
    geom_tile(aes(fill = factor(present)), color = "white", width = 0.9, height = 0.9) +
    scale_fill_manual(values = c("0" = "white", "1" = "gray"), guide = "none") +
    theme_bw() +
    labs(x = "Clusters", y = "Aetiology") + # title = "Cluster vs. Aetiology"
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position="none", plot.margin = unit(c(5.5, 5.5, 0, 5.5), "pt")) + custom_theme
  
  design <- "aaa
             aaa
             bbb"
  
  p <- patchwork::wrap_plots(a, b, ncol = 1, design = design)
  
  return(p)
}


#-------------------------------------------------------------------------------


plot_custom_exposures = function(
    x,
    types=get_types(x),
    samples=get_samples(x),
    clusters=get_cluster_labels(x),
    sample_name=FALSE,
    color_palette=NULL,
    sort_by=NULL, 
    signature_list
) {
  
  exposures = lapply(types, function(t)
    get_exposure(x, types=types, samples=samples,
                 clusters=clusters, add_groups=TRUE)[[t]] %>%
      dplyr::mutate(type=t)) %>%
    do.call(rbind, .)
  
  if (is.null(color_palette)) cls = gen_palette(x, types=sort(types))
  
  # merging signatures where their exposure in all the samples are below the threshold
  to_keep = exposures$sigs %>% unique(); sigs_levels = NULL
  
  exp <- exposures %>% dplyr::mutate(sigs=ifelse( (sigs %in% signature_list), sigs, "Others"))
  
  cls["Others"] <- "gainsboro"
  
  
  p <- plot_exposures_aux(exposures=exp, cls=cls,
                          titlee="Exposure",
                          sample_name=FALSE,
                          sort_by=NULL,
                          sigs_levels=sigs_levels) +
    scale_fill_manual(values=cls, breaks=to_keep)
  
  return(p)
}






