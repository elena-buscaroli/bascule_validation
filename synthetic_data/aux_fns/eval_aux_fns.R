compute_quantiles = function(all_stats, colname) {
  qq = all_stats[[colname]] %>% quantile(c(0.33,0.66,1.))
  
  all_stats %>% 
    dplyr::mutate("{colname}_cat":=dplyr::case_when(
    .data[[colname]] <= qq[[1]] ~ paste0(0,"-",qq[[1]]),
    .data[[colname]] > qq[[1]] & .data[[colname]] <= qq[[2]] ~ paste0(qq[[1]]+1,"-",qq[[2]]),
    .data[[colname]] > qq[[2]] & .data[[colname]] <= qq[[3]] ~ paste0(qq[[2]]+1,"-",qq[[3]])
  ))
}




make_boxplot = function(all_stats, colname) {
  all_stats %>% dplyr::select(N, G, dplyr::contains(colname)) %>%
    reshape2::melt(id=c("N","G"), variable.name="type") %>%
    dplyr::mutate(type=stringr::str_replace_all(type, paste0(colname,"_"),"")) %>%
    ggplot() +
    geom_boxplot(aes(x=as.factor(N), y=value)) +
    ggh4x::facet_nested(type ~ G) +
    theme_bw()
}



make_plots_stats = function(all_stats) {
  stats_tmp = all_stats
    # dplyr::select(N, G, seed, penalty, dplyr::contains("cosine_fixed")) %>%
    # dplyr::rowwise() %>%
    # dplyr::mutate(cosine_fixed_SBS=ifelse(length(cosine_fixed_SBS)>0,
    #                                       list(unlist(cosine_fixed_SBS)),
    #                                       cosine_fixed_SBS),
    #               cosine_fixed_DBS=ifelse(length(cosine_fixed_DBS)>0,
    #                                       list(unlist(cosine_fixed_DBS)),
    #                                       cosine_fixed_DBS))

  sim1 = make_boxplot(stats_tmp %>% tidyr::unnest(cosine_fixed_SBS), "cosine_fixed_SBS") + labs(title="cosine_fixed_SBS")
  sim2 = make_boxplot(stats_tmp %>% tidyr::unnest(cosine_fixed_DBS), "cosine_fixed_DBS") + labs(title="cosine_fixed_DBS")
  sim = patchwork::wrap_plots(sim1, sim2, ncol=2)

  ari = nmi = NULL
  cosine_expos = make_boxplot(all_stats, "cosine_expos") + labs(title="cosine_expos")
  cosine_sigs = make_boxplot(all_stats, "cosine_sigs") + labs(title="cosine_sigs")
  mse_counts = make_boxplot(all_stats, "mse_counts") + labs(title="mse_counts")
  try( {ari = make_boxplot(all_stats, "ari") + labs(title="ari")} )
  try( {nmi = make_boxplot(all_stats, "nmi") + labs(title="nmi")} )

  k_ratio = all_stats %>% dplyr::select(N, G, seed, dplyr::contains("K_")) %>%
    reshape2::melt(id=c("N","G","seed"), variable.name="type") %>%
    tidyr::separate(col="type", into=c("else","what","type")) %>% dplyr::mutate("else"=NULL) %>%
    tidyr::pivot_wider(names_from="what", values_from="value") %>%
    dplyr::mutate(value=found/true) %>%
    ggplot() +
    geom_violin(aes(x=as.factor(N), y=value)) +
    ggh4x::facet_nested(type ~ G) +
    labs(title="K_ratio") +
    theme_bw()

  if (!is.null(ari))
    return(patchwork::wrap_plots(mse_counts, cosine_expos, cosine_sigs,
                                 k_ratio, sim1, sim2, ari, nmi, ncol=4))

  return(
    patchwork::wrap_plots(mse_counts, cosine_expos, cosine_sigs, k_ratio, sim1, sim2, ncol=3)
  )
}



make_plots_stats_compare = function(all_stats, boxplot=TRUE) {
  
  boxplot_compare = function(all_stats, colname) {
    all_stats %>% dplyr::select(N, G, penalty, K_true_SBS, K_cat, dplyr::contains(colname)) %>%
      reshape2::melt(id=c("N","G","penalty","K_true_SBS","K_cat"), variable.name="type") %>%
      dplyr::mutate(type=stringr::str_replace_all(type, paste0(colname,"_"),"")) %>%
      dplyr::filter(type=="SBS") %>% 
      # dplyr::mutate(value=filter_lims(value)) %>% dplyr::filter(!is.na(value)) %>%
      ggplot() +
      geom_boxplot(aes(x=factor(N), y=value, fill=penalty)) +
      ggh4x::facet_nested( ~ K_cat) +
      scale_fill_manual(values=c("tan2","dodgerblue3","#ca472f")) +
      theme_bw()
  }
  
  all_stats = all_stats %>% compute_quantiles(colname="K_true")
  
  if (boxplot) fn = boxplot_compare else fn = line_compare
  
  cosine_expos = fn(all_stats, "cosine_expos") + labs(title="cosine_expos")
  cosine_expos_missing = fn(all_stats, "cosine_expos_missing") + 
    labs(title="cosine_expos_missing")
  mse_expos_missing = fn(all_stats, "mse_expos_missing") + 
    labs(title="mse_expos_missing")
  cosine_sigs = fn(all_stats, "cosine_sigs") + labs(title="cosine_sigs")
  mse_counts = fn(all_stats, "mse_counts") + labs(title="mse_counts")
  
  k_ratio = all_stats %>% 
    # dplyr::select(-K_assigned_SBS, -K_assigned_DBS) %>% 
    dplyr::select(N, G, seed, penalty, dplyr::contains("K_")) %>%
    reshape2::melt(id=c("N","G","seed","penalty"), variable.name="type") %>%
    tidyr::separate(col="type", into=c("else","what","type")) %>% dplyr::mutate("else"=NULL) %>%
    dplyr::filter(type=="SBS") %>% 
    tidyr::pivot_wider(names_from="what", values_from="value") %>%
    dplyr::mutate(value=found/true) %>%
    dplyr::mutate(value=filter_lims(value)) %>% dplyr::filter(!is.na(value)) %>% 
    ggplot() +
    
    geom_violin(aes(x=as.factor(true), y=value, fill=penalty), draw_quantiles=c(0.5)) +
    geom_hline(yintercept=1, lty="dotted", color="grey40") +
    ggh4x::facet_nested(~ N) +
    labs(title="K_ratio") + # scale_y_continuous(breaks=seq(0, 2, by=0.2)) +
    theme_bw()
  
  return(
    patchwork::wrap_plots(mse_counts, cosine_expos, cosine_expos_missing,
                          mse_expos_missing, cosine_sigs, k_ratio, 
                          ncol=2, guides="collect")
  )
}



stats_single_data = function(fname, names_fits=list("NoPenalty"="fit.0", "PenaltyN"="fit.N")) {
  cat(paste0(fname, "\n"))
  simul_fit = readRDS(fname)
  x.simul = simul_fit$dataset
  
  fits = lapply(names_fits, function(fitname) 
    simul_fit[[fitname]] %>% merge_clusters()) %>% 
    setNames(names(names_fits))
  
  idd = strsplit(fname, "/")[[1]]; idd = idd[[length(idd)]]
  
  stats = lapply(names(fits), function(ff) {
    print(ff)
    eval_single_fit_matched(x.fit=fits[[ff]], x.simul=x.simul, fname=fname) %>%
      dplyr::bind_rows() %>% 
      dplyr::mutate(penalty=ff)
  }) %>% dplyr::bind_rows()
  
  return(
    tibble::tibble(fname=fname,
                   "N"=(stringr::str_replace_all(idd, "N", "") %>% strsplit(split="[.]"))[[1]][2] %>%
                     as.numeric(),
                   "G"=(stringr::str_replace_all(idd, "G", "") %>% strsplit(split="[.]"))[[1]][3] %>%
                     as.numeric(),
                   "seed"=(stringr::str_replace_all(idd, "s", "") %>% strsplit(split="[.]"))[[1]][4] %>%
                     as.numeric(),
                   "idd"=idd) %>%
      dplyr::select(N, G, seed, idd, dplyr::everything()) %>%
      dplyr::bind_cols(stats)
  )
}



eval_single_fit_matched = function(x.fit, x.simul, fname=NULL, cutoff=0.8) {
  add_unassigned = rep(FALSE, length.out=length(get_types(x.fit))) %>% setNames(get_types(x.fit))
  for (tid in get_types(x.fit)) {
    if("unassigned" %in% colnames(get_exposure(x.fit, matrix=T)[[tid]])) {
      x.fit$nmf[[tid]]$exposure = x.fit$nmf[[tid]]$exposure %>% dplyr::filter(sigs!="unassigned")
      add_unassigned[[tid]] = TRUE
    }
  }
  
  x.fit = x.fit %>% rename_dn_expos()
  assigned_missing_all = get_assigned_missing(x=x.fit, x.simul=x.simul, cutoff=cutoff)
  
  # for (tid in get_types(x.fit)) {
  #   expos = get_exposure(x.fit, matrix=T)[[tid]]
  #   if (add_unassigned[[tid]]) {
  #     print("INSIDE IF")
  #     expos$unassigned = 1 - rowSums(expos)
  #     x.fit$nmf[[tid]]$exposure = expos %>% wide_to_long(what="exposures")
  #   }
  # }


  ari_nmi = ari_nmi_KM = ari_nmi_KL = ari_nmi_JS = list(NA, NA)
  if (have_groups(x.fit)) {
    ari_nmi = compute_ari_nmi(groups_simul=get_cluster_assignments(x.simul) %>% dplyr::arrange(samples) %>% dplyr::pull(clusters), 
                              groups_fit=get_cluster_assignments(x.fit) %>% dplyr::arrange(samples) %>% dplyr::pull(clusters))
    
    clustering_fname = fname %>% stringr::str_replace_all("simul_fit", "clustering")
    
    KM_groups = KL.KM_groups = JS.spect_groups = NULL
    # if (!is.null(fname) & file.exists(clustering_fname)) {
    #   clustering_fits = readRDS(clustering_fname)
    #   KM_groups = clustering_fits$KMeans
    #   KL.KM_groups = clustering_fits$KL_KMeans
    #   JS.spect_groups = clustering_fits$JS_spectral
    # }
    if (is.null(KM_groups)) {
      KM_groups = run_clustering(x.fit, method="kmeans")
      cat("Kmeans done.\n")
    }
    if (is.null(KL.KM_groups)) {
      KL.KM_groups = run_clustering(x.fit, method="kl_kmeans")
      cat("KL-Kmeans done.\n")
    }
    if (is.null(JS.spect_groups)) {
      JS.spect_groups = run_clustering(x.fit, method="js_spectral")
      cat("Spectral clustering done.\n")
    }
      
    clustering_fits = list(KMeans=KM_groups,
                           KL_KMeans=KL.KM_groups,
                           JS_spectral=JS.spect_groups)
    
    # if (!is.null(fname)) saveRDS(clustering_fits, file=clustering_fname)
    
    ari_nmi_KM = compute_ari_nmi(groups_simul=get_cluster_assignments(x.simul) %>% dplyr::arrange(samples) %>% dplyr::pull(clusters), 
                                 groups_fit=KM_groups %>% dplyr::arrange(samples) %>% dplyr::pull(clusters))
    ari_nmi_KL = compute_ari_nmi(groups_simul=get_cluster_assignments(x.simul) %>% dplyr::arrange(samples) %>% dplyr::pull(clusters), 
                                 groups_fit=KL.KM_groups %>% dplyr::arrange(samples) %>% dplyr::pull(clusters))
    ari_nmi_JS = compute_ari_nmi(groups_simul=get_cluster_assignments(x.simul) %>% dplyr::arrange(samples) %>% dplyr::pull(clusters), 
                                 groups_fit=JS.spect_groups %>% dplyr::arrange(samples) %>% dplyr::pull(clusters))
  }

  lapply(get_types(x.fit), function(tid) {
    sigs.fit = get_signatures(x.fit, matrix=T)[[tid]]; sigs.simul = get_signatures(x.simul, matrix=T)[[tid]]
    sigs_fixed.fit = get_fixed_signatures(x.fit, matrix=T)[[tid]]
    sigs_dn.fit = get_denovo_signatures(x.fit, matrix=T)[[tid]]
    
    expos.fit = get_exposure(x.fit, matrix=T)[[tid]]; expos.simul = get_exposure(x.simul, matrix=T)[[tid]]
    
    assigned_missing = assigned_missing_all[[tid]]
    assigned = assigned_missing$assigned_tp
    unassigned = c(assigned_missing$missing_fn, assigned_missing$added_fp)
    
    mse_counts = compute.mse(m_true=get_input(x.simul, matrix=T)[[tid]], 
                             m_inf=get_input(x.fit, reconstructed=T, matrix=T)[[tid]])
    mse_expos = compute.mse(m_true=expos.simul, m_inf=expos.fit,
                            assigned_missing=assigned_missing)
    mse_expos_missing = compute.mse(m_true=expos.simul, m_inf=expos.fit,
                                    assigned_missing=assigned_missing, keep_missing=T)
    
    cosine_sigs = compute.cosine(m_true=sigs.simul, m_inf=sigs.fit,
                                 assigned_missing=assigned_missing,
                                 what="sigs")
    cosine_expos = compute.cosine(m_true=expos.simul, m_inf=expos.fit,
                                  assigned_missing=assigned_missing,
                                  what="expos")
    cosine_expos_missing = compute.cosine(m_true=expos.simul, m_inf=expos.fit,
                                          assigned_missing=assigned_missing,
                                          what="expos", keep_missing=T)

    res = tibble::tibble(
      "assigned_missing"=list(assigned_missing),
      "input_sigs"=list(get_input_signames(x.fit)[[tid]]),
      "fixed_sigs"=list(get_fixed_signames(x.fit)[[tid]]),
      "dn_sigs"=list(get_denovo_signames(x.fit)[[tid]]),

      "K_true"=length(get_signames(x.simul)[[tid]]),
      "K_found"=length(get_signames(x.fit)[[tid]]),
      "K_assigned"=length(assigned),
      
      "mse_counts"=mse_counts,
      "mse_expos"=mse_expos,
      "mse_expos_missing"=mse_expos_missing,
      
      "cosine_sigs"=cosine_sigs,
      "cosine_expos"=cosine_expos,
      "cosine_expos_missing"=cosine_expos_missing,
  
      "groups_found"=length(get_cluster_labels(x.fit)),
      "ari"=ari_nmi[[1]],
      "nmi"=ari_nmi[[2]],
      "ari_KM"=ari_nmi_KM[[1]],
      "nmi_KM"=ari_nmi_KM[[2]],
      "ari_KL"=ari_nmi_KL[[1]],
      "nmi_KL"=ari_nmi_KL[[2]],
      "ari_JS"=ari_nmi_JS[[1]],
      "nmi_JS"=ari_nmi_JS[[2]],
      
      "type"=tid
    ) %>% 
      dplyr::rowwise() %>% 
      dplyr::mutate(K_input_found=length(intersect(assigned, input_sigs)),
                    K_input_true=length(input_sigs),
                    K_dn_found=length(intersect(assigned, dn_sigs)),
                    K_dn_true=K_true - length(input_sigs),

                    K_input_ratio=K_input_found / K_input_true,
                    K_dn_ratio=K_dn_found / K_dn_true,
                    K_ratio=K_assigned / K_true) %>% 
      dplyr::ungroup()
    
    return(res)
  }) %>% setNames(get_types(x.fit))
}




# Clustering ####

# method in "kmeans", "kl_kmeans" or "js_spectral"
run_clustering = function(x.fit, method, B=50) {
  max_g = x.fit$clustering$pyro$params$init_params$pi %>% length()
  expos = get_exposure(x.fit, matrix=T) %>% dplyr::bind_cols()

  if (method == "kmeans") {
    best_K = kmeans_bestK(expos, kmin=2, kmax=max_g)
    
    if (best_K > 1) {
      res_tmp = kmeans(expos, centers=best_K, nstart=25)
      fit_obj = res_tmp
    } else {
      res_tmp = list(cluster=rep(1, nrow(expos)) %>% setNames(rownames(expos)))
      fit_obj = NULL
    }

  } else if (method == "kl_kmeans") {
    best_K = kl_kmeans_bestK(expos, kmin=2, kmax=max_g)
    
    if (best_K > 1) {
      res_tmp = kl_kmeans(expos, best_K)
      fit_obj = res_tmp$fit_obj
    } else {
      res_tmp = list(cluster=rep(1, nrow(expos)) %>% setNames(rownames(expos)))
      fit_obj = NULL
    }
    
  } else if (method == "js_spectral") {
    dist_matrix = js_dist_matrix(expos)
    sim_matrix = similarity_matrix(dist_matrix)
    best_K = js_spectral_bestK(sim_matrix, kmin=2, kmax=max_g, js_dist_matrix=dist_matrix)
    
    if (best_K > 1) {
      res_tmp = js_spectral(sim_matrix, best_K)
      fit_obj = res_tmp$fit_obj
    } else {
      res_tmp = list(cluster=rep(1, nrow(expos)) %>% setNames(rownames(expos)))
      fit_obj = NULL
    }
    
  }
  
  return(tibble::tibble(samples=rownames(expos), clusters=res_tmp$cluster, obj=list(fit_obj)))
}


## Kmeans #####

kmeans_bestK = function(input_mat, kmin, kmax) {
  k_range = kmin:kmax
  sil_values = numeric(length(k_range))
  
  km_dist_matrix = dist(input_mat)
  
  for (i in seq_along(k_range)) {
    k_i = k_range[i]
    
    set.seed(i)
    kmeans_model = kmeans(as.matrix(input_mat), centers=k_i, nstart=25)
    sil = cluster::silhouette(kmeans_model$cluster, km_dist_matrix)
    sil_values[i] = mean(sil[, 3])
  }
  
  return(k_range[which.max(sil_values)])
}

## KL clustering #####

library(flexclust)

kl_distance_row = function(x, centers) {
  if (ncol(x) != ncol(centers)) 
    stop(sQuote("x"), " and ", sQuote("centers"), " must have the same number of columns")
  z = matrix(0, nrow=nrow(x), ncol=nrow(centers))
  for (j in 1:nrow(centers)) {
    for (nn in 1:nrow(x)) {
      z[nn, j] = seewave::kl.dist(x[nn, ], centers[j, ], base=2)$D
    }
  }
  z
}

kl_family = flexclust::kccaFamily(dist=kl_distance_row, cent=function(x) colMeans(x), name="KL_dist")

kl_kmeans = function(input_mat, k) {
  set.seed(k)
  res = flexclust::kcca(as.matrix(input_mat), k=k, family=kl_family)
  return(list(cluster=res@cluster, fit_obj=res))
}

kl_kmeans_bestK = function(input_mat, kmin, kmax) {
  k_range = kmin:kmax
  sil_values = numeric(length(k_range))
  
  kl_dist_matrix = kl_distance_row(as.matrix(input_mat), centers=as.matrix(input_mat))
  cat("KL distance matrix done\n")
  
  for (i in seq_along(k_range)) {
    k_i = k_range[i]
    
    kcca_model = kl_kmeans(input_mat, k=k_i) # flexclust::kcca(as.matrix(input_mat), k=k_i, family=kl_family)
    sil = cluster::silhouette(kcca_model$cluster, as.dist(kl_dist_matrix))
    sil_values[i] = mean(sil[, 3])
  }
  
  return(k_range[which.max(sil_values)])
}


## JS divergence spectral clustering ####

library(kernlab)
library(proxy)

js_divergence = function(p, q) {
  m = (p + q) / 2
  kl_p_m = sum(p * log(p / (m + 1e-10)))
  kl_q_m = sum(q * log(q / (m + 1e-10)))
  0.5 * (kl_p_m + kl_q_m)
}

js_dist_matrix = function(input_mat) {
  as.matrix(proxy::dist(input_mat, method=js_divergence))
}

similarity_matrix = function(dist_matrix) {
  sigma = mean(dist_matrix)  # Scale factor
  exp(-dist_matrix^2 / (2 * sigma^2))
}

js_spectral = function(sim_matrix, k) {
  set.seed(k)
  res = kernlab::specc(as.kernelMatrix(sim_matrix), centers=k)
  return(list(cluster=res@.Data, fit_obj=res))
}

js_spectral_bestK = function(sim_matrix, kmin, kmax, js_dist_matrix) {
  k_range = kmin:kmax
  sil_values = numeric(length(k_range))
  
  for (i in seq_along(k_range)) {
    k_i = k_range[i]

    sc_model = js_spectral(sim_matrix, k=k_i) # kernlab::specc(as.kernelMatrix(sim_matrix), centers=k_i)
    sil = cluster::silhouette(sc_model$cluster, as.dist(js_dist_matrix))
    sil_values[i] = mean(sil[, 3])
  }
  
  return(k_range[which.max(sil_values)])
}


