compute_quantiles = function(all_stats, colname) {
  qq=all_stats[[colname]] %>% quantile(c(0.33,0.66,1.))
  
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


make_plots_stats=function(all_stats) {
  stats_tmp = all_stats

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
  
  idd = strsplit(fname, "/")[[1]]; idd=idd[[length(idd)]]
  
  stats = lapply(names(fits), function(ff) {
    print(ff)
    eval_single_fit_matched(x.fit=fits[[ff]], x.simul=x.simul, fname=fname) %>%
      dplyr::bind_rows() %>% 
      dplyr::mutate(penalty=ff)
  }) %>% dplyr::bind_rows()
  
  source = nMuts = NULL
  idd_tmp = stringr::str_remove_all(idd, ".Rds")
  idd_tmp = stringr::str_remove_all(idd_tmp, paste0(".", dataset_id))
  try({
    source = strsplit(idd_tmp, "[.]")[[1]][5]
    nMuts = strsplit(idd_tmp, "[.]")[[1]][6]
  })
  
  return(
    tibble::tibble(fname=fname,
                   "N"=(stringr::str_replace_all(idd, "N", "") %>% strsplit(split="[.]"))[[1]][2] %>%
                     as.numeric(),
                   "G"=(stringr::str_replace_all(idd, "G", "") %>% strsplit(split="[.]"))[[1]][3] %>%
                     as.numeric(),
                   "seed"=(stringr::str_replace_all(idd, "s", "") %>% strsplit(split="[.]"))[[1]][4] %>%
                     as.numeric(),
                   "source"=source,
                   "nMuts"=nMuts,
                   "idd"=idd) %>%
      dplyr::select(N, G, seed, idd, dplyr::everything()) %>%
      dplyr::bind_cols(stats)
  )
}


eval_single_fit_matched = function(x.fit, x.simul, fname=NULL, cutoff=0.8) {
  add_unassigned = rep(FALSE, length.out=length(get_types(x.fit))) %>% setNames(get_types(x.fit))
  for (tid in get_types(x.fit)) {
    if("unassigned" %in% colnames(get_exposure(x.fit, matrix=T)[[tid]])) {
      x.fit$nmf[[tid]]$exposure=x.fit$nmf[[tid]]$exposure %>% dplyr::filter(sigs!="unassigned")
      add_unassigned[[tid]]=TRUE
    }
  }
  
  x.fit = x.fit %>% rename_dn_expos()
  assigned_missing_all = get_assigned_missing(x=x.fit, x.simul=x.simul, cutoff=cutoff)
  
  ari_nmi = ari_nmi_KM = ari_nmi_KL = ari_nmi_JS = list(NA, NA)
  clustering_fits = NULL
  sil_ch_TRUE = sil_ch_DP = sil_ch_KM = sil_ch_KL = sil_ch_JS = list(NA, NA)
  if (have_groups(x.fit)) {
    ari_nmi = compute_ari_nmi(groups_simul=get_cluster_assignments(x.simul) %>% dplyr::arrange(samples) %>% dplyr::pull(clusters), 
                              groups_fit=get_cluster_assignments(x.fit) %>% dplyr::arrange(samples) %>% dplyr::pull(clusters))
    
    clustering_fname = fname %>% stringr::str_replace_all("simul_fit", "clustering")
    
    KM_groups = KL.KM_groups = JS.spect_groups = NULL
    if (is.null(KM_groups)) {
      KM_groups = run_clustering(x.fit, method="kmeans", B=10)
      cat("Kmeans done.\n")
    }
    if (is.null(KL.KM_groups)) {
      KL.KM_groups = run_clustering(x.fit, method="kl_kmeans", B=10)
      cat("KL-Kmeans done.\n")
    }
    if (is.null(JS.spect_groups)) {
      JS.spect_groups = run_clustering(x.fit, method="js_spectral", B=10)
      cat("Spectral clustering done.\n")
    }
      
    clustering_fits = list(KMeans=KM_groups,
                           KL_KMeans=KL.KM_groups,
                           JS_spectral=JS.spect_groups)
    
    ari_nmi_KM = compute_ari_nmi(groups_simul=get_cluster_assignments(x.simul) %>% dplyr::arrange(samples) %>% dplyr::pull(clusters), 
                                 groups_fit=KM_groups %>% dplyr::arrange(samples) %>% dplyr::pull(clusters))
    ari_nmi_KL = compute_ari_nmi(groups_simul=get_cluster_assignments(x.simul) %>% dplyr::arrange(samples) %>% dplyr::pull(clusters), 
                                 groups_fit=KL.KM_groups %>% dplyr::arrange(samples) %>% dplyr::pull(clusters))
    ari_nmi_JS = compute_ari_nmi(groups_simul=get_cluster_assignments(x.simul) %>% dplyr::arrange(samples) %>% dplyr::pull(clusters), 
                                 groups_fit=JS.spect_groups %>% dplyr::arrange(samples) %>% dplyr::pull(clusters))
    
    expos = get_exposure(x.fit, matrix=T) %>% dplyr::bind_cols()
    d_matrix = js_dist_matrix(as.matrix(expos))
    sil_ch_TRUE = compute_sil_ch(expos, d_matrix, x.simul$clustering$clusters$clusters %>% as.factor() %>% as.numeric())
    sil_ch_DP = compute_sil_ch(expos, d_matrix, stringr::str_remove_all(x.fit$clustering$clusters$clusters, "G") %>% as.numeric())
    sil_ch_KM = compute_sil_ch(expos, d_matrix, KM_groups %>% dplyr::arrange(samples) %>% dplyr::pull(clusters))
    sil_ch_KL = compute_sil_ch(expos, d_matrix, KL.KM_groups %>% dplyr::arrange(samples) %>% dplyr::pull(clusters))
    sil_ch_JS = compute_sil_ch(expos, d_matrix, JS.spect_groups %>% dplyr::arrange(samples) %>% dplyr::pull(clusters))
  }

  lapply(get_types(x.fit), function(tid) {
    sigs.fit = get_signatures(x.fit, matrix=T)[[tid]]; sigs.simul=get_signatures(x.simul, matrix=T)[[tid]]
    sigs_fixed.fit = get_fixed_signatures(x.fit, matrix=T)[[tid]]
    sigs_dn.fit = get_denovo_signatures(x.fit, matrix=T)[[tid]]
    
    expos.fit = get_exposure(x.fit, matrix=T)[[tid]]; expos.simul=get_exposure(x.simul, matrix=T)[[tid]]
    
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
      "scores_TRUE"=list(sil_ch_TRUE),
      "scores_DP"=list(sil_ch_DP),
      "scores_KM"=list(sil_ch_KM),
      "scores_KL"=list(sil_ch_KL),
      "scores_JS"=list(sil_ch_JS),
      "clustering_fits"=list(clustering_fits),
      
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
    best_K = kmeans_bestK(expos, kmax=max_g, B=B)
    
    if (best_K > 1) {
      res_tmp = kmeans(expos, centers=best_K, nstart=25)
      fit_obj = res_tmp
    } else {
      res_tmp = list(cluster=rep(1, nrow(expos)) %>% setNames(rownames(expos)))
      fit_obj = NULL
    }

  } else if (method == "kl_kmeans") {
    best_K = kl_kmeans_bestK(input_mat=expos, kmax=max_g, B=B)
    
    if (best_K > 1) {
      res_tmp = kl_kmeans(expos, best_K)
      fit_obj = res_tmp$fit_obj
    } else {
      res_tmp = list(cluster=rep(1, nrow(expos)) %>% setNames(rownames(expos)))
      fit_obj = NULL
    }
    
  } else if (method == "js_spectral") {
    best_K = js_spectral_bestK(input_mat=expos, kmax=max_g, B=B)
    
    if (best_K > 1) {
      res_tmp = js_spectral(input_mat=expos, k=best_K)
      fit_obj = res_tmp$fit_obj
    } else {
      res_tmp = list(cluster=rep(1, nrow(expos)) %>% setNames(rownames(expos)))
      fit_obj = NULL
    }
  }
  
  return(tibble::tibble(samples=rownames(expos), clusters=res_tmp$cluster, obj=list(fit_obj)))
}


## Kmeans #####

kmeans_bestK = function(input_mat, kmax, B) {
  gap_stats = cluster::clusGap(input_mat, FUNcluster=kmeans, K.max=kmax, B=B)
  return(cluster::maxSE(gap_stats$Tab[, "gap"], gap_stats$Tab[, "SE.sim"], method="Tibs2001SEmax"))
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
  if (k == 1) 
    return(list(cluster=rep(1, length.out=nrow(input_mat)), fit_obj=NA))
  
  set.seed(k)
  res = flexclust::kcca(as.matrix(input_mat), k=k, family=kl_family)
  return(list(cluster=res@cluster, fit_obj=res))
}

kl_kmeans_bestK = function(input_mat, kmax, B) {
  gap_stats = gap_stat_custom(input_mat=input_mat, kmax=kmax, B=B,
                              cluster_fn=kl_kmeans, distance_fn=kl_family@dist)
  return(gap_stats$best_k)
  
  # gap_stats = cluster::clusGap(input_mat, FUNcluster=kl_kmeans, K.max=kmax, B=10, spaceH0="original")
  # return(cluster::maxSE(gap_stats$Tab[, "gap"], gap_stats$Tab[, "SE.sim"], method="Tibs2001SEmax"))
}


## JS divergence spectral clustering ####

library(kernlab)
library(proxy)

js_divergence = function(p, q) {
  m = (p + q) / 2
  kl_p_m = sum(p * log(p / (m + 1e-10), base=2))
  kl_q_m = sum(q * log(q / (m + 1e-10), base=2))
  0.5 * (kl_p_m + kl_q_m)
}

js_dist_matrix = function(input_mat) {
  as.matrix(proxy::dist(input_mat, method=js_divergence))
}

similarity_matrix = function(dist_matrix) {
  sigma = mean(dist_matrix)  # Scale factor
  exp(-dist_matrix^2 / (2 * sigma^2))
}

js_spectral = function(k, input_mat=NULL, sim_matrix=NULL) {
  N = max(nrow(input_mat), nrow(sim_matrix))
  if (k == 1)
    return(list(cluster=rep(1, length.out=N), fit_obj=NA))
  
  if (is.null(sim_matrix)) {
    d_matrix = js_dist_matrix(input_mat)
    sim_matrix = similarity_matrix(d_matrix)
  }
  
  set.seed(k)
  res = kernlab::specc(as.kernelMatrix(sim_matrix), centers=k)
  return(list(cluster=res@.Data, fit_obj=res))
}

js_spectral_bestK = function(input_mat, kmax, B) {
  gap_stats = gap_stat_custom(input_mat=input_mat, kmax=kmax, B=B,
                              cluster_fn=js_spectral, distance_fn=js_divergence)
  return(gap_stats$best_k)

  # gap_stats = cluster::clusGap(sim_matrix, FUNcluster=js_spectral, K.max=kmax, B=10, spaceH0="original")
  # return(cluster::maxSE(gap_stats$Tab[, "gap"], gap_stats$Tab[, "SE.sim"], method="Tibs2001SEmax"))
}



## Scores #####

compute_sil_ch = function(input_mat, d_matrix, labels) {
  list("sil"=silhouette_js(d_matrix, labels),
       "ch"=js_ch_index(input_mat, labels))
}


## JS Silhouette score #####

silhouette_js = function(d_matrix, labels) {
  if (length(unique(labels)) == 1) labels[1] = labels[1] + 1
  
  sil = cluster::silhouette(labels, d_matrix)
  return(mean(sil[, 3]))
}



## JS Calinski-Harabasz Index #####

js_ch_index = function(input_mat, labels) {
  
  if (length(unique(labels)) == 1) labels[1]=labels[1] + 1
  
  unique_labels = unique(labels)
  k = length(unique_labels)
  N = nrow(input_mat)
  
  global_centroid = colMeans(input_mat)
  
  cluster_centroids = lapply(unique_labels, function(label) {
    colMeans(input_mat[labels == label, , drop=FALSE])
  }) %>% setNames(unique_labels)

  bss_js = sum(sapply(unique_labels, function(label) {
    n_cluster = sum(labels == label)
    centroid = cluster_centroids[[as.character(label)]]
    n_cluster * js_divergence(centroid, global_centroid)^2
  }))
  
  wss_js = sum(sapply(1:N, function(i) {
    cluster_id = labels[i]
    js_divergence(input_mat[i, ], cluster_centroids[[as.character(cluster_id)]])^2
  }))
  
  js_ch_index = (bss_js / (k - 1)) / (wss_js / (N - k))
  
  return(js_ch_index)
}




## Custom gap statistics #####

gap_stat_custom = function(input_mat, kmax, cluster_fn, distance_fn, B=10, seed=123) {
  set.seed(seed)
  n = nrow(input_mat)
  d = ncol(input_mat)
  gap_values = numeric(kmax)
  log_wks_null = matrix(0, nrow=B, ncol=kmax)
  log_wks_obs = numeric(kmax)
  
  compute_dispersion = function(null_input_mat, cluster_labs, distance_fn) {
    total_disp = 0
    for (k in unique(cluster_labs)) {
      cluster_points = null_input_mat[cluster_labs == k, , drop=FALSE]
      if (nrow(cluster_points) > 1) {
        centroid = as.matrix(colMeans(cluster_points)) %>% t()
        dists = apply(cluster_points, 1, function(x) distance_fn(as.matrix(x) %>% t(), centroid))
        total_disp = total_disp + sum(dists)
      }
    }
    return(total_disp)
  }
  
  for (k in 1:kmax) {
    clusters_obs = cluster_fn(input_mat=input_mat, k=k)
    Wk_obs = compute_dispersion(input_mat, clusters_obs$cluster, distance_fn)
    log_wks_obs[k] = log(Wk_obs)
    
    for (b in 1:B) {
      set.seed(b+k)
      null_input_mat = t(apply(input_mat, 1, function(x) { dirmult::rdirichlet(n=1, alpha=x*10) + 1e-10 }))
      clusters_null = cluster_fn(input_mat=null_input_mat, k=k)
      Wk_null = compute_dispersion(null_input_mat, clusters_null$cluster, distance_fn)
      log_wks_null[b, k] = log(Wk_null)
      cat(".")
    }
    
    gap_values[k] = mean(log_wks_null[, k]) - log_wks_obs[k]
    cat("\n")
  }
  
  # 1-SE rule: pick the smallest k such that Gap(k) >= Gap(k+1) - s_{k+1}
  sdk = apply(log_wks_null, 2, sd) * sqrt(1 + 1/B)
  best_k = 1
  for (k in 1:(length(gap_values) - 1)) {
    if (gap_values[k] >= gap_values[k + 1] - sdk[k + 1]) {
      best_k = k
      break
    }
  }
  
  return(list(gap=gap_values,
              logWk_obs=log_wks_obs,
              logWk_null=log_wks_null,
              best_k=best_k))
}

