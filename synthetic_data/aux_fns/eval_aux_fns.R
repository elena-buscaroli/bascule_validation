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
    eval_single_fit_matched(x.fit=fits[[ff]], x.simul=x.simul) %>%
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



eval_single_fit_matched = function(x.fit, x.simul, cutoff=0.8) {
  x.fit = x.fit %>% rename_dn_expos()
  assigned_missing_all = get_assigned_missing(x=x.fit, x.simul=x.simul, cutoff=cutoff)
  
  # clustering stuff ####
  ari_nmi = ari_nmi_KM = list(NA, NA)
  if (have_groups(x.fit)) {
    ari_nmi = compute_ari_nmi(groups_simul=get_cluster_assignments(x.simul) %>% dplyr::arrange(samples) %>% dplyr::pull(clusters), 
                              groups_fit=get_cluster_assignments(x.fit) %>% dplyr::arrange(samples) %>% dplyr::pull(clusters))
    
    KM_groups = run_clustering(x.fit, method="kmeans")
    # KL.KM_groups = run_clustering(x.fit, method="kl_kmeans")
    # JS.SC_groups = run_clustering(x.fit, method="js_spectral")
    ari_nmi_KM = compute_ari_nmi(groups_simul=get_cluster_assignments(x.simul) %>% dplyr::arrange(samples) %>% dplyr::pull(clusters), 
                                 groups_fit=KM_groups %>% dplyr::arrange(samples) %>% dplyr::pull(clusters))
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
run_clustering = function(x.fit, method) {
  max_g = x.fit$clustering$pyro$params$init_params$pi %>% length()
  expos = get_exposure(x.fit, matrix=T) %>% dplyr::bind_cols()

  if (method == "kmeans") {
    gap_stats = cluster::clusGap(expos, FUNcluster=kmeans, K.max=max_g, nstart=25)
    best_K = cluster::maxSE(gap_stats$Tab[, "gap"], gap_stats$Tab[, "SE.sim"], method="Tibs2001SEmax")
    res = kmeans(expos, centers=best_K, nstart=25)
  } else if (method == "kl_kmeans") {
    gap_stats = cluster::clusGap(expos, FUNcluster=kl_kmeans, K.max=max_g)
    best_K = cluster::maxSE(gap_stats$Tab[, "gap"], gap_stats$Tab[, "SE.sim"], method="Tibs2001SEmax")
    res = kl_kmeans(expos, best_K)
  } else if (method == "js_spectral") {
    gap_stats = cluster::clusGap(expos, FUNcluster=js_spectral, K.max=max_g)
    best_K = cluster::maxSE(gap_stats$Tab[, "gap"], gap_stats$Tab[, "SE.sim"], method="Tibs2001SEmax")
    js_spectral(expos, best_K)
  }
  
  return(tibble::tibble(samples=names(res$cluster), clusters=res$cluster))
  
  # gap_stats = cluster::clusGap(expos, FUNcluster=kmeans, K.max=max_g, nstart=25)
  # best_K = cluster::maxSE(gap_stats$Tab[, "gap"], gap_stats$Tab[, "SE.sim"], method="Tibs2001SEmax")
  # 
  # km = kmeans(expos, centers=best_K, nstart=25)
  # return(tibble::tibble(samples=names(km$cluster), clusters=km$cluster))
}

## KL clustering #####

# kl_dist = function(p, q) {
#   p = p / sum(p)
#   q = q / sum(q)
#   kl_div = sum(p * log(p / q), na.rm=TRUE)
#   return(kl_div)
# }

## KL distance matrix
# kl_distance_matrix = function(input_mat) {
#   n = nrow(input_mat)
#   dist_mat = matrix(0, n, n)
#   
#   for (i in 1:(n - 1)) {
#     for (j in (i + 1):n) {
#       dist_mat[i, j] = kl_dist(input_mat[i, ], input_mat[j, ])
#       dist_mat[j, i] = dist_mat[i, j]  # Symmetric matrix
#     }
#   }
#   as.dist(dist_mat)  # Convert matrix to dist object
# }

kl_distance_matrix = function(input_mat) {
  distance_mat = matrix(nrow=nrow(input_mat), ncol=nrow(input_mat))
  input_mat = as.matrix(input_mat)
  for (i in 1:nrow(input_mat)) {
    for (j in 1:nrow(input_mat)) {
      distance_mat[i,j] = seewave::kl.dist(input_mat[i,], input_mat[j,], base=2)$D
    }
  }
  rownames(distance_mat) = colnames(distance_mat) = rownames(input_mat)
  return(distance_mat)
} 

# Custom K-Means with KL divergence
kl_kmeans = function(input_mat, k) {
  kl_dist_mat = kl_distance_matrix(input_mat)
  res = cluster::pam(kl_dist_mat, k=k)
  return(res)
}


## JS divergence spectral clustering ####

library(kernlab)  # For Spectral Clustering
library(proxy)    # For custom distance functions

# Define Jensen-Shannon distance
js_divergence = function(p, q) {
  m = (p + q) / 2
  kl_p_m = sum(p * log(p / (m + 1e-10)))
  kl_q_m = sum(q * log(q / (m + 1e-10)))
  0.5 * (kl_p_m + kl_q_m)
}

# Compute pairwise Jensen-Shannon distances
js_dist_matrix = function(input_mat) {
  as.matrix(proxy::dist(input_mat, method=js_divergence))
}

# Convert to similarity matrix
similarity_matrix = function(dist_matrix) {
  sigma = mean(dist_matrix)  # Scale factor
  exp(-dist_matrix^2 / (2 * sigma^2))
}

# Spectral clustering function
js_spectral = function(input_mat, k) {
  dist_matrix = js_dist_matrix(input_mat)
  sim_matrix = similarity_matrix(dist_matrix)
  clusters = specc(as.kernelMatrix(sim_matrix), centers=k)
  return(clusters)
}

