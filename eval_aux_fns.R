compute_quantiles = function(all_stats, colname) {
  qq = all_stats[[colname]] %>% quantile(c(0.33,0.66,1.))
  
  all_stats %>% 
    dplyr::mutate("{colname}_cat":=dplyr::case_when(
    .data[[colname]] <= qq[[1]] ~ paste0(0,"-",qq[[1]]),
    .data[[colname]] > qq[[1]] & .data[[colname]] <= qq[[2]] ~ paste0(qq[[1]]+1,"-",qq[[2]]),
    .data[[colname]] > qq[[2]] & .data[[colname]] <= qq[[3]] ~ paste0(qq[[2]]+1,"-",qq[[3]])
  ))
}


## function to remove boxplot outliers
filter_lims = function(x){
  l = boxplot.stats(x)$stats[1]
  u = boxplot.stats(x)$stats[5]
  x_tmp = x
  for (i in 1:length(x)) x_tmp[i] = ifelse(x[i]>l & x[i]<u, x[i], NA)
  if (all(is.na(x_tmp))) return(x)
  return(x_tmp)
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
    eval_single_fit_matched(fits[[ff]], x.simul) %>%
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
  assigned_missing_all = get_assigned_missing(x.fit=x.fit, x.simul=x.simul, cutoff=cutoff)
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
    
    ari_nmi = ari_nmi_km = ari_nmi_km_em = list(NA, NA)
    if (have_groups(x.fit)) ari_nmi = compute_ari_nmi(x.simul=x.simul, x.fit=x.fit)
    
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
      "type"=tid
    )
    
    # colnames(res) = paste(colnames(res), tid, sep="_")
    return(res)
  }) %>% setNames(get_types(x.fit))
}


