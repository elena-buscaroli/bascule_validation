make_boxplot = function(all_stats, colname) {
  all_stats %>% dplyr::select(N, G, dplyr::contains(colname)) %>%
    reshape2::melt(id=c("N","G"), variable.name="type") %>%
    dplyr::mutate(type=stringr::str_replace_all(type, paste0(colname,"_"),"")) %>%
    ggplot() +
    geom_boxplot(aes(x=as.factor(N), y=value)) +
    ggh4x::facet_nested(type ~ G) +
    theme_bw()
}


eval_single_fit_matched = function(x.fit, x.simul, cutoff=0.8) {
  x.fit = x.fit %>% rename_dn_expos()
  assigned_missing_all = get_assigned_missing(x.fit=x.fit, x.simul=x.simul, cutoff=cutoff)
  lapply(get_types(x.fit), function(tid) {
    print(tid)
    sigs.fit = get_signatures(x.fit, matrix=T)[[tid]]; sigs.simul = get_signatures(x.simul, matrix=T)[[tid]]
    sigs_fixed.fit = get_fixed_signatures(x.fit, matrix=T)[[tid]]
    sigs_dn.fit = get_denovo_signatures(x.fit, matrix=T)[[tid]]

    expos.fit = get_exposure(x.fit, matrix=T)[[tid]]; expos.simul = get_exposure(x.simul, matrix=T)[[tid]]

    assigned_missing = assigned_missing_all[[tid]]
    assigned = assigned_missing$assigned_tp
    unassigned = c(assigned_missing$missing_fn, assigned_missing$added_fp)

    mse_counts = compute.mse(m_true=get_input(x.simul, matrix=T)[[tid]], m_inf=get_input(x.fit, reconstructed=T, matrix=T)[[tid]])
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

    cosine_fixed = list()
    if (!is.null(sigs_dn.fit)) {
      cosine_fixed = lapply(rownames(sigs_dn.fit), function(sid) {
        cosine_sim = lsa::cosine(t(sigs.fit))[sid, rownames(sigs_fixed.fit)]
      }) %>% setNames(rownames(sigs_dn.fit))
    }

    ari_nmi = ari_nmi_km = ari_nmi_km_em = list(NA, NA)
    if (have_groups(x.fit)) {
      # kmeans1 = get_initial_object(x.fit, what="clustering")
      # kmeans2 = get_initial_object(x.fit, what="clustering") %>% run_kmeans()

      ari_nmi = compute_ari_nmi(x.simul=x.simul, x.fit=x.fit)
      # ari_nmi_km1 = compute_ari_nmi(x.simul=x.simul, x.fit=kmeans1)
      # ari_nmi_km2 = compute_ari_nmi(x.simul=x.simul, x.fit=kmeans2)
    }

    res = tibble::tibble(
      "assigned_missing"=list(assigned_missing),

      "sigs_found"=length(assigned),
      # "shared_found"=sum(assigned %in% rare_common$shared),
      # "private_found"=sum(assigned %in% c(rare_common$private_rare, rare_common$private_common)),
      "K_true"=length(get_signames(x.simul)[[tid]]),
      "K_found"=length(get_signames(x.fit)[[tid]]),
      "K_assigned"=length(assigned_missing_all[[tid]][["assigned_tp"]]),

      "mse_counts"=mse_counts,
      "mse_expos"=mse_expos,
      "mse_expos_missing"=mse_expos_missing,

      "cosine_sigs"=cosine_sigs,
      "cosine_expos"=cosine_expos,
      "cosine_expos_missing"=cosine_expos_missing,
      "cosine_fixed"=list(cosine_fixed),

      "groups_found"=length(get_cluster_labels(x.fit)),
      "ari"=ari_nmi[[1]],
      "nmi"=ari_nmi[[2]]
      # "ari_km1"=ari_nmi_km1[[1]],
      # "nmi_km1"=ari_nmi_km1[[2]],
      # "ari_km2"=ari_nmi_km2[[1]],
      # "nmi_km2"=ari_nmi_km2[[2]]
    )

    colnames(res) = paste(colnames(res), tid, sep="_")
    return(res)
  }) %>% setNames(get_types(x.fit))
}



make_plots_stats = function(stats) {
  stats_tmp = stats %>%
    dplyr::select(N, G, seed, penalty, dplyr::contains("cosine_fixed")) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(cosine_fixed_SBS=ifelse(length(cosine_fixed_SBS)>0,
                                          list(unlist(cosine_fixed_SBS)),
                                          cosine_fixed_SBS),
                  cosine_fixed_DBS=ifelse(length(cosine_fixed_DBS)>0,
                                          list(unlist(cosine_fixed_DBS)),
                                          cosine_fixed_DBS))
    # tidyr::unnest(cosine_fixed_SBS)

  sim1 = make_boxplot(stats_tmp %>% tidyr::unnest(cosine_fixed_SBS), "cosine_fixed_SBS") + labs(title="cosine_fixed_SBS")
  sim2 = make_boxplot(stats_tmp %>% tidyr::unnest(cosine_fixed_DBS), "cosine_fixed_DBS") + labs(title="cosine_fixed_DBS")
  sim = patchwork::wrap_plots(sim1, sim2, ncol=2)

  ari = nmi = NULL
  cosine_expos = make_boxplot(stats, "cosine_expos") + labs(title="cosine_expos")
  cosine_sigs = make_boxplot(stats, "cosine_sigs") + labs(title="cosine_sigs")
  mse_counts = make_boxplot(stats, "mse_counts") + labs(title="mse_counts")
  try( {ari = make_boxplot(stats, "ari") + labs(title="ari")} )
  try( {nmi = make_boxplot(stats, "nmi") + labs(title="nmi")} )

  k_ratio = stats %>% dplyr::select(N, G, seed, dplyr::contains("K_")) %>%
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
  qq = all_stats$K_true_SBS %>% quantile(c(0.33,0.66,1.))
  
  all_stats = all_stats %>% dplyr::mutate(K_cat=dplyr::case_when(
    K_true_SBS <= qq[[1]] ~ paste0(0,"-",qq[[1]]),
    K_true_SBS > qq[[1]] & K_true_SBS <= qq[[2]] ~ paste0(qq[[1]]+1,"-",qq[[2]]),
    K_true_SBS > qq[[2]] & K_true_SBS <= qq[[3]] ~ paste0(qq[[2]]+1,"-",qq[[3]])
  ))
  
  filter_lims = function(x){
    l = boxplot.stats(x)$stats[1]
    u = boxplot.stats(x)$stats[5]
    x_tmp = x
    for (i in 1:length(x)) x_tmp[i] = ifelse(x[i]>l & x[i]<u, x[i], NA)
    if (all(is.na(x_tmp))) return(x)
    return(x_tmp)
  }
  
  boxplot_compare = function(all_stats, colname) {
    all_stats %>% dplyr::select(N, G, penalty, K_true_SBS, K_cat, dplyr::contains(colname)) %>%
      reshape2::melt(id=c("N","G","penalty","K_true_SBS","K_cat"), variable.name="type") %>%
      dplyr::mutate(type=stringr::str_replace_all(type, paste0(colname,"_"),"")) %>%
      dplyr::filter(type=="SBS") %>% 
      dplyr::mutate(value=filter_lims(value)) %>% dplyr::filter(!is.na(value)) %>%
      ggplot() +
      # geom_boxplot(aes(x=factor(K_true_SBS), y=value, fill=penalty), na.rm=T) +
      # ggh4x::facet_nested(penalty ~ N) +
      
      geom_boxplot(aes(x=factor(N), y=value, fill=penalty)) +
      ggh4x::facet_nested( ~ K_cat) +
      scale_fill_manual(values=c("tan2","dodgerblue3","#ca472f")) +
      theme_bw()
  }
  
  line_compare = function(all_stats, colname) {
    all_stats %>% dplyr::select(N, G, penalty, K_true_SBS, dplyr::contains(colname)) %>%
      reshape2::melt(id=c("N","G","penalty","K_true_SBS"), variable.name="type") %>%
      dplyr::mutate(type=stringr::str_replace_all(type, paste0(colname,"_"),"")) %>%
      dplyr::filter(type=="SBS") %>% 
      dplyr::mutate(value=filter_lims(value)) %>% dplyr::filter(!is.na(value)) %>% 
      ggplot() +
      geom_jitter(aes(x=K_true_SBS, y=value, color=penalty),
                  na.rm=T, size=.7, height=0, width=0.2, alpha=0.5) +
      geom_smooth(aes(x=K_true_SBS, y=value, color=penalty),
                  na.rm=T, se=F) +
      ggh4x::facet_nested( ~ N) +
      scale_x_continuous(breaks=unique(all_stats$K_true_SBS)) +
      theme_bw()
  }
  
  if (boxplot) fn = boxplot_compare else fn = line_compare
  
  col_palette = ggsci::pal_simpsons()(length(unique(all_stats$penalty)))
  
  cosine_expos = fn(all_stats, "cosine_expos") + labs(title="cosine_expos")
  cosine_expos_missing = fn(all_stats, "cosine_expos_missing") + 
    labs(title="cosine_expos_missing")
  mse_expos_missing = fn(all_stats, "mse_expos_missing") + 
    labs(title="mse_expos_missing")
  cosine_sigs = fn(all_stats, "cosine_sigs") + labs(title="cosine_sigs")
  mse_counts = fn(all_stats, "mse_counts") + labs(title="mse_counts")
  
  k_ratio = all_stats %>% 
    dplyr::select(-K_assigned_SBS, -K_assigned_DBS) %>% 
    dplyr::select(N, G, seed, penalty, dplyr::contains("K_")) %>%
    reshape2::melt(id=c("N","G","seed","penalty"), variable.name="type") %>%
    tidyr::separate(col="type", into=c("else","what","type")) %>% dplyr::mutate("else"=NULL) %>%
    dplyr::filter(type=="SBS") %>% 
    tidyr::pivot_wider(names_from="what", values_from="value") %>%
    dplyr::mutate(value=found/true) %>%
    dplyr::mutate(value=filter_lims(value)) %>% dplyr::filter(!is.na(value)) %>% 
    ggplot() +
    
    # geom_jitter(aes(x=true, y=value, color=penalty), size=.7, height=0, width=0.2, alpha=0.5) +
    geom_smooth(aes(x=true, y=value, color=penalty), se=F) + 
    scale_x_continuous(breaks=unique(all_stats$K_true_SBS)) +
    
    # geom_violin(aes(x=as.factor(true), y=value, fill=penalty), draw_quantiles=c(0.5)) +
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
      dplyr::bind_cols() %>% dplyr::mutate(penalty=ff)
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

