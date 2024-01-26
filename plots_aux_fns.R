# Plot aux fn #####
plt = function(df, fn, pal, facet, grps_cols=NULL) {
  p = df %>% 
    ggplot() +
    fn +
    scale_fill_manual(values=pal) +
    scale_color_manual(values=pal) +
    ggh4x::facet_nested(as.formula(facet), scales="free_y") +
    theme_bw()
  if (!is.null(grps_cols)) p + coord_cartesian(ylim=get_l_u(df, grps_cols))
  else p
}


# N of signatures #####
plot_K = function(all_stats, fill="", facet="type ~ metric") {
  list_cols = get_colnames_islist(all_stats)
  id_cols = get_id_cols(all_stats)
  grps_cols = get_grps_cols(fill, facet)
  
  vln_fn = geom_violin(aes(x=factor(N), y=value), draw_quantiles=c(.5), 
                       position=position_dodge(width=1.2))
  pal = c()
  if (fill != "") {
    vln_fn = geom_violin(aes(x=factor(N), y=value, fill=get(fill)), 
                         draw_quantiles=c(.5), lwd=.3, 
                         position=position_dodge(width=.7))
    pal = c("tan2","#8FBC8B","thistle2")
  }
  
  p = all_stats %>% 
    dplyr::select(-dplyr::all_of(list_cols)) %>% 
    reshape2::melt(id=id_cols, variable.name="metric") %>% 
    dplyr::filter(grepl("_ratio$", metric)) %>%
    plt(fn=vln_fn, pal=pal, facet=facet, grps_cols=grps_cols) +
    geom_hline(yintercept=1, color="grey60", linetype="dashed") +
    ylim(0,NA)
  
  if (fill=="") p + geom_jitter(aes(x=factor(N), y=value), size=.5, alpha=.7)
  else p
}


# Metrics #####
plot_performance = function(all_stats, fill="", facet="type ~ variable") {
  list_cols = get_colnames_islist(all_stats)
  id_cols = get_id_cols(all_stats)
  metrics = get_metrics()
  grps_cols = get_grps_cols(fill, facet)
  
  bxplt_fn = geom_boxplot(aes(x=factor(N), y=value), outlier.shape=NA); pal = c()
  if (fill != "") {
    bxplt_fn = geom_boxplot(aes(x=factor(N), y=value, fill=get(fill)), outlier.shape=NA, lwd=0.3)
    pal = c("tan2","#8FBC8B","thistle2")
  }

  cosine = all_stats %>% 
    dplyr::select(-dplyr::all_of(list_cols)) %>% 
    reshape2::melt(id=id_cols, variable.name="metric") %>% 
    dplyr::filter(grepl("^cosine", metric)) %>% 
    tidyr::separate("metric", into=c("metric","variable"), extra="merge", sep="_") %>% 
    plt(fn=bxplt_fn, pal=pal, facet=facet, grps_cols=grps_cols) +
    ylab("Cosine similarity")
  
  mse = all_stats %>% 
    dplyr::select(-dplyr::all_of(list_cols)) %>% 
    reshape2::melt(id=id_cols, variable.name="metric") %>% 
    dplyr::filter(grepl("^mse", metric)) %>% 
    tidyr::separate("metric", into=c("metric","variable"), extra="merge", sep="_") %>% 
    plt(fn=bxplt_fn, pal=pal, facet=facet, grps_cols=grps_cols) +
    theme_bw() + ylab("MSE")
  
  patchwork::wrap_plots(mse, cosine, guides="collect") & theme(legend.position="bottom")
}


# Clustering #####
plot_performance_clustering = function(all_stats, fill="penalty", facet="~metric") {
  list_cols = get_colnames_islist(all_stats)
  id_cols = get_id_cols(all_stats)
  grps_cols = get_grps_cols(fill, facet)
  
  bxplt_fn = geom_boxplot(aes(x=factor(N), y=value, fill=get(fill)), outlier.shape=NA)
  pal = c("tan2","steelblue2","forestgreen")

  all_stats_sub = all_stats %>% 
    dplyr::select(-dplyr::all_of(list_cols)) %>%
    reshape2::melt(id=id_cols, variable.name="metric") %>% 
    dplyr::select(-type) %>% unique() %>% 
    dplyr::filter(metric %in% c("ari","nmi"))
  
  ylim = get_l_u(all_stats_sub, grps_cols)
    
  all_stats_sub %>% 
    plt(fn=bxplt_fn, pal=pal, facet=facet, grps_cols=grps_cols) +
    geom_hline(yintercept=1, color="grey70", linetype="dashed") +
    theme(legend.position="bottom") + guides(fill=guide_legend(title=fill))
}


# Getters and aux #####
## function to remove boxplot outliers
filter_lims = function(x) {
  l = boxplot.stats(x)$stats[1]
  u = boxplot.stats(x)$stats[5]
  x_tmp = replace(x, x<l | x>u, NA)
  if (all(is.na(x_tmp))) return(x)
  return(x_tmp)
}


get_l_u = function(all_stats, grps_cols) {
  l_u = all_stats %>% 
    dplyr::group_by(dplyr::pick(dplyr::all_of(grps_cols))) %>% 
    dplyr::reframe(l=boxplot.stats(value)$stats[1],
                   u=boxplot.stats(value)$stats[5])
  return(c(min(l_u$l), max(l_u$u)))
}


get_colnames_islist = function(all_stats) {
  lapply(colnames(all_stats), 
         function(x) if (is.list(all_stats[[x]])) return(x)) %>% unlist()
}


get_id_cols = function(all_stats) {
  ids = c("N","G","seed","idd","fname","type","penalty")
  if ("K_true_cat" %in% colnames(all_stats)) return(c(ids, "K_true_cat"))
  return(ids)
}


get_metrics = function() {
  c("mse_counts","cosine_expos","cosine_expos_missing","cosine_sigs")
}


get_grps_cols = function(fill, facet) {
  c("N",fill,strsplit(facet,"~| |[+]")[[1]]) %>% purrr::discard(.p=function(x) x=="")
}


