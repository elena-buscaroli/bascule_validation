plot_K = function(all_stats, fill="", facet="") {
  list_cols = lapply(colnames(all_stats), 
                     function(x) if (is.list(all_stats[[x]])) return(x)) %>% unlist()
  id_cols = c("N","G","seed","idd","fname","type","penalty")
  
  vln_fn = geom_violin(aes(x=factor(N), y=value)); pal = c()
  if (fill != "") {
    vln_fn = geom_violin(aes(x=factor(N), y=value, fill=get(fill), color=get(fill)))
    pal = c("tan2","#8FBC8B","thistle2")
  }
  
  if (facet=="") facet = "type ~ metric"
  
  all_stats %>% 
    compute_quantiles(colname="K_true") %>% 
    dplyr::filter(penalty=="Autoguide") %>% 
    dplyr::select(-dplyr::all_of(list_cols)) %>% 
    reshape2::melt(id=c(id_cols,"K_true_cat"), variable.name="metric") %>% 
    dplyr::filter(grepl("_ratio$", metric)) %>% 
    ggplot() +
    vln_fn +
    geom_hline(yintercept=1, color="red4", linetype="dashed") +
    facet_grid(as.formula(facet), scales="free_y") + 
    scale_fill_manual(values=pal) +
    scale_color_manual(values=pal) +
    theme_bw() + ylim(0,NA)
}

