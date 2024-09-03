library(ggplot2)
library(ggrepel)
library(ggpubr)
library(reshape2)

sig_cls_organ_all <<- list(
  breast = list(G0=c("SBS2", "SBS13", "DBS2", "DBS13", "DBS11"), 
                G1=c("SBS3", "DBS2", "DBS13"), 
                G10=c("SBS1", "SBS3", "SBS2", "SBS13", "SBSD11", "DBS11"), 
                G11=c("SBS1", "SBS3", "SBS2", "SBS13", "SBSD11", "DBS14", "DBS13"), 
                G13=c("SBS1", "SBS3", "SBS2", "SBS13", "SBSD11", "DBS13")),
  lung = list(G1=c("SBS4", "DBS2"), 
              G3=c("SBS31", "DBS5"), 
              G13=c("SBS17b", "DBS13"), 
              G14=c("SBS1", "SBS5", "DBS1"), 
              G0=c("SBS1", "SBS5", "DBS13"), 
              G12=c("SBS1", "SBS5", "DBS6"), 
              G8=c("SBS1", "SBS5", "DBS20")),
  colorectal = list(G1=c("SBS1", "SBS3", "SBS5", "SBS18", "SBS35", "SBS44", "SBSD8", "DBS25"), 
                    G3=c("SBS1", "SBS3", "SBS5", "SBS18", "SBS35", "SBS44", "SBSD8", "DBS5"), 
                    G6=c("SBS1", "SBS3", "SBS5", "SBS18", "SBS35", "SBS44", "SBSD8", "DBS8"), 
                    G12=c("SBS1", "SBS3", "SBS5", "SBS18", "SBS35", "SBS44", "SBSD8", "DBS13"), 
                    G10=c("SBS44", "SBSD7", "SBSD12", "DBS14"), 
                    G9=c("SBS10a", "DBS3"), 
                    G0=c("SBS3", "SBS18", "DBS2"), 
                    G5=c("SBS15", "SBSD7", "DBSD3"))
)


# input:
#  - x: exposure (wide format)
#  - type: basilica | serena
#  - exposure_thr: threshold to eliminate exposure values
# output:
#  - dataframe of signatures, their mean exposure, frequency and type
remove_sparsity = function(x, type, exposure_thr=0.01) {
  aaa = wide_to_long(x, what="exposures") %>% subset( value > exposure_thr )
  bbb = data.frame(
    mean=tapply(aaa$value, aaa$sigs, mean), 
    freq=tapply(aaa$value, aaa$sigs, length), 
    type=type
  ) %>% tibble::rownames_to_column("name")
  return(bbb)
}



library(dplyr)
# input
#   basilica_exposure --> wide dataframe
#   serena_exposure ----> wide dataframe
#   basilica_singles ---> unmapped basilica signatures (character vector)
#   serena_singles -----> unmapped serena signatures (character vector)
# output
#   ggplot object
plot_unmapped = function(basilica_exp, serena_exp, basilica_singles, serena_singles, exposure_thr=0.01) {
  
  if ( (basilica_singles %>% length) > 0 ) {
    a = remove_sparsity(
      ( basilica_exp %>% select( all_of(basilica_singles) ) ), 
      "Basilica", 
      exposure_thr
    )
  } else {
    a = NULL
  }
  
  if ( (serena_singles %>% length) > 0 ) {
    b = remove_sparsity(
      (serena_exp %>% select( all_of(serena_singles[serena_singles %in% colnames(serena_exp)]) )), 
      "Serena", 
      exposure_thr
    )
  } else {
    b = NULL
  }
  
  df = rbind(a, b)
  
  #df = df %>% replace(is.na(.), 0)
  
  p = ggplot(
    data=df, 
    aes(x=freq, y=mean, label=name)) + 
    geom_point(aes(color=factor(type)), size=2) + 
    #ggtitle(paste0("Unexplained Signatures (exposure > ", exposure_thr*100, "%)")) + 
    #xlab("Number of Samples") + 
    #ylab("Mean Exposure") + 
    labs(
      title="Unexplained Signatures Mean Exposures", 
      subtitle=paste0("(exposures > ", exposure_thr*100, "%)"), 
      caption="samples with exposure below the threshold are removed", 
      x="Number of Samples", 
      y="Mean Exposure", 
      color='Tools'
    ) + 
    geom_label_repel(
      aes(label=name),
      box.padding=0.35, 
      point.padding=0.5,
      segment.color='grey50', 
      max.overlaps=25
    ) +
    expand_limits(y=c(0,1)) + 
    theme_bw()
  
  return(p)
}



plot_exposure_comp = function(obj, map) {
  
  basilica_exp = obj$BasilicaExposure
  serena_exp = obj$SerenaExposure
  basilica_serena = map$basilica_serena
  
  p_list = list()
  for (i in 1:nrow(basilica_serena)) {
    df = data.frame(
      b_pair=basilica_exp[, basilica_serena[i,]$Basilica], 
      s_pair=serena_exp[, basilica_serena[i,]$Serena]
    )
    p_list[[length(p_list)+1]] = ggplot(df, aes(x=b_pair, y=s_pair)) + 
      geom_point() + 
      labs(
        title="", 
        x=paste0(basilica_serena[i,]$Basilica, " exposure (basilica)"), 
        y=paste0(basilica_serena[i,]$Serena, " exposure (serena)")
      ) + 
      #xlab(paste0(basilica_serena[i,]$Basilica, " exposure (basilica)")) +
      #ylab(paste0(basilica_serena[i,]$Serena, " exposure (serena)")) + 
      ylim(c(0,1)) + 
      theme_bw()
  }
  
  p = patchwork::wrap_plots(p_list) & patchwork::plot_annotation(
    title="Exposure Alignment of Mapped Signatures", 
    subtitle="(Basilica vs. Serena)"
    )
  
  return(p)
}




plot_mapping_barplot = function(x_orig, cls_catalogues,
                                degasperi_sbs="real_data/datasets/signatures_sbs.Rds",
                                degasperi_dbs="real_data/datasets/signatures_dbs.Rds") {
  degasperi <<- readRDS(degasperi_sbs) %>% 
    wide_to_long(what="beta") %>% dplyr::mutate(type="SBS") %>% 
    dplyr::bind_rows(
      readRDS(degasperi_dbs) %>% 
        wide_to_long(what="beta") %>% dplyr::mutate(type="DBS")
    )
  cosmic <<- COSMIC_sbs %>% 
    wide_to_long(what="beta") %>% dplyr::mutate(type="SBS") %>% 
    dplyr::bind_rows(
      COSMIC_dbs %>% 
        wide_to_long(what="beta") %>% dplyr::mutate(type="DBS")
    )
  
  reference_cat_cosmic = lapply(X=unique(cosmic$type), FUN=function(tid) 
    cosmic %>% dplyr::filter(type==tid) %>% dplyr::select(-type) %>% long_to_wide(what="beta")
  ) %>% setNames(unique(cosmic$type))
  
  # signatures in degasperi but not in cosmic
  keep_degasperi = setdiff(unique(degasperi$sigs), unique(cosmic$sigs))
  reference_cat_degasperi = lapply(X=unique(cosmic$type), FUN=function(tid) 
    degasperi %>% 
      # dplyr::filter(sigs%in%keep_degasperi) %>% 
      dplyr::filter(type==tid) %>% dplyr::select(-type) %>% 
      long_to_wide(what="beta")
  ) %>% setNames(unique(cosmic$type))
  
  assigned_missing_cosmic = get_assigned_missing(x_orig, reference_cat=reference_cat_cosmic, 
                                                 cutoff=0.8)
  x_mapped_cosmic = x_orig %>% convert_dn_names(reference_cat=reference_cat_cosmic, cutoff=0.8)
  
  assigned_missing_degasperi = get_assigned_missing(x_mapped_cosmic, reference_cat=reference_cat_degasperi, 
                                                    cutoff=0.)
  
  assigned_cosmic = lapply(names(assigned_missing_cosmic), function(tid) {
    assigned = assigned_missing_cosmic[[tid]]$assigned_tp
    data.frame(list(bascule=assigned, reference=names(assigned)), row.names=NULL) %>% 
      dplyr::mutate(type=tid, catalogue="COSMIC v3.4") %>% 
      dplyr::mutate(signature_type=ifelse(bascule==reference, "Fixed", "De novo"))
  }) %>% setNames(names(assigned_missing_cosmic))
  
  assigned_degasperi = lapply(names(assigned_missing_degasperi), function(tid) {
    assigned = assigned_missing_degasperi[[tid]]$assigned_tp
    data.frame(list(bascule=assigned, reference=names(assigned)), row.names=NULL) %>% 
      dplyr::mutate(type=tid, catalogue="Degasperi et al.") %>% 
      dplyr::mutate(signature_type=ifelse(bascule%in%unlist(get_denovo_signames(x_orig)), 
                                          "De novo", "Fixed")) %>% 
      dplyr::filter(signature_type=="De novo")
  }) %>% setNames(names(assigned_missing_degasperi))
  
  assigned <<- lapply(get_types(x_orig), function(tid) {
    rbind(assigned_cosmic[[tid]], assigned_degasperi[[tid]])
  }) %>% setNames(get_types(x_orig))
  
  
  whole_catalogue = list(
    SBS=rbind(reference_cat_cosmic$SBS[assigned_cosmic$SBS$reference,],
              reference_cat_degasperi$SBS[assigned_degasperi$SBS$reference,]),
    DBS=rbind(reference_cat_cosmic$DBS[assigned_cosmic$DBS$reference,],
              reference_cat_degasperi$DBS[assigned_degasperi$DBS$reference,])
  )
  
  similarity_df_all = lapply(get_types(x_orig), function(tid) {
    lsa::cosine(rbind(whole_catalogue[[tid]],
                      get_denovo_signatures(x_orig, matrix=T)[[tid]]) %>% 
                  t() %>% as.matrix()) %>% as.data.frame() %>% 
      tibble::rownames_to_column(var="bascule") %>% 
      reshape2::melt(variable.name="reference", value.name="cosine") %>% 
      dplyr::filter(bascule%in%assigned[[tid]]$bascule,
                    reference%in%assigned[[tid]]$reference) %>% 
      dplyr::right_join(assigned[[tid]])
  }) %>% dplyr::bind_rows() %>% tibble::as_tibble() %>% 
    dplyr::mutate(bascule=fct_reorder(bascule, cosine, max)) %>% 
    dplyr::mutate(signature_type=factor(signature_type, levels=c("Fixed","De novo")))
  
  similarity_df = similarity_df_all %>% dplyr::filter(type=="SBS")
  similarity_df %>% 
    dplyr::mutate(b_number=1:nrow(similarity_df)) %>% 
    ggplot() +
    geom_bar(aes(x=cosine, y=bascule, 
                 fill=catalogue, alpha=ifelse(cosine>=0.8,1,0.3)), 
             stat="identity", width=0.5) +
    geom_vline(xintercept=0.8, linewidth=0.2, linetype="solid") +
    facet_grid(signature_type ~ ., scales="free", space="free") +
    scale_fill_manual(values=cls_catalogues, name="Catalogue") +
    scale_alpha_continuous(limits=c(0,1)) +
    ylab("BASCULE") + xlab("Cosine similarity") +
    guides(y.sec=ggh4x::guide_axis_manual(breaks=similarity_df$bascule, 
                                          labels=similarity_df$reference, 
                                          title="Known"),
           alpha="none") +
    theme_bw() + 
    theme(panel.grid=element_blank(), 
          axis.ticks.y=element_blank(),
          legend.position="bottom")
  
}


# input : mapped function output
# output: ggplot object - heatmap of similarity between basilica, serena and refrence signatures
plot_heatmp_comp = function(mapped) {
  
  fixed_common = mapped$fixed_common %>% dplyr::mutate(rowType="Fixed", colType="Common (serena)")
  fixed_rare = mapped$fixed_rare %>% dplyr::mutate(rowType="Fixed", colType="Rare (serena)")
  fixed_reference = mapped$fixed_reference %>% dplyr::mutate(rowType="Fixed", colType="COSMIC")
  
  denovo_common = mapped$denovo_common %>% dplyr::mutate(rowType="Denovo", colType="Common (serena)")
  denovo_rare = mapped$denovo_rare %>% dplyr::mutate(rowType="Denovo", colType="Rare (serena)")
  denovo_reference = mapped$denovo_reference %>% dplyr::mutate(rowType="Denovo", colType="COSMIC")
  
  all = rbind(
    fixed_common, 
    fixed_rare, 
    fixed_reference, 
    denovo_common, 
    denovo_rare, 
    denovo_reference
  )
  
  all$colType = factor(all$colType, levels=c("Common (serena)", "Rare (serena)", "COSMIC"))
  
  p = ggplot(data=all, aes(x=first, y=second, fill=values)) + 
    geom_tile(color="black") + 
    facet_grid(colType ~ rowType, scales="free", space="free") + 
    geom_text(aes(label=round(values, digits=2))) + 
    scale_fill_gradient(low="gray", high="dodgerblue1") + 
    labs(
      title="Signatures Mapping", 
      x="Basilica Signatures", 
      y="[COSMIC + Serena] Signatures", 
      subtitle="Basilica vs. [COSMIC + Serena]", 
      caption="COSMIC v3.4 - October 2023", 
      fill="Cosine Similarity") + 
    theme_bw() + 
    theme(
      legend.title=element_text(angle=90, vjust=1), 
      axis.text.x=element_text(angle=45, vjust=1, hjust=1)
      )
  
  return(p)
}



plot_scores_nmf = function(x, types=get_types(x), remove_outliers=FALSE) {
  
  scores = get_scores(x) %>% 
    dplyr::select(-value) %>% unique() %>% 
    
    dplyr::filter(score_id=="bic", type%in%types, parname=="K") %>%
    dplyr::select(-score_id, -parname) %>% 
    
    dplyr::group_by(type, value_fit) %>%
    dplyr::summarise(score=min(score)) %>%
    
    dplyr::group_by(type) %>%
    dplyr::mutate(is.min=score==min(score),
                  label=replace(is.min, is.min==T, "Best fit"), 
                  label=replace(label, label==F, NA))

  if (remove_outliers)
    scores=scores %>%
    dplyr::group_by(score_id, parname, type) %>%
    dplyr::mutate(is.outlier=score %in% boxplot.stats(score)$out) %>%
    dplyr::filter(!is.outlier)
  
  scores %>%
    ggplot(aes(x=as.integer(value_fit), y=score)) +
    geom_point(aes(shape=label), size=2) +
    geom_line(linetype="solid") +
    ggrepel::geom_label_repel(aes(label=label), box.padding=0.05, size=2,
                             na.rm=T, show.legend=F, 
                             direction="y", point.size=10) +
    ggh4x::facet_nested_wrap(type ~ ., scales="free", nrow=length(types)) +
    theme_bw() + labs(title="Model selection scores") +
    xlab("No. of Denovo Signatures") + ylab("BIC") +
    scale_x_continuous(breaks=scores %>% dplyr::pull(value_fit) %>% unique()) + 
    theme(legend.title=element_blank(), legend.position="none")
}



plot_centroids_etiology = function(
    x, 
    sbs_aetiology_path, 
    dbs_aetiology_path, 
    exposure_thr=0.05, 
    quantile_thr=0.9, 
    custom_theme=NULL
    #manual=F, 
    #clusters=get_cluster_labels(x), 
    #sig_cls=NULL, 
    #types=get_types(x)
  ) {
  
  aetiology = rbind(
    read.table(file=sbs_aetiology_path, sep=',', header=TRUE), 
    read.table(file=dbs_aetiology_path, sep=',', header=TRUE)
  )
  
    
  scores = get_clusters_score(x, get_types(x), exposure_thr, quantile_thr) %>% subset(significance == T)
  
  # add denovo signatures to aetiology dataframe as UNKNOWN aetiology
  #(scores %>% dplyr::pull(signature) %>% unique) %in% (aetiology %>% dplyr::pull(signature))
  dn_sigs = setdiff(
    (scores %>% dplyr::pull(signature) %>% unique), 
    (aetiology %>% dplyr::pull(signature))
  )
  aetiology = rbind( aetiology, data.frame(signature=dn_sigs, aetiology=rep("UNKNOWN", length(dn_sigs))) )
  
  # quality control
  if ( !(all( (scores %>% dplyr::pull(signature) %>% unique) %in% (aetiology %>% dplyr::pull(signature)) ) == TRUE) ) {
    warning("not all signatures in aetiology")
  }
  
  # need tidyverse package to laod
  df = list(aetiology, scores) %>% reduce(full_join, by='signature') %>% select(signature, aetiology, cluster) %>% na.omit()
  
  # Create a dataframe with counts of cluster and aetiology combinations
  heatmap_data = df %>%
    group_by(cluster, aetiology) %>%
    summarise(count=n(), .groups='drop') %>%
    mutate(present=ifelse(count > 0, 1, 0))
  
  #a = plot_centroids(x_mapped, get_types(x_mapped), cls=NULL, sort_by=NULL, exposure_thr=0, quantile_thr=0.9) + 
  #  theme(
  #    text=element_text(size=20), 
  #    legend.title=element_text(size=18)
  #  )
  
  a = plot_centroids(
    x, types=get_types(x), cls=NULL, sort_by=NULL, exposure_thr=exposure_thr, quantile_thr=quantile_thr
  ) + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), 
            plot.margin=unit(c(5.5, 5.5, 0, 5.5), "pt")) + 
    labs(y="Density") + custom_theme
  
  
  # Plot the heatmap
  b = ggplot(heatmap_data, aes(x=cluster, y=aetiology)) +
    geom_tile(aes(fill=factor(present)), color="white", width=0.9, height=0.9) +
    scale_fill_manual(values=c("0"="white", "1"="gray"), guide="none") +
    theme_bw() +
    labs(x="Clusters", y="Aetiology") + # title="Cluster vs. Aetiology"
    theme(axis.text.x=element_text(angle=45, hjust=1), legend.position="none", plot.margin=unit(c(5.5, 5.5, 0, 5.5), "pt")) + custom_theme
  
  design = "aaa
             aaa
             bbb"
  
  p = patchwork::wrap_plots(a, b, ncol=1, design=design)
  
  return(p)
}



# plot centroids by giving the significant signatures in each cluster
custom_centroid_plot = function(x, sig_cls, col_palette,
                                sbs_aetiology_path, 
                                dbs_aetiology_path) {

  # CENTROIDS
  cluster_levels = gtools::mixedsort(get_cluster_labels(x))
  centroids = plot_centroids(x, signatures_list=unlist(sig_cls), 
                             cls=col_palette, sample_levels=cluster_levels) +
    xlab("Clusters") + ylab("Relative exposures") +
    theme(plot.margin=unit(c(t=5.5,r=5.5,l=5.5,b=0),"pt"), 
          panel.grid=element_blank()) +
    scale_y_continuous(breaks=c(0,1))
  
  # ETIOLOGY
  aetiology = rbind(read.table(file=sbs_aetiology_path, sep=',', header=TRUE), 
                    read.table(file=dbs_aetiology_path, sep=',', header=TRUE)) %>% 
    dplyr::mutate(aetiology=replace(aetiology,aetiology=="UNKNOWN","Unknown")) %>% 
    dplyr::rename(sigs=signature)
  
  # check if all signatures are present in etiology signatures list
  aetiology = aetiology %>% 
    dplyr::bind_rows(list(sigs=unlist(sig_cls) %>% setdiff(aetiology$sigs), 
                          aetiology="Unknown") %>% expand.grid()) %>% 
    dplyr::filter(sigs%in%unlist(sig_cls)) %>% 
    dplyr::left_join(
      lapply(names(sig_cls), function(cl_id) { 
        data.frame(clusters=cl_id, sigs=sig_cls[[cl_id]]) 
        }) %>% dplyr::bind_rows()
    )

  # Create a dataframe with counts of cluster and aetiology combinations
  heatmap_data = aetiology %>%
    group_by(clusters, aetiology) %>%
    summarise(count=n()) %>% dplyr::ungroup() %>% 
    mutate(is.present=ifelse(count > 0, 1, 0))
  

  # Plot the heatmap
  row_levels = c("Unknown", sort(heatmap_data$aetiology, decreasing=T) %>% purrr::discard_at("Unknown")) %>% unique()
  heatmap = ggplot(heatmap_data, 
                   aes(x=factor(clusters, levels=cluster_levels), 
                       y=factor(aetiology, levels=row_levels))) +
    geom_tile(aes(fill=factor(is.present)), 
              color="grey5", width=0.9, height=0.5) +
    scale_fill_manual(values=c("0"="white", "1"="gray50"), guide="none") +
    scale_y_discrete(position="right") +
    theme_bw() +
    labs(x="Clusters", y="Aetiology") + ylab("") +
    theme(plot.margin=unit(c(t=0,r=5.5,l=5.5,b=5.5),"pt")) +
    theme(axis.ticks.y=element_blank(), 
          panel.spacing.y=unit(100, "pt"),
          panel.grid.major=element_blank())
  
  return(list("centroids"=centroids, "heatmap"=heatmap))
}


renormalize_sigs = function(catalogue, thr) {
  catalogue %>% 
    dplyr::mutate(value=replace(value, value<thr, 0)) %>%
    dplyr::group_by(sigs) %>%
    dplyr::mutate(value=value/(sum(value)+1e-10)) %>%
    dplyr::ungroup()
}

plot_mirrored_sigs = function(denovo_cat, ref_cat, 
                              dn_name, ref_name,
                              cat_name, type, thr=0.02) {
  breaks = c("De novo", "COSMIC v3.4", "Degasperi et al.")
  title_name = strsplit(cat_name, split=" ")[[1]][1]
  
  input_df = denovo_cat %>% 
    dplyr::filter(sigs==dn_name) %>% 
    reformat_contexts(what=type) %>% 
    dplyr::mutate(Catalogue=breaks[1]) %>% 
    
    renormalize_sigs(thr=thr) %>%  
    
    dplyr::bind_rows(
      ref_cat %>% dplyr::filter(sigs %in% c(ref_name)) %>% 
        reformat_contexts(what=type) %>% dplyr::select(-type) %>% 
        dplyr::mutate(value=-value) %>% 
        dplyr::mutate(Catalogue=cat_name)
    )
  
  max_abs = max(max(input_df$value), max(abs(input_df$value)))

  input_df %>% {
    ggplot(.) +
      geom_bar(aes(x=context, y=value, fill=Catalogue), show.legend=T, 
               stat="identity", position="identity") +
      facet_grid(~variant) +
      scale_fill_manual(values=cls_catalogues, 
                        breaks=breaks,
                        limits=breaks) +
      scale_y_continuous(breaks=pretty(.[["value"]]),
                         labels=abs(pretty(.[["value"]]))) +
      theme_bw() +
      theme(axis.text=element_blank(),
            axis.ticks=element_blank(),
            panel.grid=element_blank()) +
      ylab("Density") + xlab("Trinucleotide contexts") +
      labs(title=paste0("Signatures ", dn_name, " and ", ref_name)) +
      ylim(-max_abs, max_abs)
    }
}


plot_custom_exposures=function(
    x,
    types=get_types(x),
    samples=get_samples(x),
    clusters=get_cluster_labels(x),
    sample_name=FALSE,
    color_palette=NULL,
    sort_by=NULL, 
    signature_list
) {
  
  exposures=lapply(types, function(t)
    get_exposure(x, types=types, samples=samples,
                 clusters=clusters, add_groups=TRUE)[[t]] %>%
      dplyr::mutate(type=t)) %>%
    do.call(rbind, .)
  
  if (is.null(color_palette)) cls=gen_palette(x, types=sort(types))
  
  # merging signatures where their exposure in all the samples are below the threshold
  to_keep=exposures$sigs %>% unique(); sigs_levels=NULL
  
  exp = exposures %>% dplyr::mutate(sigs=ifelse( (sigs %in% signature_list), sigs, "Others"))
  
  cls["Others"] = "gainsboro"
  
  
  p = plot_exposures_aux(exposures=exp, cls=cls,
                          titlee="Exposure",
                          sample_name=FALSE,
                          sort_by=NULL,
                          sigs_levels=sigs_levels) +
    scale_fill_manual(values=cls, breaks=to_keep)
  
  return(p)
}



get_color_palette = function(cosmic, degasperi, sig_cls_organ_all) {
  
  ref_names = unique(unlist(lapply(sig_cls_organ_all, unlist))) %>% 
    purrr::keep(function(i) i %in% c(cosmic$sigs, degasperi$sigs))
  set.seed(123)
  cls_ref = unique(c(
    yarrr::piratepal(palette="info2", mix.col="yellow", mix.p=0) %>% purrr::discard_at("pink"),
    yarrr::piratepal(palette="appletv", mix.col="yellow", mix.p=0),
    yarrr::piratepal(palette="nemo", mix.col="yellow", mix.p=0),
    yarrr::piratepal(palette="espresso", mix.col="yellow", mix.p=0))) %>% 
    sample(size=length(ref_names))
  
  seeds = c(33,2211,432) %>% setNames(names(sig_cls_organ_all))
  cls_dn = cls = list()
  for(organ_id in names(sig_cls_organ_all)) {
    sig_cls_organ = sig_cls_organ_all[[tolower(organ_id)]]
    
    other_names = cls_dn %>% lapply(names) %>% unlist() %>% setNames(NULL)
    dn_names = setdiff(unlist(sig_cls_organ), ref_names)
    n_cls = dn_names %>% length()
    
    set.seed(seeds[[organ_id]])
    cls_dn[[organ_id]] = yarrr::piratepal(palette="basel", mix.col="yellow", mix.p=0) %>% 
      purrr::discard_at(c("pink", other_names, names(ref_names))) %>% 
      sample(size=n_cls)
    
    cls[[organ_id]] = c(cls_ref %>% setNames(ref_names), 
            cls_dn[[organ_id]] %>% setNames(dn_names))
    
    cls[[organ_id]] = cls[[organ_id]][!is.na(names(cls[[organ_id]]))]
  }
  
  return(cls)
}





