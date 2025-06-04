theme_legend = theme(legend.text=element_text(size=7.4),
                     legend.title=element_text(size=10.6),
                     legend.key.size=unit(0.2,"cm"),
                     legend.key.height=unit(0.2,"cm"),
                     legend.key.width=unit(0.2,"cm"))

theme_text = theme(axis.title=element_text(size=7.4),
                   axis.text=element_text(size=5.35), 
                   plot.title=element_text(size=10.6),
                   plot.subtitle=element_text(size=7.4),
                   strip.text=element_text(size=5.35), 
                   plot.caption=element_text(size=5.35))


# list of relevant signatures in each tumor type and cluster
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


get_degasperi = function(degasperi_sbs="~/Google Drive/My Drive/work/bascule_shared/zenodo/real_data_validation/data/degasperi_sbs.Rds",
                         degasperi_dbs="~/Google Drive/My Drive/work/bascule_shared/zenodo/real_data_validation/data/degasperi_dbs.Rds") {
  readRDS(degasperi_sbs) %>% 
    bascule:::wide_to_long(what="beta") %>% dplyr::mutate(type="SBS") %>% 
    dplyr::bind_rows(
      readRDS(degasperi_dbs) %>% 
        bascule:::wide_to_long(what="beta") %>% dplyr::mutate(type="DBS")
    )
}


get_cosmic = function() {
  COSMIC_sbs %>% 
    bascule:::wide_to_long(what="beta") %>% dplyr::mutate(type="SBS") %>% 
    dplyr::bind_rows(
      COSMIC_dbs %>% 
        bascule:::wide_to_long(what="beta") %>% dplyr::mutate(type="DBS")
    )
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
