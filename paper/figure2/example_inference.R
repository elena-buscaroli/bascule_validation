library(magrittr)
library(ggplot2)
devtools::load_all("~/GitHub/basilica/")

fit_simul = readRDS("~/Dropbox/dropbox_shared/2022. Basilica/simulations/fits/fits_dn.matched.2011/simul_fit.N500.G3.s11.matched.2011.Rds")
# fit_simul = readRDS("~/Dropbox/dropbox_shared/2022. Basilica/simulations/fits/fits_dn.matched.2011/simul_fit.N150.G3.s22.matched.2011.Rds")


bas_mapped = fit_simul$fit.0.auto %>% 
  convert_dn_names(reference_cat=get_signatures(fit_simul$dataset, matrix=T), cutoff=0.7) %>% 
  merge_clusters()

clusters_new = c("G3","G2","G1") %>% setNames(c("G0","G1","G3"))
get_exposure(bas_mapped, add_groups=T)[["SBS"]] %>% 
  dplyr::rowwise() %>% 
  dplyr::mutate(clusters=clusters_new[clusters], method="Predicted") %>% 
  dplyr::ungroup() %>% 
  
  dplyr::bind_rows(
    get_exposure(fit_simul$dataset, add_groups=T)[["SBS"]] %>% 
      dplyr::mutate(method="Simulated") 
  ) %>% 
  
  ggplot() +
  geom_bar(aes(x=samples, y=value, fill=sigs), stat="identity") +
  facet_grid(method ~ clusters, scales="free_x", space="free_x") +
  scale_fill_manual(values=gen_palette(n=5)) +
  
  theme_bw() +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), 
        panel.grid.major.x=element_blank(), panel.grid.major.y=element_blank(),
        legend.position="bottom") +
  labs(fill="Signatures") + ylab("") + xlab("Samples")



plot_exposures(bas_mapped, types="SBS")
plot_exposures(fit_simul$dataset, types="SBS")

