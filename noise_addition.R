library(pso)
library(dplyr)
library(ggplot2)
devtools::load_all('~/Workspace/breaktools/')
load("~/Workspace/Everything/DSB_Paper/data/data_matices.rda")
load("~/Workspace/Everything/DSB_Paper/data/training_parameters.rda")


mtx_source <- do.call(cbind, data_matrix)
mtx_target <- do.call(cbind, repliseq1_data_matrix)


both_ldf = dplyr::bind_rows(
  reshape2::melt(unname(mtx_source)) %>% dplyr::mutate(dataset="source"),
  reshape2::melt(unname(mtx_target)) %>% dplyr::mutate(dataset="target")) %>%
  dplyr::rename(fraction="Var1", bin="Var2") %>%
  dplyr::group_by(dataset, bin) %>%
  dplyr::filter(all(!is.na(value))) %>%
  dplyr::mutate(fraction_highest=fraction[which.max(value)]) %>%
  dplyr::ungroup()

noise_df = both_ldf %>%
  dplyr::group_by(dataset, fraction_highest, fraction) %>%
  dplyr::summarise(vsd=sd(value), value=mean(value)) %>%
  dplyr::group_by(fraction_highest, fraction) %>%
  dplyr::summarise(
    noise_mean=value[dataset=="target"]-value[dataset=="source"],
    noise_sd=pmax(0.01, vsd[dataset=="target"]-vsd[dataset=="source"]))

noise_add = function(mtx, noise_df) {
  mtx_colnames = colnames(mtx)
  colnames(mtx) = as.character(1:ncol(mtx))
  mtx_colmap = mtx_colnames 
  names(mtx_colmap) = colnames(mtx) 
  mtx_noisy = reshape2::melt(unname(mtx)) %>%
    dplyr::rename(fraction="Var1", bin="Var2") %>%
    dplyr::group_by(bin) %>%
    dplyr::mutate(
      fraction_highest=ifelse(all(is.na(value)), NA_integer_, fraction[which.max(value)])) %>%
    dplyr::ungroup() %>%
    dplyr::left_join(noise_df, by=c("fraction_highest", "fraction")) %>%
    dplyr::mutate(value_noise=rnorm(value, mean=noise_mean, sd=noise_sd), value_fitted=value+value_noise) %>%
    #dplyr::mutate(value_fitted=pmin(pmax(value_fitted, 0), 1)) %>%
    dplyr::group_by(bin) %>%
    dplyr::mutate(value_fitted=value_fitted/max(value_fitted)) %>%
    dplyr::ungroup() %>%
    reshape2::dcast(fraction~bin, value.var="value_fitted") %>%
    tibble::column_to_rownames("fraction") %>%
    as.matrix()
  colnames_sorted = colnames(mtx)[sort(as.integer(colnames(mtx)))]
  rownames(mtx_noisy) = c(8,7,6,5,4,3,2,1)
  mtx_new = mtx_noisy[,colnames_sorted] 
  colnames(mtx_new) = mtx_colmap[colnames(mtx_new)]
  mtx_new
}


data_noise = lapply(data_matrix, noise_add, noise_df)


i= 202
dplyr::bind_rows(
  reshape2::melt(data_matrix[[i]]) %>% dplyr::mutate(set="Zhao"),
  reshape2::melt(data_noise[[i]]) %>% dplyr::mutate(set="Zhao + Noise"),
  reshape2::melt(repliseq1_data_matrix[[i]]) %>% dplyr::mutate(set="Wei Run 01")) %>%
  dplyr::mutate(Var1=factor(Var1), training_id=gsub("_", ":", names(data_matrix)[i])) %>%
  ggplot() +
  geom_tile(aes(x=Var2, y=Var1, fill=value)) +
  facet_grid(set~training_id, scales="free", space="free") +
  scico::scale_fill_scico(palette = "bilbao") +
  theme_bw() +
  theme(axis.title.y=element_blank(), axis.ticks.y=element_blank(), panel.grid.major=element_blank(), panel.grid.minor=element_blank())







