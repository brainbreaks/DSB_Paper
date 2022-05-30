library(rtracklayer)
library(readr)
library(dplyr)
library(ggplot2)
library(cowplot)
library(stringr)
library(dplyr)
library(scico)
devtools::load_all('~/Workspace/breaktools/')
setwd("~/Workspace/DSB_paper")

#
# Load data
#
repliseq1_df = readr::read_tsv("~/Workspace/Datasets/B400_RS_001/results/DMSO_repliseq50000.tsv")
repliseq2_df = readr::read_tsv("~/Workspace/Datasets/B400_RS_002/results/DMSO_repliseq50000.tsv")

# x = repliseq1_df %>%
#   dplyr::inner_join(repliseq2_df, by=c("repliseq_chrom", "repliseq_start", "repliseq_end", "repliseq_fraction"))

repliseq_df = readr::read_tsv("~/Workspace/Datasets/zhao_bmc_repliseq_2020/results/zhao_mESC_repliseq50000.tsv") %>%
  dplyr::mutate(repliseq_fraction=ceiling(repliseq_fraction/2)) %>%
  dplyr::group_by(repliseq_chrom, repliseq_start, repliseq_end, repliseq_fraction) %>%
  dplyr::summarise(repliseq_value=sum(repliseq_value, na.rm=T), repliseq_value_abs=sum(repliseq_value_abs)) %>%
  dplyr::group_by(repliseq_chrom, repliseq_start, repliseq_end) %>%
  dplyr::mutate(repliseq_value=repliseq_value/max(repliseq_value)) %>%
  dplyr::ungroup()
rfd_ranges = rtracklayer::import.bw("~/Workspace/Datasets/petryk_okseq_mESC/GSM3290342_OK-seq_smooth_results_w1000_s30_d30_z1_filter.bw") %>%
  as.data.frame() %>%
  dplyr::rename(rfd_chrom="seqnames", rfd_start="start", rfd_end="end", rfd_score="score") %>%
  dplyr::select(rfd_chrom, rfd_start, rfd_end, rfd_score) %>%
  df2ranges(rfd_chrom, rfd_start, rfd_end)
iz_ranges = readr::read_tsv("~/Workspace/Datasets/petryk_okseq_mESC/GSM3290342_Ok_IZ.txt") %>%
  tidyr::extract(col=break_ID, into=c("iz_chrom", "iz_start", "iz_end"), regex="([^:]+):(\\d+)-(\\d+)") %>%
  dplyr::mutate(iz_start=as.numeric(iz_start), iz_end=as.numeric(iz_end)) %>%
  dplyr::select(iz_chrom, iz_start, iz_end) %>%
  df2ranges(iz_chrom, iz_start, iz_end)
colfork_ranges = as.data.frame(rtracklayer::import.bed("data/zhao_colliding_forks.bed")) %>%
  dplyr::rename(colfork_chrom="seqnames", colfork_start="start", colfork_end="end", colfork_width="width") %>%
  dplyr::select(-strand)%>%
  dplyr::mutate(colfork_id=paste0(colfork_chrom, ":", colfork_start, "-", colfork_end)) %>%
  df2ranges(colfork_chrom, colfork_start, colfork_end)

#
# Find termination zones
#
tz_df = iz_ranges %>%
  innerJoinByOverlaps(colfork_ranges) %>%
  dplyr::group_by(colfork_chrom, colfork_start, colfork_end) %>%
  dplyr::filter(dplyr::n()==2) %>%
  dplyr::summarize(colfork_min_pos=min(iz_end), colfork_max_pos=max(iz_start)) %>%
  dplyr::ungroup() %>%
  df2ranges(colfork_chrom, colfork_min_pos, colfork_max_pos) %>%
  innerJoinByOverlaps(rfd_ranges) %>%
  dplyr::arrange(colfork_chrom, colfork_start, colfork_end, rfd_chrom, rfd_start) %>%
  dplyr::group_by(colfork_chrom, colfork_start, colfork_end) %>%
  dplyr::do((function(z) {
    zz <<- z
    rfd_model = loess(rfd_score ~ rfd_start, data=z, span=0.3)
    rfd_pred = data.frame(rfd_start=seq(z$rfd_start[which.max(z$rfd_score)], z$rfd_start[which.min(z$rfd_score)], 1e3)) %>%
      dplyr::mutate(rfd_chrom=z$colfork_chrom[1], rfd_end=rfd_start+1e3)
    rfd_pred$rfd_score = predict(rfd_model, rfd_pred)
    rfd_pred_ranges = rfd_pred %>%
      dplyr::filter(abs(rfd_score) < 0.01) %>%
      df2ranges(rfd_chrom, rfd_start, rfd_end) %>%
      GenomicRanges::reduce(min.gapwidth=1e4)

    # ggplot(z, aes(x=rfd_start, y=rfd_score)) +
    #   geom_line() +
    #   geom_smooth(method="loess", span=.3) +
    #   geom_vline(xintercept=start(rfd_pred_ranges)[1])

    data.frame(tz_chrom=z$colfork_chrom[1], tz_start=start(rfd_pred_ranges)[1], tz_end=end(rfd_pred_ranges)[1], tz_alternative_n=length(rfd_pred_ranges))
  })(.)) %>%
  dplyr::ungroup() %>%
  dplyr::select(dplyr::starts_with("tz_"))

#
# Split into training regions
#
training_step = 1e5
training_margin = 1e5
training_width = 1e6
training_binsize = 5e4
training_regions_df = as.data.frame(GenomicRanges::reduce(colfork_ranges)) %>%
  dplyr::select(-strand) %>%
  dplyr::rename(colfork_reduced_chrom="seqnames", colfork_reduced_start="start", colfork_reduced_end="end", colfork_reduced_width="width") %>%
  dplyr::mutate(seq_start=colfork_reduced_start-training_margin, seq_end=colfork_reduced_end-training_width+training_margin) %>%
  dplyr::mutate(seq_start=round(seq_start/training_binsize)*training_binsize, seq_end=round(seq_end/training_binsize)*training_binsize) %>%
  dplyr::filter(seq_end>=seq_start) %>%
  dplyr::group_by(colfork_reduced_chrom, colfork_reduced_start, colfork_reduced_end) %>%
  dplyr::summarize(training_chrom=colfork_reduced_chrom[1], training_start=seq(seq_start, seq_end, by=training_step), training_end=training_start+training_width) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(training_id=paste0(training_chrom, "_", format(training_start, scientific=F, justify="none", trim=T), "-", format(training_end, scientific=F, justify="none", trim=T)))

#
# Create training and label data.frames
#
data_df = training_regions_df %>%
  df2ranges(training_chrom, training_start, training_end) %>%
  innerJoinByOverlaps(repliseq_df %>% df2ranges(repliseq_chrom, repliseq_start, repliseq_end), minoverlap=training_binsize) %>%
  dplyr::group_by(training_chrom, training_start, training_end) %>%
  dplyr::mutate(repliseq_rel_start=repliseq_start-training_start) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(training_id, repliseq_chrom, repliseq_start, repliseq_end) %>%
  dplyr::mutate(repliseq_value=repliseq_value/max(repliseq_value)) %>%
  dplyr::ungroup()

label_df = training_regions_df %>%
  dplyr::group_by(training_id, training_chrom, training_start, training_end) %>%
  dplyr::summarize(training_tz_start=seq(training_start, training_end-1, training_binsize), training_tz_end=training_tz_start+training_step) %>%
  dplyr::mutate(training_tz_rel_start=training_tz_start-training_start) %>%
  df2ranges(training_chrom, training_tz_start, training_tz_end) %>%
  leftJoinByOverlaps(tz_df %>% df2ranges(tz_chrom, tz_start, tz_end)) %>%
  dplyr::group_by(training_id, training_chrom, training_start, training_end, tz_start, tz_end, training_tz_rel_start, training_tz_start, training_tz_end) %>%
  dplyr::summarize(training_tz_present=as.integer(!is.na(tz_alternative_n))) %>%
  dplyr::arrange(training_tz_rel_start) %>%
  dplyr::inner_join(data_df %>% dplyr::group_by(repliseq_chrom, repliseq_rel_start) %>% dplyr::summarize(repliseq_value.mean=weighted.mean(repliseq_fraction, repliseq_value)), by=c("training_chrom"="repliseq_chrom", "training_tz_rel_start"="repliseq_rel_start")) %>%
  dplyr::group_by(training_id, training_start, training_end, tz_start, tz_end) %>%
  dplyr::mutate(training_tz_present=ifelse(repliseq_value.mean==max(repliseq_value.mean) & training_width-training_tz_rel_start>=training_binsize*2, training_tz_present, 0)) %>%
  dplyr::ungroup()


repliseq_compare_df = dplyr::bind_rows(
  repliseq_df %>% dplyr::mutate(repliseq_value_abs=repliseq_value_abs/2, Run="Zhao"),
  repliseq1_df %>% dplyr::mutate(Run="RS1"),
  repliseq2_df %>% dplyr::mutate(Run="RS2"))
ggplot(repliseq_compare_df %>% dplyr::filter(repliseq_value_abs<=5000)) +
  ggridges::geom_density_ridges(aes(y=Run, x=repliseq_value_abs)) +
  labs(x="Number of reads per bin/fraction", y="")
ggplot(repliseq_compare_df %>% dplyr::filter(repliseq_value<=1)) +
  ggridges::geom_density_ridges(aes(y=Run, x=repliseq_value)) +
  labs(x="Relative fraction abundance bin", y="")


#
# Plot training samples distribution
#
training_regions_ggplot = label_df %>%
  dplyr::group_by(training_id) %>%
  dplyr::summarize(training_tz_count=as.character(sum(training_tz_present, na.rm=T))) %>%
  dplyr::group_by(training_tz_count) %>%
  dplyr::summarize(training_samples_count=dplyr::n()) %>%
  dplyr::bind_rows(data.frame(training_tz_count="Total", training_samples_count=label_df %>% dplyr::distinct(training_id) %>% nrow()))
ggplot(training_regions_ggplot) +
  geom_bar(aes(x=training_tz_count, y=training_samples_count), stat="identity") +
  labs(y="Training samples", x="No. of termination zones in sample") +
  scale_y_continuous(breaks = scales::pretty_breaks(n=20))+
  ggpubr::theme_pubclean()

#
# Plot some examples
#
dplyr::bind_rows(
  data_df %>% dplyr::mutate(set="Data", x=repliseq_rel_start, y=repliseq_fraction, value=repliseq_value),
  label_df %>% tidyr::crossing(data.frame(y=1:2)) %>% dplyr::mutate(set="Labels", x=training_tz_rel_start, value=training_tz_present)) %>%
  dplyr::filter(training_id %in% training_regions_df$training_id[c(1,500,1000,1500)]) %>%
  ggplot() +
  geom_tile(aes(x=as.numeric(x), y=y, fill=value)) +
  facet_grid(set~training_id, scales="free", space="free") +
  scico::scale_fill_scico(palette = "bilbao") +
  scale_x_continuous(labels=scales::label_number(accuracy=1, scale=1e-3, suffix="K"))  +
  theme_bw() +
  theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.title.x=element_blank(), axis.ticks.y=element_blank(), panel.grid.major=element_blank(), panel.grid.minor=element_blank())


#
# Final matrix
#
data_matrix = data_df %>%
  dplyr::mutate(repliseq_rel_start=format(repliseq_rel_start, scientific=F, justify="none", trim=T)) %>%
  reshape2::dcast(training_id+repliseq_fraction ~ repliseq_rel_start, value.var="repliseq_value") %>%
  dplyr::arrange(dplyr::desc(repliseq_fraction)) %>%
  split(.$training_id) %>%
  lapply(FUN=function(z) z %>% `rownames<-`(NULL) %>% tibble::column_to_rownames("repliseq_fraction") %>% select(-training_id) %>% as.matrix())
label_matrix = label_df %>%
  dplyr::mutate(training_tz_rel_start=format(training_tz_rel_start, scientific=F, justify="none", trim=T), training_tz_rel_start=factor(training_tz_rel_start, unique(training_tz_rel_start))) %>%
  reshape2::dcast(training_id ~ training_tz_rel_start, value.var="training_tz_present", drop=F) %>%
  split(.$training_id) %>%
  lapply(FUN=function(z) z %>% select(-training_id) %>% as.matrix()) %>%
  `[`(names(training_matrix))

#
# PLoting final data for testing
#
i = 500
dplyr::bind_rows(
  reshape2::melt(training_matrix[[i]]) %>% dplyr::mutate(set="Data"),
  reshape2::melt(label_matrix[[i]]) %>% dplyr::mutate(set="Labels")) %>%
  dplyr::mutate(Var1=factor(Var1), training_id=gsub("_", ":", names(training_matrix)[i])) %>%
  ggplot() +
  geom_tile(aes(x=Var2, y=Var1, fill=value)) +
  facet_grid(set~training_id, scales="free", space="free") +
  scico::scale_fill_scico(palette = "bilbao") +
  theme_bw() +
  theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank(), panel.grid.major=element_blank(), panel.grid.minor=element_blank())

