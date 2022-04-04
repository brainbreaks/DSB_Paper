# setwd("/mnt/e/Workspace/Emily")
# library(randomForest)
library(dplyr)
library(reshape2)
library(readr)
library(GenomicRanges)
library(GenomicFeatures)
library(ggplot2)
library(tidyr)
library(baseline)
library(smoother)
library(dbscan)
library(ggpmisc)
library(rtracklayer)
library(ggbeeswarm)
library(Rtsne)
library(pracma)
library(ggiraph)
library(ggannotate)
devtools::load_all('~/Workspace/breaktools/')
setwd("~/Workspace/Everything")

correct_annotation_chrom = function() {
  rdc_df = readr::read_tsv("~/Workspace/Datasets/HTGTS/rdc_pnas_mm10.tsv")
  replication_df = readr::read_tsv("data/replication_subsets.tsv") %>%
    dplyr::rename(replication_chrom.old="replication_chrom") %>%
    dplyr::inner_join(rdc_df %>% dplyr::select(rdc_cluster, replication_chrom=rdc_chrom), by="rdc_cluster") %>%
    dplyr::select(rdc_cluster, replication_chrom, replication_start, replication_end, replication_strand)

  readr::write_tsv(replication_df, file="data/replication_subsets.tsv")
}

annotate_rdc_repliseq_drop_duplicates = function() {
  rdc_df = dplyr::bind_rows(readr::read_tsv("~/Workspace/Datasets/HTGTS/rdc_pnas_mm10.tsv"), readr::read_tsv("~/Workspace/Datasets/HTGTS/rdc_macs2_mm10.tsv"))
  replication_df = readr::read_tsv("data/replication_subsets.tsv") %>%
    dplyr::mutate(annotation_rdc_cluster=rdc_cluster) %>%
    dplyr::select(-rdc_cluster)
  replication_ranges = replication_df %>%
    dplyr::mutate(seqnames=replication_chrom, start=pmin(replication_start, replication_end), end=pmax(replication_start, replication_end)) %>%
    GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns=T)
  rdc_regions_df = replication_df %>%
    dplyr::inner_join(rdc_df, by=c("annotation_rdc_cluster"="rdc_cluster")) %>%
    dplyr::group_by(annotation_rdc_cluster, replication_chrom) %>%
    dplyr::summarize(rdc_region_start=min(pmin(replication_start, replication_end)), rdc_region_end=max(pmax(replication_start, replication_end)))

  # replication_df %>% dplyr::filter(grepl("RDC_(087|088|089)", annotation_rdc_cluster))
  # rdc_df %>% dplyr::filter(grepl("RDC_(087|088|089)", rdc_cluster))
  # rdc_regions_df %>% dplyr::filter(grepl("RDC_(087|088|089)", annotation_rdc_cluster))

  rdc_regions_ranges = rdc_regions_df %>%
    dplyr::mutate(chrom=replication_chrom, start=rdc_region_start, end=rdc_region_end) %>%
    GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns=T)
  rdc_regions_overlapping_ranges = as(GenomicRanges::coverage(rdc_regions_ranges), "GRanges")
  # rdc_regions_overlapping_ranges = rdc_regions_overlapping_ranges[rdc_regions_overlapping_ranges$score>1]
  rdc_regions_overlapping_ranges$range_id = 1:length(rdc_regions_overlapping_ranges)
  rdc_regions_overlapping_ranges$region_start = GenomicRanges::start(rdc_regions_overlapping_ranges)
  rdc_regions_overlapping_ranges$region_end = GenomicRanges::end(rdc_regions_overlapping_ranges)
  rdc_regions_overlapping_df = as.data.frame(IRanges::mergeByOverlaps(rdc_regions_overlapping_ranges, rdc_regions_ranges))  %>%
    dplyr::select(-dplyr::matches("_ranges\\."), overlaps=score) %>%
    dplyr::group_by(range_id, replication_chrom, region_start, region_end, overlaps) %>%
    dplyr::summarize(main_rdc_cluster = annotation_rdc_cluster[which.max(rdc_region_end)])
  rdc_regions_overlapping_ranges = rdc_regions_overlapping_df %>%
    dplyr::mutate(chrom=replication_chrom, start=region_start, end=region_end) %>%
    GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns=T)
  replication_reduced_df = as.data.frame(IRanges::mergeByOverlaps(replication_ranges, rdc_regions_overlapping_ranges)) %>%
    dplyr::select(-dplyr::matches("_ranges\\."))  %>%
    dplyr::filter(main_rdc_cluster==annotation_rdc_cluster) %>%
    dplyr::select(replication_chrom, replication_start, replication_end, replication_strand)
  replication_reduced_ranges = replication_reduced_df %>%
    dplyr::mutate(seqnames=replication_chrom, strand=replication_strand, start=pmin(replication_start, replication_end), end=pmax(replication_start, replication_end)) %>%
    GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns=T)
  replication_reducedMore_df = as.data.frame(GenomicRanges::reduce(replication_reduced_ranges, drop.empty.ranges=T)) %>%
    dplyr::mutate(replication_chrom=seqnames, replication_strand=strand) %>%
    dplyr::mutate(replication_start=ifelse(replication_strand=="+", start, end)) %>%
    dplyr::mutate(replication_end=ifelse(replication_strand=="+", end, start)) %>%
    dplyr::select(replication_chrom, replication_strand, replication_start, replication_end) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(replication_chrom, replication_start, replication_end)

  readr::write_tsv(replication_reducedMore_df, file="data/replication_reduced_subsets.tsv")
}

annotate_rdc_repliseq = function()
{
  #
  # Load repliseq
  #
  repliseq_df = readr::read_tsv("~/Workspace/Datasets/zhao_bmc_repliseq_2020/preprocessed/repliseq.tsv") %>% dplyr::filter(repliseq_celltype=="npc")
  repliseq_ranges = GenomicRanges::makeGRangesFromDataFrame(repliseq_df %>% dplyr::mutate(seqnames=repliseq_chrom, start=repliseq_start, end=repliseq_end), keep.extra.columns=T)
  repliseqTime_df = readr::read_tsv("~/Workspace/Datasets/zhao_bmc_repliseq_2020/preprocessed/repliseqTime.tsv") %>% dplyr::filter(repliseqTime_celltype=="npc") %>%
    dplyr::mutate(repliseqTime_id=1:n()) %>%
    dplyr::group_by(repliseqTime_celltype, repliseqTime_chrom) %>%
    dplyr::mutate(repliseqTime_smooth=smoother::smth.gaussian(repliseqTime_avg, window=50))
  repliseqTime_ranges = repliseqTime_df %>%
    dplyr::mutate(seqnames=repliseqTime_chrom, start=repliseqTime_start, end=repliseqTime_end) %>%
    GenomicRanges::makeGRangesFromDataFrame(ignore.strand=T, keep.extra.columns=T)


  #
  # Load RDC
  #
  rdc_df = readr::read_tsv("~/Workspace/Datasets/HTGTS/rdc_macs2_mm10.tsv") %>%
    dplyr::mutate(rdc_region_start=rdc_start-4e6, rdc_region_end=rdc_end+4e6)
  rdc_ranges = rdc_df %>%
    dplyr::mutate(seqnames=rdc_chrom, start=rdc_region_start, end=rdc_region_end) %>%
    GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns=T, ignore.strand=T)

  # :RDC_094 RDC_105

  #
  # Join datasets
  #
  #""
  rdc_subsets_df.f = data.frame(rdc_cluster=c("MACS_014"))
  rdc2repliseq_df = as.data.frame(IRanges::mergeByOverlaps(rdc_ranges, repliseq_ranges))  %>%
    dplyr::select(-dplyr::matches("_ranges\\.")) %>%
    dplyr::inner_join(rdc_subsets_df.f, by="rdc_cluster")
  rdc2repliseqTime_df = as.data.frame(IRanges::mergeByOverlaps(rdc_ranges, repliseqTime_ranges))  %>%
    dplyr::select(-dplyr::matches("_ranges\\.")) %>%
    dplyr::inner_join(rdc_subsets_df.f, by="rdc_cluster")


  # Add already annotated values
  replication_df = dplyr::bind_rows(readr::read_tsv("data/replication_subsets.tsv"), readr::read_tsv("data/replication_reduced_subsets.tsv") %>% dplyr::mutate(rdc_cluster="RDC_COMBINED")) %>%
    dplyr::mutate(seqnames=replication_chrom, start=pmin(replication_start, replication_end), end=pmax(replication_start, replication_end), annotated_rdc_cluster=rdc_cluster)
  replication_ranges = replication_df %>%
    dplyr::select(-dplyr::matches("^rdc_cluster$")) %>%
    GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns=T)
  rdc2replication_df = as.data.frame(IRanges::mergeByOverlaps(rdc_ranges, replication_ranges))  %>%
    dplyr::select(-dplyr::matches("_ranges\\.")) %>%
    dplyr::group_by(rdc_cluster) %>%
    dplyr::mutate(rdc_number=as.numeric(as.factor(annotated_rdc_cluster)))

  i = 1
  cl = unique(rdc2repliseq_df$rdc_cluster)[i]
  repliseq_dff = rdc2repliseq_df %>%
    dplyr::filter(rdc_cluster==cl) %>%
    dplyr::mutate(repliseq_value=repliseq_value-min(repliseq_value, na.rm=T)+0.01) %>%
    dplyr::mutate(repliseq_value=repliseq_value^1.5) %>%
    dplyr::group_by(repliseq_start, repliseq_end) %>%
    dplyr::mutate(repliseq_value=repliseq_value/max(repliseq_value)) %>%
    dplyr::ungroup()
    # dplyr::mutate(repliseq_value=repliseq_value/quantile(repliseq_value, 0.95, na.rm=T), repliseq_value=pmin(repliseq_value, 1, na.rm=T))
  repliseqTime_dff = rdc2repliseqTime_df %>%
    dplyr::filter(rdc_cluster==cl)
  rdc_dff = rdc_df %>% dplyr::filter(rdc_cluster==cl)
  rdc2replication_dff = rdc2replication_df %>% dplyr::filter(cl==rdc_cluster)


  p = ggplot(repliseqTime_dff)
  if(nrow(rdc2replication_dff)>0) {
    # pdf("reports/x.pdf", width=6*8.27, height=6*11.69, paper="a4")
    p = p +
      geom_text(aes(x=rdc_region_start, y=-rdc_number, label=annotated_rdc_cluster), data=rdc2replication_dff) +
      geom_segment(aes(x=replication_start, xend=replication_end, y=y-rdc_number, yend=y-rdc_number, color=replication_strand), size=0.5, arrow=grid::arrow(length=unit(3,"pt")), data=rdc2replication_dff %>% dplyr::mutate(y=ifelse(replication_strand=="+", 0.05, -0.05))) +
      geom_vline(aes(xintercept=replication_start), size=0.1, data=rdc2replication_dff) +
      geom_vline(aes(xintercept=replication_end), size=0.1, data=rdc2replication_dff)
    # dev.off()
  }

  p = p +
    geom_rect(aes(ymin=-Inf, ymax=Inf, xmin=rdc_start, xmax=rdc_end), alpha=0.3, fill="#CCCCCC", data=rdc_dff, size=0.04) +
    geom_tile(aes(x=repliseq_start, y=repliseq_fraction), fill=scales::colour_ramp(c("#FFFFFF00", "#FF0000FF"))(repliseq_dff$repliseq_value), data=repliseq_dff) +
    geom_line(aes(x=repliseqTime_start, y=repliseqTime_avg), color="#000000") +
    ggtitle(cl)
  p
  ggannotate::ggannotate()

data = data.frame(x = c(72334443.8282381, 73214478.6042807,
71136098.6012864, 71959960.9448157, 71510581.4847088), y = c(4.41051203917385,
12.0148635231654, 4.67733138948934, 5.94472330348794, 3.27652980033299
), label = c("x", "x", "x", "x", "x"))

  data = rbind(data, data = data.frame(x = c(66804461.1674285, 67005544.8172542,
65634519.932079, 66237770.8815561, 65927005.2409164, 67846440.0801617
), y = c(12.6240392213426, 14.4595311587111, 14.2063598570051,
14.3329455078581, 12.7506248721956, 12.9405033484751), label = c("x",
"x", "x", "x", "x", "x")))

  data = data %>% dplyr::arrange(x)
  with(data, data.frame(rdc_cluster=cl, x_start=x[-length(x)], x_end=x[-1], y_start=y[-length(x)], y_end=y[-1])) %>%
    dplyr::mutate(replication_start=ifelse(y_start>y_end, x_end, x_start), replication_chrom=repliseq_dff$repliseq_chrom[1], replication_end=ifelse(y_start>y_end, x_start, x_end), replication_strand=ifelse(y_start>y_end, "-", "+")) %>%
    dplyr::mutate(replication_start=round(replication_start), replication_end=round(replication_end)) %>%
    dplyr::select(rdc_cluster=rdc_cluster, replication_chrom, replication_start, replication_end, replication_strand) %>%
    dplyr::arrange(pmin(replication_start, replication_end), replication_strand) %>%
    readr::write_tsv("data/x.tsv")

print("x")
}
