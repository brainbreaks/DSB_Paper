setwd("~/Workspace/DSB_Paper")

library(readr)
library(dplyr)
library(ggplot2)
library(cowplot)
library(stringr)
devtools::load_all('~/Workspace/breaktools/')

APH_concentration = function()
{
  params = macs2_params(extsize=5e4, exttype="symmetrical", llocal=1e7, minqvalue=0.01, effective_size=1.87e9, maxgap=2e5, minlen=1e5)

  #
  # Find and mark offtargets
  #
  offtargets_df = readr::read_tsv("~/Workspace/Datasets/HTGTS/offtargets_pnas_mm10.tsv")
  offtargets_df %>% df2ranges(offtarget_chrom, offtarget_start, offtarget_end) %>%
  rtracklayer::export.bed( con="reports/APH_concentration_allnorm/islands.bed")

  #
  # Read TLX
  #
  samples_df = tlx_read_samples("~/Workspace/Datasets/HTGTS/samples/All_samples.tsv", "~/Workspace/Datasets/HTGTS") %>%
    dplyr::filter(!control & grepl("concentration", experiment)) %>%
    dplyr::mutate(group=paste0(group_short, " (", bait_chrom, ")"))

  tlx_all_df = tlx_read_many(samples_df, threads=30)
  tlx_df = tlx_all_df %>%
    tlx_extract_bait(bait_size=19, bait_region=12e6) %>%
    tlx_remove_rand_chromosomes() %>%
    tlx_mark_dust()
  tlx_df = tlx_df %>% tlx_mark_offtargets(offtargets_df, offtarget_region=1e5, bait_region=1000) %>%
    dplyr::group_by(tlx_sample) %>%
    dplyr::filter(dplyr::n()>2000) %>%
    dplyr::ungroup()


  #
  # All coverages separately
  #
  libfactors_df = tlx_libfactors(tlx_all_df, group="group", normalize_within="group", normalize_between="all", normalization_target="smallest")
  tlxcov_all_df = tlx_df %>%
    dplyr::filter(tlx_is_bait_chrom & !tlx_is_bait_junction & !tlx_is_offtarget) %>%
    tlx_coverage(group="all", extsize=params$extsize, exttype=params$exttype, libfactors_df=libfactors_df, ignore.strand=T)
  tlxcov_all_ranges = tlxcov_all_df %>% df2ranges(tlxcov_chrom, tlxcov_start, tlxcov_end)

  #
  # Find clusters using MACS3
  #
  macs_results = tlxcov_macs2(tlxcov_all_df, group="all", params)
  islands_ranges = macs_results[["islands"]] %>% df2ranges(island_chrom, island_start, island_end)
  rdc_df = macs_results[["islands"]] %>%
    dplyr::filter(island_snr>=10) %>%
    dplyr::arrange(as.numeric(gsub("chr", "", island_chrom)), island_start) %>%
    dplyr::mutate(island_name=stringr::str_glue("MACS_{stringr::str_pad(i, 3, pad='0')}", i=1:dplyr::n())) %>%
    dplyr::select(rdc_chrom=island_chrom, rdc_start=island_start, rdc_end=island_end, rdc_cluster=island_name)
  rtracklayer::export.bed(islands_ranges, con="reports/APH_concentration_allnorm/islands.bed")


  #
  # Load RDC
  #
  rdc_df = readr::read_tsv("~/Workspace/Datasets/HTGTS/rdc_macs2_mm10.tsv")
  rdc_ranges = rdc_df %>%
    dplyr::mutate(rdc_pos=(rdc_start+rdc_end)/2) %>%
    dplyr::mutate(rdc_region_start=rdc_pos-2e6, rdc_region_end=rdc_pos+2e6) %>%
    df2ranges(rdc_chrom, rdc_region_start, rdc_region_end)

  #
  # Find RDC center
  #
  rdc2tlxcov_df = innerJoinByOverlaps(rdc_ranges, tlxcov_all_ranges)
  rdc_center_df = rdc2tlxcov_df %>%
    dplyr::group_by(rdc_cluster, rdc_chrom) %>%
    dplyr::summarize(rdc_center=stats::weighted.mean(tlxcov_start/2+tlxcov_end/2, tlxcov_pileup*(tlxcov_end-tlxcov_start), na.rm=T), rdc_region_start=rdc_center-1e6, rdc_region_end=rdc_center+1e6, tlxcov_pileup=max(tlxcov_pileup)) %>%
    dplyr::filter(tlxcov_pileup>=5)
  rdc_center_range = rdc_center_df %>% df2ranges(rdc_chrom, rdc_region_start, rdc_region_end)

  #
  # Each sample separately
  #
  libfactors_df = tlx_libfactors(tlx_df, group="group", normalize_within="group", normalize_between="all", normalization_target="smallest")
  tlxcov_strand_df = tlx_df %>%
    dplyr::filter(tlx_is_bait_chrom & !tlx_is_bait_junction) %>%
    tlx_coverage(group="group", extsize=params$extsize, exttype=params$exttype, libfactors_df=libfactors_df, ignore.strand=F)
  tlxcov_strand_ranges = tlxcov_strand_df %>% df2ranges(tlxcov_chrom, tlxcov_start, tlxcov_end)

  #
  # Write strand-stratified bedgraph
  #
  tlxcov_write_bedgraph(tlxcov_strand_df %>% dplyr::mutate(tlx_group=gsub(" \\(.*", "", tlx_group)), path="reports/APH_concentration_allnorm", group="group")


  rdc2tlxcov_weight_df = rdc_center_range %>%
    innerJoinByOverlaps(tlxcov_strand_ranges) %>%
    dplyr::group_by(rdc_cluster, tlx_strand, tlx_group) %>%
    dplyr::summarize(tlxcov_strand_center=stats::weighted.mean(tlxcov_start/2+tlxcov_end/2, tlxcov_pileup*(tlxcov_end-tlxcov_start), na.rm=T)) %>%
    dplyr::group_by(rdc_cluster, tlx_group) %>%
    # dplyr::summarise(tlxcov_strand_shift=tlxcov_strand_center[tlx_strand=="+"]-tlxcov_strand_center[tlx_strand=="-"]) %>%
    dplyr::summarise(tlxcov_strand_shift=abs(diff(tlxcov_strand_center))) %>%
    dplyr::group_by(rdc_cluster) %>%
    dplyr::mutate(tlxcov_strand_relshift=tlxcov_strand_shift/max(tlxcov_strand_shift)) %>%
    dplyr::mutate(concentration=gsub("APH ([0-9.]+) .*", "\\1", tlx_group), chrom=gsub(".*(chr[0-9]+).*", "\\1", tlx_group))
  ggplot(rdc2tlxcov_weight_df, aes(y=tlxcov_strand_shift, x=chrom, fill=concentration)) +
    geom_boxplot() +
    geom_jitter()
    coord_cartesian(ylim=c(-5, 5))



  tlx_df.f = tlx_df %>%
    dplyr::filter(tlx_is_bait_chrom & !tlx_is_bait_junction) %>%
    dplyr::group_by(tlx_sample, Rname) %>%
    dplyr::mutate(dbscan_cluster=dbscan::dbscan(matrix(Junction), minPts=20, eps=100)$cluster) %>%
    dplyr::ungroup() %>%
    dplyr::filter(dbscan_cluster==0)
  tlx_write_bedgraph(tlx_df.f, path="reports/concentration-bedgraph-5e5", group="group", exttype="symmetrical", split_strand=T)

  tlx_df %>%
    dplyr::group_by(tlx_sample) %>%
    dplyr::filter(dplyr::n() >= 2000) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(dbscan_cluster=dbscan::dbscan(matrix(Junction), minPts=20, eps=100)$cluster) %>%
    dplyr::filter(dbscan_cluster==0) %>%
    tlx_write_bedgraph(path="reports/bedgraph-1e5", group="group", exttype="symmetrical", extsize=1e5)

}