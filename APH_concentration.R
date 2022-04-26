setwd("~/Workspace/DSB_Paper")

library(readr)
library(dplyr)
library(ggplot2)
library(cowplot)
library(stringr)
devtools::load_all('~/Workspace/breaktools/')

APH_concentration = function()
{
  #
  # Read TLX
  #
  params = macs2_params(extsize=5e4, exttype="symmetrical", llocal=1e7, minqvalue=0.01, effective_size=1.87e9, maxgap=2e5, minlen=1e5)
  samples_df = tlx_read_samples("~/Workspace/Datasets/HTGTS/samples/All_samples.tsv", "~/Workspace/Datasets/HTGTS") %>%
    dplyr::filter(!control & grepl("concentration", experiment)) %>%
    dplyr::mutate(group=paste0(group_short, " (", bait_chrom, ")"))

  tlx_df = tlx_read_many(samples_df, threads=30)
  tlx_df = tlx_remove_rand_chromosomes(tlx_df)
  tlx_df = tlx_extract_bait(tlx_df, bait_size=19, bait_region=12e6)
  tlx_df = tlx_mark_dust(tlx_df)
  tlx_df = tlx_df %>%
    dplyr::group_by(tlx_sample) %>%
    dplyr::filter(dplyr::n()>2000) %>%
    dplyr::ungroup()


  #
  # All coverages separately
  #
  libfactors_df = tlx_libfactors(tlx_df%>% dplyr::filter(tlx_is_bait_chrom & !tlx_is_bait_junction), group="group", normalize_within="group", normalize_between="group", normalization_target="smallest")
  tlxcov_all_df = tlx_df %>%
    dplyr::filter(tlx_is_bait_chrom & !tlx_is_bait_junction) %>%
    tlx_coverage(group="all", extsize=params$extsize, exttype=params$exttype, libfactors_df=libfactors_df, ignore.strand=T)
  tlxcov_all_ranges = tlxcov_all_df %>% df2ranges(tlxcov_chrom, tlxcov_start, tlxcov_end)

  #
  # Load RDC
  #
  rdc_df = readr::read_tsv("~/Workspace/Datasets/HTGTS/rdc_macs2_mm10.tsv")
  rdc_ranges = rdc_df %>%
    dplyr::mutate(rdc_pos=(rdc_start+rdc_end)/2) %>%
    dplyr::mutate(rdc_region_start=rdc_pos-2e6, rdc_region_end=rdc_pos+2e6) %>%
    df2ranges(rdc_chrom, rdc_region_start, rdc_region_end)
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
  tlxcov_write_bedgraph(tlxcov_strand_df %>% dplyr::mutate(tlx_group=gsub(" \\(.*", "", tlx_group)), path="treatment_allnorm", group="group")
  tlxcov_strand_ranges = tlxcov_strand_df %>% df2ranges(tlxcov_chrom, tlxcov_start, tlxcov_end)


  rdc2tlxcov_df = innerJoinByOverlaps(rdc_center_range, tlxcov_strand_ranges)
  rdc2tlxcov_weight_df = rdc2tlxcov_df %>%
    dplyr::group_by(rdc_cluster, tlx_strand, tlx_group) %>%
    dplyr::summarize(tlxcov_wcenter=stats::weighted.mean(tlxcov_start/2+tlxcov_end/2, tlxcov_pileup*(tlxcov_end-tlxcov_start), na.rm=T)) %>%
    dplyr::group_by(rdc_cluster, tlx_group) %>%
    dplyr::summarise(tlxcov_wcenter_shift=max(tlxcov_wcenter)-min(tlxcov_wcenter)) %>%
    dplyr::group_by(rdc_cluster) %>%
    dplyr::mutate(tlxcov_wcenter_relshift=tlxcov_wcenter_shift) %>%
    dplyr::mutate(concentration=gsub("APH ([0-9.]+) .*", "\\1", tlx_group), chrom=gsub(".*(chr[0-9]+).*", "\\1", tlx_group))
  ggplot(rdc2tlxcov_weight_df, aes(y=tlxcov_wcenter_relshift, x=chrom, fill=concentration)) +
    geom_boxplot() +
    geom_jitter()



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