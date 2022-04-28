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
  samples_df = tlx_read_samples("~/Workspace/Datasets/HTGTS/samples/All_samples.tsv", "~/Workspace/Datasets/HTGTS") %>%
    dplyr::filter(!control & grepl("concentration", experiment))

  tlx_all_df = tlx_read_many(samples_df, threads=30)
  tlx_df = tlx_all_df %>%
    tlx_extract_bait(bait_size=19, bait_region=12e6) %>%
    tlx_remove_rand_chromosomes() %>%
    tlx_mark_dust() %>%
    tlx_calc_copynumber(bowtie_index="~/Workspace/genomes/mm10/mm10", max_hits=100, threads=24) %>%
    dplyr::filter(tlx_copynumber==1 & !tlx_duplicated) %>%
    dplyr::group_by(tlx_sample) %>%
    dplyr::filter(dplyr::n()>5000) %>%
    dplyr::ungroup()

  libfactors_bait_df = tlx_libfactors(tlx_df %>% dplyr::mutate(tlx_group=bait_chrom), group="group", normalize_within="group", normalize_between="none", normalization_target="smallest")

  #
  # Indentify offtarget candidates
  #
  tlx_offtarget_df = tlx_df %>% dplyr::mutate(tlx_group=bait_chrom)
  offtargets_params = macs2_params(extsize=50, exttype="opposite", llocal=1e7, minqvalue=0.05, effective_size=1.87e9, maxgap=1e3, minlen=10)
  tlxcov_offtargets_df = tlx_offtarget_df %>%
    dplyr::filter(tlx_is_bait_chrom & !tlx_is_bait_junction) %>%
    tlx_coverage(group="group", extsize=offtargets_params$extsize, exttype=offtargets_params$exttype, libfactors_df=libfactors_bait_df, ignore.strand=T)
  macs_offtargets = tlxcov_macs2(tlxcov_offtargets_df, group="group", offtargets_params)

  # Export debuging info
  tlx_write_bed(tlx_offtarget_df, "reports/APH_concentration/offtargets", "all", mode="alignment")
  macs_offtargets$islands %>%
    dplyr::mutate(score=1) %>%
    dplyr::select(island_chrom, island_start, island_end, score) %>%
    readr::write_tsv("reports/APH_concentration/offtargets-islands.bed", col_names=F)
  macs_offtargets$qvalues %>%
    dplyr::select(qvalue_chrom, qvalue_start, qvalue_end, qvalue_score) %>%
    readr::write_tsv("reports/APH_concentration/offtargets-qvalues.bedgraph", col_names = F)
  tlxcov_write_bedgraph(tlxcov_offtargets_df, "reports/APH_concentration/offtargets", "all")


  #
  # Find and mark offtargets
  #
  offtargets_df = readr::read_tsv("~/Workspace/Datasets/HTGTS/offtargets_pnas_mm10.tsv")
  offtargets_extended_df = dplyr::bind_rows(
    offtargets_df,
    macs_offtargets$islands %>%
      dplyr::inner_join(offtargets_df %>% dplyr::distinct(offtarget_bait_chrom, offtarget_bait_start, offtarget_bait_end, offtarget_bait_strand), by=c("island_chrom"="offtarget_bait_chrom")) %>%
      dplyr::mutate(offtarget_is_primary=F) %>%
      dplyr::select(
        offtarget_bait_chrom=island_chrom, offtarget_bait_start, offtarget_bait_end, offtarget_bait_strand,
        offtarget_chrom=island_chrom, offtarget_start=island_start, offtarget_end=island_end, offtarget_is_primary)
  )

  tlx_clean_df = tlx_df %>%
    tlx_mark_offtargets(offtargets_extended_df, offtarget_region=1e5, bait_region=1e4) %>%
    dplyr::filter(!tlx_is_offtarget & tlx_is_bait_chrom)

  #
  # Use all concentrations to call RDC
  #
  params_clean = macs2_params(extsize=5e4, exttype="symmetrical", llocal=1e7, minqvalue=0.01, effective_size=1.87e9, maxgap=2e5, minlen=1e5)
  tlxcov_clean_df = tlx_clean_df %>%
    dplyr::mutate(tlx_group=bait_chrom) %>%
    tlx_coverage(group="group", extsize=params_clean$extsize, exttype=params_clean$exttype, libfactors_df=libfactors_bait_df, ignore.strand=T)
  macs_clean = tlxcov_macs2(tlxcov_clean_df, group="group", params_clean)

  # Export debuging info
  macs_clean$islands %>% df2ranges(island_chrom, island_start, island_end) %>% rtracklayer::export.bed("reports/APH_concentration/all-rdc.bed")
  tlxcov_write_bedgraph(tlxcov_clean_df, "reports/APH_concentration/all", "all")
  tlx_clean_df %>%
    dplyr::mutate(score=1) %>%
    dplyr::select(Rname, Rstart, Rend, Qname, score, tlx_strand) %>%
    readr::write_tsv("reports/APH_concentration_allnorm/all-reference.bed", col_names=F)

  # table(island=macs_clean[["islands"]]$island_chrom, group=macs_clean[["islands"]]$tlx_group)


  #
  # Find clusters using MACS3
  #
  # macs_results = tlxcov_macs2(tlxcov_all_df, group="all", params)
  # islands_ranges = macs_results[["islands"]] %>% df2ranges(island_chrom, island_start, island_end)
  # macs_results[["qvalues"]] %>%
  #   dplyr::filter(qvalue_chrom=="chr5") %>%
  #   dplyr::select(qvalue_chrom, qvalue_start, qvalue_end, qvalue_score) %>%
  #   readr::write_tsv("reports/APH_concentration_allnorm/islands.bedgraph", col_names = F)
  # rdc_df = macs_results[["islands"]] %>%
  #   # dplyr::filter(island_snr>=10) %>%
  #   dplyr::arrange(as.numeric(gsub("chr", "", island_chrom)), island_start) %>%
  #   dplyr::mutate(island_name=stringr::str_glue("MACS_{stringr::str_pad(i, 3, pad='0')}", i=1:dplyr::n())) %>%
  #   dplyr::select(rdc_chrom=island_chrom, rdc_start=island_start, rdc_end=island_end, rdc_cluster=island_name)
  # rtracklayer::export.bed(islands_ranges, con="reports/APH_concentration_allnorm/islands.bed")

  #
  # Load RDC
  #
  # rdc_df = readr::read_tsv("~/Workspace/Datasets/HTGTS/rdc_macs2_mm10.tsv")
  rdc_df = macs_clean$islands %>%
    # dplyr::filter(island_snr>=10) %>%
    dplyr::select(rdc_chrom=island_chrom, rdc_start=island_start, rdc_end=island_end, rdc_cluster=island_name)

  #
  # Find RDC center
  #
  rdc_ranges = rdc_df %>% df2ranges(rdc_chrom, (rdc_start+rdc_end)/2-2e6, (rdc_start+rdc_end)/2+2e6)
  tlxcov_clean_ranges = tlxcov_clean_df %>% df2ranges(tlxcov_chrom, (tlxcov_start+tlxcov_end)/2-2e6, (tlxcov_start+tlxcov_end)/2+2e6)
  rdc2tlxcov_df = innerJoinByOverlaps(rdc_ranges, tlxcov_clean_ranges)
  rdc_center_ranges = rdc2tlxcov_df %>%
    dplyr::group_by(rdc_cluster, rdc_chrom) %>%
    dplyr::summarize(rdc_center=stats::weighted.mean(tlxcov_start/2+tlxcov_end/2, tlxcov_pileup*(tlxcov_end-tlxcov_start), na.rm=T), rdc_region_start=rdc_center-1e6, rdc_region_end=rdc_center+1e6, tlxcov_pileup=max(tlxcov_pileup)) %>%
    dplyr::filter(tlxcov_pileup>=5) %>%
    df2ranges(rdc_chrom, rdc_region_start, rdc_region_end)

  #
  # Each sample separately
  #
  params_baitconcentration = macs2_params(extsize=100e3, exttype="symmetrical", llocal=1e7, minqvalue=0.01, effective_size=1.87e9, maxgap=2e5, minlen=1e5)
  tlx_baitconcentration_df = tlx_clean_df %>% dplyr::mutate(tlx_group=paste0(tlx_group, " (", bait_chrom, ")"))
  libfactors_baitconcentration_df = tlx_libfactors(tlx_baitconcentration_df, group="group", normalize_within="group", normalize_between="none", normalization_target="smallest")
  tlxcov_baitconcentration_strand_df = tlx_clean_df %>%
    tlx_coverage(group="group", extsize=params_baitconcentration$extsize, exttype=params_baitconcentration$exttype, libfactors_df=libfactors_baitconcentration_df, ignore.strand=F)
  tlxcov_baitconcentration_strand_ranges = tlxcov_baitconcentration_strand_df %>% df2ranges(tlxcov_chrom, tlxcov_start, tlxcov_end)

  #
  # Write strand-stratified bedgraph
  #
  tlxcov_baitconcentration_strand_df %>%
    dplyr::mutate(tlx_group=gsub(" \\(.*", "", tlx_group)) %>%
    tlxcov_write_bedgraph(path="reports/APH_concentration/concentration", group="group")

  x = rdc_ranges %>%
    innerJoinByOverlaps(tlxcov_baitconcentration_strand_ranges) %>%
    dplyr::filter(rdc_chrom=="chr6" & rdc_start<77e6 & rdc_end>77e6 & tlx_group=="APH 0.2 uM 96h")

  devtools::load_all('~/Workspace/breaktools/')
  ccs = rdc_ranges %>%
    innerJoinByOverlaps(tlxcov_baitconcentration_strand_ranges) %>%
    dplyr::group_by(rdc_cluster, rdc_chrom, tlx_group) %>%
    summarize(tlx_strand_crosscorrelation(dplyr::cur_data(), step=5000)) %>%
    dplyr::filter(crosscorrelation_value>0.5)
  ggplot(ccs, aes(y=crosscorrelation_rellag, x=1, fill=tlx_group)) +
    geom_boxplot(outlier.shape = NA, outlier.alpha = 0) +
    geom_point(aes(text=rdc_cluster), position=position_jitter())


  plot(y=ccf.sense$acf, x=ccf.sense$lag*100)

  rdc2tlxcov_weight_df = rdc_ranges %>%
    innerJoinByOverlaps(tlxcov_baitconcentration_strand_ranges) %>%
    dplyr::group_by(rdc_cluster, rdc_chrom, tlx_strand, tlx_group) %>%
    dplyr::summarize(tlxcov_strand_center=stats::weighted.mean(tlxcov_start/2+tlxcov_end/2, tlxcov_pileup*(tlxcov_end-tlxcov_start), na.rm=T)) %>%
    dplyr::group_by(rdc_cluster, rdc_chrom, tlx_group) %>%
    # dplyr::summarise(tlxcov_strand_shift=tlxcov_strand_center[tlx_strand=="+"]-tlxcov_strand_center[tlx_strand=="-"]) %>%
    dplyr::summarise(tlxcov_strand_shift=abs(diff(tlxcov_strand_center))) %>%
    dplyr::group_by(rdc_cluster, rdc_chrom) %>%
    dplyr::mutate(tlxcov_strand_relshift=tlxcov_strand_shift/tlxcov_strand_shift[grepl("APH 0.2", tlx_group)]) %>%
    dplyr::mutate(concentration=gsub("APH ([0-9.]+) .*", "\\1", tlx_group), chrom=gsub(".*(chr[0-9]+).*", "\\1", tlx_group))
  # library(ggiraph)
  g = ggplot(rdc2tlxcov_weight_df, aes(y=tlxcov_strand_shift, x=chrom, fill=concentration)) +
    geom_boxplot(outlier.shape = NA, outlier.alpha = 0) +
    geom_point(aes(text=rdc_cluster), position=position_jitter())
  plotly::ggplotly(g)

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