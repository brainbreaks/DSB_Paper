setwd("~/Workspace/Everything")

library(readr)
library(dplyr)
library(ggplot2)
library(cowplot)
library(stringr)
devtools::load_all("~/Workspace/breaktools/")



detect_rdc = function()
{
  params = macs2_params(extsize=1e5, exttype="symmetrical", llocal=1e7, minqvalue=0.01, effective_size=1.87e9, maxgap=2e5, minlen=1e5)

  #
  # Load gene annotations
  #
  genes_df = gtf_read('~/Workspace/genomes/mm10/annotation/mm10.ncbiRefSeq.gtf.gz')
  genes_ranges = genes_df %>%
    dplyr::filter(gene_length>=1e5) %>%
    dplyr::mutate(seqnames=gene_chrom, start=gene_start, end=gene_end) %>%
    GenomicRanges::makeGRangesFromDataFrame(, keep.extra.columns=T)

  #
  # Load RDC
  #
  samples_df = tlx_read_samples("~/Workspace/Datasets/HTGTS/samples/All_samples.tsv", "~/Workspace/Datasets/HTGTS/TLX") %>%
    dplyr::filter(!control & (grepl("promoter/enhancer", experiment) & alleles==2 | grepl("concentration", experiment) & concentration==0.4)) %>%
    dplyr::mutate(group=paste0("All (", bait_chrom, ")"))

  # samples_df = tlx_read_samples("~/Workspace/Datasets/HTGTS/samples/All_samples.tsv", "~/Workspace/Datasets/HTGTS/TLX") %>%
  #   dplyr::filter(!control & grepl("concentration", experiment)) %>%
  #   dplyr::mutate(group=paste0(group_short, " (", bait_chrom, ")"))



  tlx_df = tlx_read_many(samples_df, threads=30)
  tlx_df = tlx_remove_rand_chromosomes(tlx_df)
  tlx_df = tlx_extract_bait(tlx_df, bait_size=19, bait_region=12e6)
  tlx_df = tlx_mark_dust(tlx_df)
  tlx_df = tlx_df %>%
    dplyr::group_by(tlx_sample) %>%
    dplyr::filter(dplyr::n()>2000) %>%
    dplyr::ungroup()
  libfactors_df = tlx_libfactors(tlx_df, group="group", normalize_within="group", normalize_between="none", normalization_target="smallest")

  tlxcov_df = tlx_df %>%
    dplyr::filter(B_Rname==bait_chrom & tlx_is_bait_chrom & !tlx_is_bait_junction) %>%
    tlx_coverage(group="group", extsize=params$extsize, exttype=params$exttype, libfactors_df=libfactors_df, ignore.strand=T)

  #
  # Find clusters using MACS3
  #
  macs_results = tlxcov_macs2(tlxcov_df, group="group", params)
  islands_df = macs_results[["islands"]] %>%
    dplyr::inner_join(samples_df %>% dplyr::distinct(group, bait_chrom), by=c("tlx_group"="group"))


  #
  # Remove clusters overlapping with offtarget regions
  #
  islands_ranges = islands_df %>%
    dplyr::mutate(seqnames=island_chrom, start=island_start, end=island_end) %>%
    GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns=T)
  offtarget_df = readr::read_tsv("~/Workspace/Datasets/HTGTS/offtargets_pnas_mm10.tsv")
  offtarget_ranges = offtarget_df %>%
    dplyr::rename(offtarget_bait_chrom="bait_chrom") %>%
    dplyr::mutate(offtarget_start=offtarget_start-5e4, offtarget_end=offtarget_end+5e4) %>%
    GenomicRanges::makeGRangesFromDataFrame(seqnames.field="offtarget_chrom", start.field="offtarget_start", end.field="offtarget_end", ignore.strand=T, keep.extra.columns=T)
  islands_df = as.data.frame(leftJoinByOverlaps(islands_ranges, offtarget_ranges)) %>%
    dplyr::select(-dplyr::matches("_ranges\\.")) %>%
    dplyr::mutate(island_is_offtarget=!is.na(offtarget_bait_chrom) & offtarget_bait_chrom==bait_chrom) %>%
    dplyr::filter(!island_is_offtarget) %>%
    dplyr::arrange(as.numeric(gsub("chr", "", island_chrom)), island_start) %>%
    dplyr::mutate(island_name=stringr::str_glue("MACS_{stringr::str_pad(i, 3, pad='0')}", i=1:dplyr::n()))

  #
  # Write final results
  #
  islands_ranges = islands_df %>%
    dplyr::mutate(seqnames=island_chrom, start=island_start, end=island_end) %>%
    GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns=T)
  island_rdc_df = as.data.frame(leftJoinByOverlaps(islands_ranges, genes_ranges)) %>%
    dplyr::select(-dplyr::matches("_ranges\\.")) %>%
    dplyr::arrange(dplyr::desc(tidyr::replace_na(gene_length, 0))) %>%
    dplyr::distinct(island_name, .keep_all=T)  %>%
    dplyr::arrange(as.numeric(gsub("[^0-9]*(\\d+).*", "\\1", island_name))) %>%
    dplyr::mutate(gene_strand=tidyr::replace_na(gene_strand, "."), gene_id=tidyr::replace_na(gene_id, "-")) %>%
    dplyr::mutate(rdc_group=1, rdc_cluster_display=paste0(island_chrom, ":", island_name, gene_strand)) %>%
    dplyr::select(rdc_chrom=island_chrom, rdc_start=island_start, rdc_end=island_end, rdc_cluster=island_name,	rdc_group, rdc_strand=gene_strand, rdc_gene=gene_id, rdc_cluster_display, rdc_stn=island_snr, rdc_qvalue_log10=island_summit_qvalue, rdc_pileup=island_summit_abs)

  readr::write_tsv(island_rdc_df, "~/Workspace/Datasets/HTGTS/rdc_macs2_mm10.tsv", col_names=T)
  island_rdc_df %>%
    dplyr::select(rdc_chrom, rdc_start, rdc_end, rdc_cluster_display, rdc_qvalue_log10, rdc_strand) %>%
    readr::write_tsv(file="~/Workspace/Datasets/HTGTS/rdc_macs2_mm10.bed", col_names=F)

  macs_results[["qvalues"]] %>%
    dplyr::select(qvalue_chrom, qvalue_start, qvalue_end, qvalue_score) %>%
    readr::write_tsv(file="~/Workspace/Datasets/HTGTS/rdc_macs2_mm10_qvalues.bdg", col_names=F)

  tlxcov_df %>%
    dplyr::select(tlxcov_chrom, tlxcov_start, tlxcov_end, tlxcov_pileup) %>%
    readr::write_tsv(file="~/Workspace/Datasets/HTGTS/rdc_macs2_mm10.bdg", col_names=F)


  islands_df %>% dplyr::filter(grepl("MACS_003", island_name))
}