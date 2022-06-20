setwd("~/Workspace/DSB_Paper")

library(readr)
library(dplyr)
library(ggplot2)
library(cowplot)
library(stringr)
library(randomcoloR)
library(ComplexHeatmap)
devtools::load_all("~/Workspace/breaktools/")



detect_rdc = function()
{
  debug=T
  dir.create("reports/detect_rdc", recursive=T)

  #
  # Load offtargets
  #
  offtargets_df = readr::read_tsv("data/offtargets_dkfz.tsv")


  #
  # Load gene annotations
  #
  genes_df = gtf_read('~/Workspace/genomes/mm10/annotation/mm10.ncbiRefSeq.gtf.gz')
  genes_ranges = genes_df %>%
    dplyr::filter(gene_length>=1e5) %>%
    df2ranges(gene_chrom, gene_start, gene_end)

  #
  # Load RDC
  #
  samples_df = tlx_read_samples("~/Workspace/Datasets/HTGTS/samples/All_samples.tsv", "~/Workspace/Datasets/HTGTS") %>%
    dplyr::filter(tlx_exists & celltype=="NPC" & organism=="mouse" & sample!="VI035" & !control & (
      grepl("(Csmd1|Ctnna2|Nrxn1) promoter/enhancer", experiment) |
      grepl("concentration", experiment) & concentration==0.4)
    )

  tlx_all_df = tlx_read_many(samples_df, threads=10) %>%
    tlx_extract_bait(bait_size=19, bait_region=12e6) %>%
    tlx_calc_copynumber(bowtie2_index="~/Workspace/genomes/mm10/mm10", max_hits=100, threads=24) %>%
    tlx_mark_offtargets(offtargets_df, offtarget_region=1e6, bait_region=1e4)

  # tlx_all_df = tlx_all_df.bck %>%
  #   dplyr::mutate(tlx_group=tlx_group)
    # dplyr::group_by(tlx_sample) %>%
    # dplyr::filter(dplyr::n()>5000) %>%
    # dplyr::ungroup()


  libfactors_df = tlx_all_df %>%
    tlx_libsizes()
    # tlx_libsizes(within=c("tlx_group", "tlx_control"), between="tlx_group") %>%
    # tlx_libfactors_within(min(library_size)/library_size)
  # libfactors_df = tlx_libfactors(tlx_df, group="group", normalize_within="group", normalize_between="none", normalization_target="smallest")

  #
  # Evaluate extsize
  #
  if(F) {
    pdf("reports/extsize_selection.pdf", width=11.69, height=8.27, paper="a4r")
    extsize_df = data.frame(extsize=c(seq(100, 500, 100), seq(1000, 9000, 1000), seq(10000, 90000, 10000), seq(1e5, 3e5, 5e5))) %>%
      dplyr::rowwise() %>%
      dplyr::do((function(z){
        zz<<-z
        params_rdc = macs2_params(extsize=z$extsize, exttype="symmetrical", llocal=1e7, minpvalue=0.01, effective_size=1.87e9, maxgap=2e5, minlen=2e5, baseline=2)
        tlxcov_rdc_df = tlx_rdc_df %>%
          dplyr::filter(!tlx_is_bait_junction & (!grepl("prom", group) | !tlx_is_bait_chrom)) %>%
          tlx_coverage(group="group", extsize=params_rdc$extsize, exttype=params_rdc$exttype, libfactors_df=libfactors_df, ignore.strand=T)
        macs_rdc = tlxcov_macs2(tlxcov_rdc_df, group="group", params=params_rdc)
        macs_rdc$islands %>%
          tidyr::crossing(as.data.frame(z))
      })(.)) %>%
      dplyr::ungroup()
    extsize_sumdf = extsize_df %>%
      dplyr::group_by(extsize, tlx_group) %>%
      dplyr::summarize(count=dplyr::n(), width=mean(island_end-island_start))
    ggplot(extsize_sumdf) +
      geom_line(aes(x=extsize, y=count)) +
      geom_vline(xintercept=5e4, color="#FF0000") +
      geom_smooth(aes(x=extsize, y=count)) +
      labs(y="RDC count") +
      facet_wrap(~tlx_group, scales="free")
    ggplot(extsize_sumdf) +
      geom_line(aes(x=extsize, y=width)) +
      geom_vline(xintercept=5e4, color="#FF0000") +
      geom_smooth(aes(x=extsize, y=width)) +
      labs(y="RDC width") +
      facet_wrap(~tlx_group, scales="free")
    dev.off()
  }

  #
  # Detect RDC
  #
  tlx_rdc_df = tlx_all_df1 %>%
    tlx_remove_rand_chromosomes() %>%
    dplyr::filter(tlx_copynumber==1 & !tlx_duplicated & !tlx_is_bait_junction) %>%
    dplyr::mutate(tlx_group=ifelse(tlx_is_bait_chrom, "Intra", "Inter"))


  params_rdc = macs2_params(extsize=5e4, exttype="symmetrical", llocal=1e7, minpvalue=0.01, effective_size=1.87e9, maxgap=2e5, minlen=2e5, baseline=2)
  # params_rdc = macs2_params(extsize=5e4, exttype="symmetrical", llocal=1e7, minpvalue=0.01, effective_size=1.87e9, maxgap=0, minlen=1e3, baseline=2)
  tlxcov_rdc_df = tlx_rdc_df %>%
    dplyr::filter(!tlx_is_offtarget & !tlx_is_bait_junction & (!grepl("prom", group) | !tlx_is_bait_chrom)) %>%
    tlx_coverage(group="group", extsize=params_rdc$extsize, exttype=params_rdc$exttype, libfactors_df=libfactors_df, ignore.strand=T)

  macs_rdc = tlxcov_macs2(tlxcov_rdc_df, group="group", params=params_rdc)
  if(debug)
  {
    macs_rdc$qvalues %>%
      dplyr::group_by(tlx_group) %>%
      dplyr::do((function(df){
        df %>%
          dplyr::select(qvalue_chrom, qvalue_start, qvalue_end, qvalue_score) %>%
          readr::write_tsv(paste0("reports/detect_rdc/qvalues-", df$tlx_group[1], ".bedgraph"), col_names = F)
      })(.))
    macs_rdc$islands %>%
      dplyr::mutate(strand="*") %>%
      dplyr::group_by(tlx_group) %>%
      dplyr::do((function(df){
        df %>%
          dplyr::select(island_chrom, island_extended_start, island_extended_end, island_name, island_baseline, strand) %>%
          readr::write_tsv(paste0("reports/detect_rdc/extislands-", df$tlx_group[1], ".bed"), col_names=F)
        df %>%
          dplyr::select(island_chrom, island_start, island_end, island_name, island_baseline, strand) %>%
          readr::write_tsv(paste0("reports/detect_rdc/islands-", df$tlx_group[1], ".bed"), col_names=F)
      })(.))

    tlxcov_rdc_df %>% tlxcov_write_bedgraph(path="reports/detect_rdc/bedgraph", group="group")
    tlx_rdc_df %>% tlx_write_bed(path="reports/detect_rdc/bed", group="group")

    # tlx_rdc_df %>%
    #   tlx_coverage(group="group", extsize=params_rdc$extsize, exttype=params_rdc$exttype, libfactors_df=libfactors_df, ignore.strand=F) %>%
    #   tlxcov_write_bedgraph(path="reports/pdetect_rdc/bedgraph", group="group")
  }


  # tlxcov_df = tlx_df %>%
  #   dplyr::filter(B_Rname==bait_chrom & tlx_is_bait_chrom & !tlx_is_bait_junction) %>%
  #   tlx_coverage(group="group", extsize=params$extsize, exttype=params$exttype, libfactors_df=libfactors_df, ignore.strand=T)

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