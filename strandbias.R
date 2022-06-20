library(readr)
library(dplyr)
library(ggplot2)
library(cowplot)
library(stringr)
devtools::load_all('~/Workspace/breaktools/')

#
main = function()
{
  baits_df = readr::read_tsv("~/Workspace/Datasets/HTGTS/wei_pnas2018_baits.tsv")

  #
  # Read TLX files
  #
  samples_df = tlx_read_samples("~/Workspace/Datasets/HTGTS/samples/All_samples.tsv", "~/Workspace/Datasets/HTGTS") %>%
    dplyr::filter(grepl("(Ctnna2|Nrxn1) promoter/enhancer|Wei|Tena|concentration", experiment) & !control & tlx_exists)

  tlx_df = tlx_read_many(samples_df, threads=20) %>%
    tlx_extract_bait(bait_size=19, bait_region=12e6) %>%
    tlx_remove_rand_chromosomes() %>%
    tlx_calc_copynumber(bowtie2_index="~/Workspace/genomes/mm10/mm10", max_hits=100, threads=24) %>%
    dplyr::filter(tlx_copynumber==1 & !tlx_duplicated & Rname!="chrY") %>%
    dplyr::group_by(experiment, run, tlx_sample) %>%
    dplyr::filter(dplyr::n()>1000) %>%
    dplyr::ungroup()


  #
  # Automatically find RDC
  #
  params = macs2_params(extsize=5e4, exttype="symmetrical", llocal=1e7, minqvalue=0.01, effective_size=1.87e9, maxgap=1e6, minlen=1e5)
  tlx_rdc_df = tlx_df %>% dplyr::filter(!tlx_control)
  libfactors_rdc_df = tlx_rdc_df %>%
    dplyr::mutate(tlx_group="all") %>%
    tlx_libfactors(normalize_within="group", normalize_between="none", normalization_target="min")
  tlxcov_rdc_df = tlx_rdc_df %>%
    dplyr::filter(!tlx_is_bait_junction) %>%
    dplyr::mutate(tlx_group=tlx_is_bait_chrom) %>%
    tlx_coverage(group="group", extsize=params$extsize, exttype=params$exttype, libfactors_df=libfactors_rdc_df, ignore.strand=T)
  macs_rdc = tlxcov_macs2(tlxcov_rdc_df, group="group", params)
  # tlxcov_write_bedgraph(tlxcov_rdc_df, path="reports/promoter-enhancer_deletion/pool", group="all")
  # macs_rdc$islands %>% dplyr::mutate(strand="*") %>% dplyr::select(island_chrom, island_start, island_end, island_name, island_baseline, strand) %>%
  #   readr::write_tsv("reports/promoter-enhancer_deletion/pool_macs.bed", col_names=F)

  tlx_rdc_df %>%
    dplyr::group_by(tlx_is_bait_junction, tlx_is_bait_chrom, Rname) %>%
    dplyr::summarize(plus=sum(tlx_strand=="+"), minus=sum(tlx_strand=="-"), breaks_fc=log2(sum(tlx_strand=="+")/sum(tlx_strand=="-")))

  pdf("reports/strandbias_boxplots.pdf", width=8.27, height=11.69, paper="a4")
  strandbias_df = macs_rdc$qvalues %>%
    dplyr::mutate(qvalue_group=as.numeric(gsub("\\((\\d+).*", "\\1", as.character(cut(qvalue_score, c(0,1,2,10,25, 50, 75, 100)))))+1) %>%
    dplyr::mutate(qvalue_group=tidyr::replace_na(qvalue_group, max(qvalue_group, na.rm=T)+1)) %>%
    # dplyr::filter(qvalue_score<1) %>%
    dplyr::group_by(qvalue_group) %>%
    dplyr::do((function(z){
      zz<<-z
      z %>%
        df2ranges(qvalue_chrom, qvalue_start, qvalue_end) %>%
        GenomicRanges::reduce() %>%
        as.data.frame() %>%
        dplyr::select(qvalue_chrom=seqnames, qvalue_start=start, qvalue_end=end) %>%
        dplyr::mutate(qvalue_width=qvalue_end-qvalue_start, qvalue_group=z$qvalue_group[1])
    })(.)) %>%
    df2ranges(qvalue_chrom, qvalue_start, qvalue_end) %>%
    innerJoinByOverlaps(tlx_rdc_df %>% df2ranges(Rname, Junction, Junction)) %>%
    dplyr::inner_join(baits_df, by="bait_chrom") %>%
    dplyr::mutate(prey_location=ifelse(tlx_bait_start>=Junction, "Up-stream", "Down-stream")) %>%
    dplyr::filter(!tlx_is_bait_junction) %>%
    dplyr::mutate(tlx_condition=dplyr::case_when(
      tlx_is_bait_junction ~ "Bait junction",
      # tlx_is_bait_chrom & tlx_bait_strand=="+"   ~ "Intra-chromosome translocation (telomere-facing primer)",
      # tlx_is_bait_chrom & tlx_bait_strand=="-"   ~ "Intra-chromosome translocation (centromere-facing primer)",
      # tlx_is_bait_chrom & bait_sgRNA_strand=="+"   ~ "Intra-chromosome translocation (telomere-facing PAM)",
      # tlx_is_bait_chrom & bait_sgRNA_strand=="-"   ~ "Intra-chromosome translocation (centromere-facing PAM)",
      tlx_is_bait_chrom                            ~ "Intra-chromosome translocation",
      !tlx_is_bait_chrom   ~ "Inter-chromosome translocation",
      T ~ "Other"
    )) %>%
    dplyr::group_by(tlx_sample, qvalue_group, prey_location, tlx_condition) %>%
    dplyr::summarize(breaks_fc=log2(sum(tlx_strand=="+")/sum(tlx_strand=="-"))) %>%
    dplyr::filter(is.finite(breaks_fc)) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(qvalue_group) %>%
    dplyr::mutate(qvalue_group_name=paste0("<", qvalue_group), qvalue_group_name=factor(qvalue_group_name, unique(qvalue_group_name)))
  ggplot(strandbias_df) +
    geom_boxplot(aes(x=qvalue_group_name, fill=prey_location, y=breaks_fc)) +
    labs(x="Significance, -log10(pvalue)", y="log2(telomeric translocations/centromeric translocations)") +
    ggpubr::theme_pubclean(base_size=16) +
    facet_wrap(~tlx_condition)
  dev.off()
}