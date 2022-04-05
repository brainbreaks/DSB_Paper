setwd("~/Workspace/DSB_Paper")

library(readr)
library(dplyr)
library(ggplot2)
library(cowplot)
library(stringr)
# library(plyranges)
devtools::load_all('~/Workspace/breaktools/')

proportions = function()
{
  params = macs2_params(extsize=1e5, exttype="symmetrical", llocal=1e7, minqvalue=0.01, effective_size=1.87e9, maxgap=2e5, minlen=1e5)

  #
  # Load genes
  #
  genes_df = gtf_read('~/Workspace/genomes/mm10/annotation/mm10.ncbiRefSeq.gtf.gz')
  genes_ranges = genes_df %>% df2ranges(gene_chrom, gene_start, gene_end)


  #
  # Load TLX
  #
  samples_df = tlx_read_samples("~/Workspace/Datasets/HTGTS/samples/All_samples.tsv", "~/Workspace/Datasets/HTGTS/TLX") %>%
    dplyr::filter((grepl("promoter/enhancer", experiment) & alleles==2 | grepl("concentration", experiment) & concentration==0.4)) %>%
    dplyr::mutate(group=paste0(ifelse(control, "DMSO", "APH"), " (", bait_chrom, ")"), treatment=ifelse(control, "DMSO", "APH"), control=F, )

  tlx_df = tlx_read_many(samples_df, threads=30)
  tlx_df = tlx_remove_rand_chromosomes(tlx_df)
  tlx_df = tlx_extract_bait(tlx_df, bait_size=19, bait_region=12e6)
  tlx_df = tlx_mark_dust(tlx_df)
  tlx_df = tlx_df %>%
    dplyr::group_by(tlx_sample) %>%
    dplyr::filter(dplyr::n()>2000) %>%
    dplyr::ungroup()
  libfactors_df = tlx_libfactors(tlx_df, group="group", normalize_within="group", normalize_between="none", normalization_target="smallest")
  tlx_df = tlx_df %>% dplyr::filter(B_Rname==bait_chrom & tlx_is_bait_chrom & !tlx_is_bait_junction)
  tlx_ranges = tlx_df %>% df2ranges(Rname, Junction, Junction)

  #
  # TLX coverage
  #
  tlxcov_df = tlx_df %>%
    tlx_coverage(group="group", extsize=params$extsize, exttype=params$exttype, libfactors_df=libfactors_df, ignore.strand=F)
  tlxcov_ranges = tlxcov_df %>% df2ranges(tlxcov_chrom, tlxcov_start, tlxcov_end)

  #
  # Replication termination site (Repli-seq data)
  #
  replication_df = readr::read_tsv("data/replication_reduced_subsets.tsv")
  replication_ranges = replication_df %>%
    df2ranges(replication_chrom, pmin(replication_start, replication_end), pmax(replication_start, replication_end))

  #
  # Load RDC
  #
  rdc_df = readr::read_tsv("~/Workspace/Datasets/HTGTS/rdc_macs2_mm10.tsv")
  rdc_ranges = rdc_df %>%
    dplyr::mutate(rdc_region_start=rdc_start-1e6, rdc_region_end=rdc_end+1e6) %>%
    df2ranges(rdc_chrom, rdc_region_start, rdc_region_end)
  rdc_filter = paste0(setdiff(rdc_df$rdc_cluster, c("MACS_018|MACS_019|MACS_033|MACS_037|MACS_038|MACS_039|MACS_042|MACS_005|MACS_007")), collapse="|")

  intersection_ranges = as.data.frame(GenomicRanges::intersect(replication_ranges, genes_ranges)) %>%
    dplyr::mutate(intersection_chrom=seqnames, intersection_start=start, intersection_end=end) %>%
    dplyr::mutate(intersection_width=intersection_end-intersection_start) %>%
    GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns=T)
  genes2replication_df = innerJoinManyByOverlaps(list(intersection_ranges, replication_ranges, genes_ranges, tlx_ranges, rdc_ranges)) %>%
    dplyr::group_by(replication_chrom, replication_start, replication_end) %>%
    dplyr::filter(mean(gene_strand==replication_strand) %in% c(0, 1)) %>%
    dplyr::group_by(gene_chrom, gene_start, gene_end) %>%
    dplyr::filter(mean(gene_strand==replication_strand) %in% c(0, 1)) %>%
    dplyr::group_by(intersection_chrom, intersection_start, intersection_end, gene_strand, replication_strand, tlx_strand, intersection_width) %>%
    dplyr::summarise(tlx_count=dplyr::n()) %>%
    dplyr::mutate(direction=ifelse(gene_strand==replication_strand, "Co-directional", "Head-on")) %>%
    dplyr::group_by(intersection_chrom, intersection_start, intersection_end, intersection_width, direction) %>%
    dplyr::summarise(tlx_proportion=ifelse(replication_strand=="-", tlx_count[tlx_strand=="+"]/tlx_count[tlx_strand=="-"], tlx_count[tlx_strand=="-"]/tlx_count[tlx_strand=="+"])) %>%
    dplyr::filter(intersection_width>=1e4)

    # dplyr::group_by(intersection_chrom, intersection_start, intersection_end, gene_strand, replication_strand, tlx_strand, intersection_width) %>%
    # dplyr::summarise(tlx_count=dplyr::n()) %>%
    # dplyr::group_by(intersection_chrom, intersection_start, intersection_end, gene_strand, replication_strand, intersection_width) %>%
    # dplyr::summarise(tlx_proportion=tlx_count[tlx_strand=="+"]/tlx_count[tlx_strand=="-"]) %>%
    # dplyr::filter(intersection_width>=1e4)

  ggplot(genes2replication_df) +
    geom_hline(yintercept=0, linetype="dashed", size=1) +
    geom_boxplot(aes(y=log2(tlx_proportion), fill=direction, x=direction)) +
    # geom_jitter(aes(y=log2(tlx_proportion), group=direction, x=1), width=0.1) +
    geom_point(aes(y=log2(tlx_proportion), group=direction, x=direction), shape=21, size=2, position=position_jitter(width=0.2)) +
    labs(fill="Replication/transcription direction", y="Junction strand proportion, log2(op/co-dir)", x="") +
    scale_fill_manual(values=c("Head-on"="#E28C8E", "Co-directional"="#6E97B2")) +
    theme_classic(base_size=18) +
    guides(fill=guide_legend(nrow=1, byrow=TRUE)) +
    theme(legend.position="bottom", axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())


  ggplot(genes2replication_df) +
    ggridges::stat_density_ridges(aes(x=log2(tlx_proportion), y=gene_strand, fill=replication_strand, height=..density..), geom="density_ridges", bandwidth=0.25, scale=1, from=-3, to=3, alpha=0.4) +
    labs(fill="Replication direction", x="Junction strand proportion, log2(+/-)", y="Transcription direction") +
    scale_fill_manual(values=c("+"="#E31A1C", "-"="#1F78B4"))

  #
  # Junctions strand proportion stratisfied by transcription direction
  #
  pdf("reports/tlx_transcription_prop.pdf", width=8.27, height=11.69)
  x = innerJoinByOverlaps(tlx_ranges, rdc_ranges) %>%
    dplyr::filter(grepl(rdc_filter, rdc_cluster)) %>%
    dplyr::group_by(rdc_cluster_display, rdc_strand, rdc_cluster, tlx_strand, treatment) %>%
    dplyr::summarise(n=dplyr::n()) %>%
    dplyr::group_by(rdc_cluster_display, rdc_strand, rdc_cluster, treatment) %>%
    dplyr::summarise(prop=log2(n[tlx_strand=="+"]/n[tlx_strand=="-"]))
  ggplot(x, aes(x=rdc_strand, y=prop, fill=treatment)) +
    geom_boxplot() +
    geom_jitter() +
    labs(y="Junctions proprortion, log2(+/-)", x="Transcription direction (RDC longest gene)", fill="Treatment")
  dev.off()



  pdf("reports/tlx_replication_prop.pdf", width=8.27, height=11.69)
  x_ranges = innerJoinByOverlaps(tlx_ranges, rdc_ranges) %>% df2ranges(Rname, Junction, Junction)
  x = innerJoinByOverlaps(x_ranges, replication_ranges) %>%
    dplyr::select(-dplyr::matches("_ranges\\.")) %>%
    dplyr::filter(grepl(rdc_filter, rdc_cluster)) %>%
    dplyr::group_by(treatment, rdc_cluster_display, rdc_strand, rdc_cluster, tlx_strand, replication_strand) %>%
    dplyr::summarise(n=dplyr::n()) %>%
    dplyr::group_by(treatment, rdc_cluster_display, rdc_strand, replication_strand, rdc_cluster) %>%
    dplyr::summarise(prop=log2(n[tlx_strand=="+"]/n[tlx_strand=="-"]))
  ggplot(x, aes(x=replication_strand, y=prop, fill=treatment)) +
    geom_boxplot() +
    geom_jitter(aes(color=rdc_strand)) +
    labs(y="Junctions proprortion, log2(+/-)")
  dev.off()

  pdf("reports/tlxcov_proportion_new.pdf", width=11.69, height=8.27)
    x = innerJoinByOverlaps(tlxcov_ranges, rdc_ranges) %>%
      dplyr::filter(grepl(rdc_filter, rdc_cluster)) %>%
      dplyr::group_by(rdc_cluster_display, rdc_strand, rdc_cluster, tlx_strand) %>%
      dplyr::summarise(tlxcov_pileup=max(tlxcov_pileup)) %>%
      dplyr::group_by(rdc_cluster_display, rdc_strand, rdc_cluster) %>%
      dplyr::summarise(prop=log2(tlxcov_pileup[tlx_strand=="+"]/tlxcov_pileup[tlx_strand=="-"]))

    ggplot(x, aes(x=rdc_strand, y=prop)) +
      geom_boxplot() +
      labs(y="Junctions density peaks proportion, log2(+/-)", y="log2( pileup(junctions[+]) / pileup(junctions[-]) )", x="Transcription") +
      geom_jitter() +
      geom_label(aes(label=paste0("log2(", round(2^prop, 2) ,")=", round(prop, 2))), data=x %>% dplyr::group_by(rdc_strand) %>% dplyr::summarise(prop=median(prop)))
  dev.off()


  pdf("reports/tlxcov_replication_prop.pdf", width=11.69, height=8.27)
    x_ranges = innerJoinByOverlaps(tlxcov_ranges, rdc_ranges) %>%
      df2ranges(tlxcov_chrom, tlxcov_start, tlxcov_end)
    x = innerJoinByOverlaps(x_ranges, replication_ranges) %>%
      dplyr::filter(grepl(rdc_filter, rdc_cluster)) %>%
      dplyr::group_by(rdc_cluster_display, rdc_strand, rdc_cluster, replication_strand, tlx_strand) %>%
      dplyr::summarise(tlxcov_pileup=max(tlxcov_pileup)) %>%
      dplyr::group_by(rdc_cluster_display, rdc_strand, rdc_cluster, replication_strand) %>%
      dplyr::summarise(prop=log2(tlxcov_pileup[tlx_strand=="+"]/tlxcov_pileup[tlx_strand=="-"]))

    ggplot(x, aes(x=replication_strand, y=prop)) +
      geom_boxplot() +
      geom_jitter() +
      labs(y="Junctions density peaks proportion, log2(+/-)")
  dev.off()
}