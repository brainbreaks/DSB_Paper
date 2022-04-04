setwd("~/Workspace/Everything")

library(readr)
library(dplyr)
library(ggplot2)
library(cowplot)
library(stringr)
# library(plyranges)
devtools::load_all('~/Workspace/breaktools/')


main = function()
{
  params = macs2_params(extsize=1e5, exttype="symmetrical", llocal=1e7, minqvalue=0.01, effective_size=1.87e9, maxgap=2e5, minlen=1e5)

  #
  # Load RDC
  #
  rdc_df = readr::read_tsv("~/Workspace/Datasets/HTGTS/rdc_macs2_mm10.tsv")
  rdc_ranges = rdc_df %>%
    dplyr::mutate(rdc_region_start=rdc_start-1e6, rdc_region_end=rdc_end+1e6) %>%
    df2ranges(rdc_chrom, rdc_region_start, rdc_region_end)

  #
  # Load genes
  #
  genes_df = gtf_read('~/Workspace/genomes/mm10/annotation/mm10.ncbiRefSeq.gtf.gz') %>%
    dplyr::filter(gene_cluster_i<=3)
  genes_ranges = genes_df %>% df2ranges(gene_chrom, gene_start, gene_end)
  rdc2genes_df = as.data.frame(IRanges::mergeByOverlaps(rdc_ranges, genes_ranges)) %>%
    dplyr::select(-dplyr::matches("_ranges\\."))

  #
  # Load TLX
  #
  samples_df = tlx_read_samples("~/Workspace/Datasets/HTGTS/samples/All_samples.tsv", "~/Workspace/Datasets/HTGTS/TLX") %>%
    dplyr::filter(!control & (grepl("promoter/enhancer", experiment) & alleles==2 | grepl("concentration", experiment) & concentration==0.4)) %>%
    dplyr::mutate(group=paste0("All (", bait_chrom, ")"))

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
  tlxcov_ranges = tlxcov_df %>%
    df2ranges(tlxcov_chrom, tlxcov_start, tlxcov_end)

  #
  # Replication termination site (Repli-seq data)
  #
  replication_df = readr::read_tsv("data/replication_reduced_subsets.tsv")
  replication_ranges = replication_df %>%
    dplyr::mutate(replication_start=pmin(replication_start, replication_end), replication_end=pmax(replication_start, replication_end)) %>%
    df2ranges(replication_chrom, replication_end, replication_end)

  rdc2replication_df = as.data.frame(IRanges::mergeByOverlaps(rdc_ranges, replication_ranges))  %>%
    dplyr::select(-dplyr::matches("_ranges\\."))

  rts_df = replication_df %>%
    dplyr::filter(replication_end>replication_start) %>%
    dplyr::select(rts_chrom=replication_chrom, rts_pos=replication_end)
  rts_ranges = rts_df %>% df2ranges(rts_chrom, rts_pos, rts_pos)

  #
  # Replication termination sites
  #
  rdc2rts_df = as.data.frame(IRanges::mergeByOverlaps(rdc_ranges, rts_ranges)) %>%
    dplyr::select(dplyr::matches("rdc"), rts_pos)
  rdc2rts_ranges = rdc2rts_df %>%
    dplyr::mutate(rts_start=rts_pos-1e5, rts_end=rts_pos+1e5) %>%
    df2ranges(rdc_chrom, rts_start, rts_end)
  rdc2rts_df = as.data.frame(IRanges::mergeByOverlaps(tlx_ranges, rdc2rts_ranges)) %>%
    dplyr::group_by(rdc_chrom, rdc_cluster, rdc_cluster_display, rdc_strand, rts_pos) %>%
    dplyr::summarise(rts_junctions_count=dplyr::n()) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(dplyr::desc(rts_junctions_count)) %>%
    dplyr::distinct(rdc_cluster, .keep_all=T)
  rdc2rts_ranges = rdc2rts_df %>%
    dplyr::mutate(rts_start=rts_pos-1e6, rts_end=rts_pos+1e6) %>%
    df2ranges(rdc_chrom, rts_start, rts_end)

  #
  # Plot
  #
  rts2rel_tlxcov_df = as.data.frame(IRanges::mergeByOverlaps(tlxcov_ranges, rdc2rts_ranges)) %>%
    dplyr::select(-dplyr::matches("_ranges\\.")) %>%
    dplyr::mutate(tlxcov_rel_start=tlxcov_start-rts_pos, tlxcov_rel_end=tlxcov_end-rts_pos) %>%
    dplyr::filter(abs(tlxcov_rel_start)<=1e6 & abs(tlxcov_rel_end)<=1e6) %>%
    dplyr::mutate(rdc_cluster_i=as.numeric(as.factor(rdc_cluster))) %>%
    dplyr::group_by(rdc_cluster) %>%
    dplyr::mutate(tlxcov_rel_pileup=tlxcov_pileup/max(tlxcov_pileup, na.rm=T)) %>%
    dplyr::ungroup()


  #
  # Select big clusters
  #
  rdc_filter = paste0(setdiff(rdc_df$rdc_cluster, c("MACS_018|MACS_019|MACS_033|MACS_037|MACS_038|MACS_039|MACS_042|MACS_005|MACS_007")), collapse="|")

  rdc2replication_ggplot = rdc2replication_df %>%
    dplyr::filter(grepl(rdc_filter, rdc_cluster)) %>%
    dplyr::inner_join(rdc2rts_df %>% dplyr::select(rdc_cluster, rts_pos), by="rdc_cluster") %>%
    dplyr::filter(pmax(replication_start, replication_end)>=rts_pos-1e6 & pmin(replication_start, replication_end)<=rts_pos+1e6) %>%
    dplyr::mutate(replication_rel_start=replication_start-rts_pos, replication_rel_end=replication_end-rts_pos)

  rdc2genes_ggplot = rdc2genes_df %>%
    dplyr::filter(grepl(rdc_filter, rdc_cluster)) %>%
    dplyr::inner_join(rdc2rts_df %>% dplyr::select(rdc_cluster, rts_pos), by="rdc_cluster") %>%
    dplyr::filter(gene_start>=rts_pos-1e6 & gene_end<=rts_pos+1e6) %>%
    dplyr::mutate(gene_rel_start=gene_start-rts_pos, gene_rel_end=gene_end-rts_pos) %>%
    dplyr::mutate(segment_start=ifelse(gene_strand=="+", gene_rel_start, gene_rel_end), segment_end=ifelse(gene_strand=="+", gene_rel_end, gene_rel_start))

  rts2rel_tlxcov_ggplot = rts2rel_tlxcov_df %>%
    dplyr::group_by(rdc_cluster) %>%
    dplyr::mutate(tlxcov_start=tlxcov_rel_start, tlxcov_end=tlxcov_rel_end, tlxcov_pileup=tlxcov_pileup/max(tlxcov_pileup)) %>%
    dplyr::ungroup()

  plist_strands = c("+", "-")
  plist = lapply(1:length(plist_strands), FUN=function(s) {
    rts2rel_tlxcov_ggplot.s = rts2rel_tlxcov_ggplot %>% dplyr::filter(rdc_strand==plist_strands[s])
    rdc2genes_ggplot.s = rdc2genes_ggplot %>% dplyr::filter(rdc_strand==plist_strands[s]) %>% dplyr::arrange(dplyr::desc(gene_length))
    rdc2replication_ggplot.s = rdc2replication_ggplot %>% dplyr::filter(rdc_strand==plist_strands[s])

    g = ggplot() +
      geom_vline(xintercept=0, size=0.3, alpha=0.1) +
      geom_tlxcov(rts2rel_tlxcov_ggplot.s) +
      geom_segment(aes(x=segment_start, xend=segment_end, y=-0.1*gene_cluster_i, yend=-0.1*gene_cluster_i), size=0.3, arrow=grid::arrow(length=unit(4,"pt")), data=rdc2genes_ggplot.s) +
      geom_segment(aes(x=segment_start, xend=segment_start, y=-0.1*gene_cluster_i-0.05, yend=-0.1*gene_cluster_i+0.05), size=0.3, data=rdc2genes_ggplot.s) +
      geom_text(aes(x=(gene_rel_end+gene_rel_start)/2, y=-0.1*(gene_cluster_i+1)+0.05, label=gene_id), data=rdc2genes_ggplot.s %>% dplyr::distinct(rdc_cluster, .keep_all=T), size=3, vjust=1, hjust=0, color="#FF6500") +
      geom_segment(aes(x=replication_rel_start, xend=replication_rel_end, y=-0.5+y, yend=-0.5+y, color=replication_strand), size=0.1, arrow=grid::arrow(length=unit(4,"pt")), data=rdc2replication_ggplot.s %>% dplyr::mutate(y=ifelse(replication_strand=="+", 0.02, -0.02))) +
      # geom_vline(aes(xintercept=replication_start), size=0.1, alpha=0.1, data=rdc2replication_ggplot.s) +
      # geom_vline(aes(xintercept=replication_end), size=0.1, alpha=0.1, data=rdc2replication_ggplot.s)
      facet_grid(paste(rdc_chrom, gsub("[^0-9]+0?", "", rdc_cluster))~., scales="free_y") +
      labs(x="Relative distance from RTS (100kb)", y="", title=plist_strands[s]) +
      scale_x_continuous(labels=scales::label_number(accuracy=1, scale=1e-5, suffix="")) +
      coord_cartesian(xlim=c(-1e6, 1e6)) +
      theme_classic(base_size=16) +
      theme(strip.background=element_blank(), strip.text.y=element_text(size=10), legend.position="bottom", axis.title.y=ggplot2::element_blank())
    if(s<length(plist_strands)) {
      g = g + theme(axis.title.x=ggplot2::element_blank())
    }
    g
  })


  pdf("reports/tlxcov_relative_new.pdf", width=8.27, height=11.69)
  plist
  dev.off()
}