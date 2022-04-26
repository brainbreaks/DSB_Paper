setwd("~/Workspace/DSB_paper")

library(readr)
library(dplyr)
library(ggplot2)
library(cowplot)
library(stringr)
# library(plyranges)
devtools::load_all("~/Workspace/breaktools/")

ggplot_rdc_breaks = function(rdc2rdc_df, rdc2tlxcov_df, rdc2tlx_df, rdc2genes_df, rdc2repliseq_df, rdc2replication_df=NULL, rdc2offtargets_df=NULL, rdc_filter=".*")
{
  rdc2rdc_dff = rdc2rdc_df %>% dplyr::filter(grepl(rdc_filter, rdc_cluster))
  rdc2tlxcov_dff = rdc2tlxcov_df %>% dplyr::filter(grepl(rdc_filter, rdc_cluster))
  rdc2tlx_dff = rdc2tlx_df %>% dplyr::filter(grepl(rdc_filter, rdc_cluster))
  rdc2genes_dff = rdc2genes_df %>% dplyr::filter(grepl(rdc_filter, rdc_cluster))
  rdc2repliseq_dff = rdc2repliseq_df %>% dplyr::filter(grepl(rdc_filter, rdc_cluster))

  rdc2genes_dff = rdc2genes_dff %>% dplyr::mutate(segment_start=ifelse(gene_strand=="+", gene_start, gene_end), segment_end=ifelse(gene_strand=="+", gene_end, gene_start))
  rdc2genes_dff.long = rdc2genes_dff %>% dplyr::filter(gene_length>=2e5)
  rdc2genes_dff.short = rdc2genes_dff %>% dplyr::filter(gene_length<2e5)
  g = ggplot()

  if(!is.null(rdc2replication_df)) {
    rdc2replication_dff = rdc2replication_df %>% dplyr::filter(grepl(rdc_filter, rdc_cluster))
    if(nrow(rdc2replication_dff)>0) {
      # pdf("reports/x.pdf", width=6*8.27, height=6*11.69, paper="a4")
      g = g +
        geom_segment(aes(x=replication_start, xend=replication_end, y=-19+y, yend=-19+y, color=replication_strand), size=0.1, arrow=grid::arrow(length=unit(1,"pt")), data=rdc2replication_dff %>% dplyr::mutate(y=ifelse(replication_strand=="+", 0.2, -0.2))) +
        geom_vline(aes(xintercept=replication_start), size=0.1, alpha=0.1, data=rdc2replication_dff) +
        geom_vline(aes(xintercept=replication_end), size=0.1, alpha=0.1, data=rdc2replication_dff)
      # dev.off()
    }
  }



  g = g +
    geom_rect(aes(ymin=-Inf, ymax=Inf, xmin=rdc_start, xmax=rdc_end), alpha=0.3, fill="#CCCCCC", data=rdc2rdc_dff, size=0.1)  +
    geom_tlxcov(rdc2tlxcov_dff, scale=15) +
    geom_tile(aes(x=repliseq_start/2+repliseq_end/2, y=repliseq_fraction-19), data=rdc2repliseq_dff, fill=scales::colour_ramp(c("#FFFFFF00", "#FF0000FF"))(rdc2repliseq_dff$repliseq_value)) +
    # geom_rect(aes(xmin=tlxcov_start, xmax=tlxcov_end-1, ymin=19, ymax=19+tlxcov_pileup*15, group=tlx_strand), fill=ifelse(rdc2tlxcov_dff$tlx_strand=="+", "#FF0000", "#0000FF"), alpha=0.5, data=rdc2tlxcov_area_df) +
    geom_step(aes(x=tlxcov_start, y=tlxcov_pileup*15, group=tlx_strand), data=rdc2tlxcov_dff, size=0.1, alpha=0.5, color=ifelse(rdc2tlxcov_dff$tlx_strand=="+", "#E31A1C", "#1F78B4")) +
    geom_segment(aes(x=Junction, xend=Junction, y=-1.2, yend=tlx_strand_sign-1.2, color=tlx_strand), data=rdc2tlx_dff %>% dplyr::mutate(tlx_strand_sign=ifelse(tlx_strand=="+", 1, -1)), size=0.1, alpha=0.6)

  if(!is.null(rdc2offtargets_df)) {
    rdc2offtargets_dff = rdc2offtargets_df %>% dplyr::filter(grepl(rdc_filter, rdc_cluster))
    if(nrow(rdc2offtargets_dff)>0) {
      g = g +
        geom_rect(aes(ymin=-25.5, ymax=-20.5, xmin=offtarget_start-1e5, xmax=offtarget_end+1e5), fill="#666666", data=rdc2offtargets_dff, size=0.1) +
        geom_rect(aes(ymin=-Inf, ymax=Inf,    xmin=offtarget_start,      xmax=offtarget_end), alpha=0.3, fill="#FF0000", data=rdc2offtargets_dff, size=0.1)
    }
  }


  g = g +
    geom_segment(aes(x=segment_start, xend=segment_end, y=-19-gene_cluster_i*4-0.5, yend=-19-gene_cluster_i*4-0.5), size=0.3, arrow=grid::arrow(length=unit(0.8,"pt")), color="#CCCCCC", data=rdc2genes_dff.short) +
    geom_segment(aes(x=segment_start, xend=segment_end, y=-19-gene_cluster_i*4-0.5, yend=-19-gene_cluster_i*4-0.5), size=0.3, arrow=grid::arrow(length=unit(0.8,"pt")), color=ifelse(rdc2genes_dff.long$gene_strand=="+", "#a14f4f", "#4f87a1"), data=rdc2genes_dff.long) +
    geom_text(aes(x=(gene_end+gene_start)/2, y=-19-gene_cluster_i*4-2, label=gene_id), data=rdc2genes_dff.long, size=1, color="#FF6500") +
    geom_hline(yintercept=-0.2, size=0.1) +
    facet_wrap(~rdc_cluster_display, scales="free_x") +
    guides(fill=guide_legend(nrow=1, byrow=TRUE)) +
    ggpubr::theme_pubclean(base_size=5) +
    theme(legend.position="bottom", axis.text.x=element_text(angle=45, hjust=1)) +
    scale_fill_manual(values=c("-"="#A6CEE3", "+"="#FB9A99", "cluster"="#00FF00")) +
    scale_color_manual(values=c("-"="#E31A1C", "+"="#1F78B4", "cluster"="#00FF00")) +
    scale_y_continuous(breaks=c(-3, seq(1, 16, 5), 24.5, 32)-19, labels=c("gene", seq(1, 16, 5), 0.5, 1)) +
    scale_x_continuous(labels=scales::label_number(accuracy=0.1, scale=1e-6, suffix="Mb")) +
    labs(x="", y="")
  g
}



repliseqClustering = function()
{
  params = macs2_params(extsize=1e5, exttype="symmetrical", llocal=1e7, minqvalue=0.01, effective_size=1.87e9, maxgap=2e5, minlen=1e5)

  #
  # Load gene annotations
  #
  genes_df = gtf_read('~/Workspace/genomes/mm10/annotation/mm10.ncbiRefSeq.gtf.gz')
  genes_ranges = genes_df %>% df2ranges(gene_chrom, gene_start, gene_end)

  #
  # Load TLX
  #
  samples_df = tlx_read_samples("~/Workspace/Datasets/HTGTS/samples/All_samples.tsv", "~/Workspace/Datasets/HTGTS") %>%
    dplyr::filter(!control & (grepl("promoter/enhancer", experiment) & alleles==2 | grepl("concentration", experiment) & concentration==0.4)) %>%
    dplyr::mutate(group=paste0("All (", bait_chrom, ")"))

  tlx_df = tlx_read_many(samples_df, threads=30)
  tlx_df = tlx_remove_rand_chromosomes(tlx_df)
  tlx_df = tlx_extract_bait(tlx_df, bait_size=19, bait_region=2e6)
  tlx_df = tlx_mark_dust(tlx_df)
  libfactors_df = tlx_libfactors(tlx_df, group="group", normalize_within="group", normalize_between="none", normalization_target="smallest")

  #
  # Search for offtargets
  #
  baits_df = readr::read_tsv("~/Workspace/Datasets/HTGTS/wei_pnas2018_baits.tsv")
  offtarget_df = readr::read_tsv("~/Workspace/Datasets/HTGTS/offtargets_pnas_mm10.tsv")
  offtarget2bait_df = offtarget_df %>%
    dplyr::filter(offtarget_is_primary==0) %>%
    dplyr::select(bait_chrom, offtarget_chrom, offtarget_start, offtarget_end, offtarget_strand, offtarget_sequence) %>%
    dplyr::inner_join(baits_df %>% dplyr::select(dplyr::matches("bait_")), by="bait_chrom")
  tlx_df = tlx_mark_offtargets(tlx_df, offtarget2bait_df=offtarget2bait_df, offtarget_region=2e4)

  #
  # Calculate TLX coverage
  #
  tlxcov_df = tlx_df %>%
    dplyr::filter(B_Rname==bait_chrom & tlx_is_bait_chrom & !tlx_is_bait_junction) %>%
    dplyr::group_by(tlx_sample) %>%
    dplyr::filter(dplyr::n() >= 2000) %>%
    dplyr::ungroup() %>%
    dplyr::filter(!tlx_control) %>%
    dplyr::mutate(tlx_group="APH") %>%
    dplyr::group_by(Rname) %>%
    dplyr::mutate(dbscan_cluster=dbscan::dbscan(matrix(Junction), minPts=20, eps=100)$cluster) %>%
    dplyr::mutate(tlx_strand=ifelse(dbscan_cluster==0, tlx_strand, "cluster")) %>%
    # dplyr::filter(Rname==tlx_group) %>%
    tlx_coverage(group="group", exttype=params$exttype, extsize=params$extsize, libfactors_df=libfactors_df, ignore.strand=F)

  #
  # Write TLX coverage
  #
  if(F) {
    tlxcov_write_bedgraph(tlxcov_df, path=paste0("reports/bedgaph-", extsize), group="group", ignore.strand=T)
    tlxcov_write_bedgraph(tlxcov_df, path=paste0("reports/bedgaph-", extsize), group="group", ignore.strand=F)
  }


  #
  # Load RDC
  #
  rdc_df = readr::read_tsv("~/Workspace/Datasets/HTGTS/rdc_macs2_mm10.tsv") %>%
    dplyr::mutate(rdc_region_start=rdc_start-2e6, rdc_region_end=rdc_end+2e6)

  #
  # Expand RDC displayed range to include other neighbouring RDC
  #

  tlxcov_ranges = tlxcov_df %>% df2ranges(tlxcov_chrom, tlxcov_start, tlxcov_end)
  rdc_wide_ranges = rdc_df %>% dplyr::mutate(rdc_original_start=rdc_start, rdc_original_end=rdc_end) %>% df2ranges(rdc_chrom, rdc_region_start, rdc_region_end)
  rdc_narrow_ranges = rdc_df %>% dplyr::mutate(rdc_original_start=rdc_start, rdc_original_end=rdc_end) %>% df2ranges(rdc_chrom, rdc_start, rdc_end)
  rdc2rdc_df = as.data.frame(IRanges::mergeByOverlaps(rdc_wide_ranges, rdc_narrow_ranges)) %>%
    dplyr::select(-dplyr::matches("_ranges\\."))
  rdc_df = rdc2rdc_df %>%
    dplyr::group_by(rdc_cluster, rdc_cluster_display) %>%
    dplyr::mutate(rdc_region_start=min(c(rdc_start, rdc_region_start)), rdc_region_end=max(c(rdc_region_end, rdc_end))) %>%
    dplyr::ungroup() %>%
    dplyr::distinct(rdc_cluster, rdc_cluster_display, .keep_all=T) %>%
    dplyr::mutate(rdc_start=rdc_original_start, rdc_end=rdc_original_end) %>%
    dplyr::select(-rdc_original_start, -rdc_original_end)
  rdc_ranges = rdc_df %>% df2ranges(rdc_chrom, rdc_region_start, rdc_region_end)

  rdc2tlxcov_df = as.data.frame(IRanges::mergeByOverlaps(rdc_ranges, tlxcov_ranges)) %>%
    dplyr::select(-dplyr::matches("_ranges\\.")) %>%
    dplyr::filter(!tlx_control & tlxcov_start>=rdc_region_start & tlxcov_end<=rdc_region_end) %>%
    dplyr::group_by(rdc_cluster, rdc_cluster_display) %>%
    dplyr::mutate(tlxcov_pileup=tlxcov_pileup/max(tlxcov_pileup[tlx_strand %in% c("+", "-")])) %>%
    dplyr::mutate(tlxcov_pileup=pmin(tlxcov_pileup, 1)) %>%
    dplyr::ungroup()

  #
  # Search for offtargets
  #
  offtargets_ranges = offtarget_df %>% df2ranges(offtarget_chrom, offtarget_start, offtarget_end, offtarget_strand)
  rdc2offtargets_df = as.data.frame(IRanges::mergeByOverlaps(rdc_ranges, offtargets_ranges)) %>%
    dplyr::select(-dplyr::matches("_ranges\\."))

  rdc2genes_df = as.data.frame(IRanges::mergeByOverlaps(rdc_ranges, genes_ranges)) %>%
    dplyr::select(-dplyr::matches("_ranges\\.")) %>%
    dplyr::filter(gene_cluster_i<=3) %>%
    dplyr::mutate(gene_start=pmax(rdc_region_start, gene_start), gene_end=pmin(rdc_region_end, gene_end))
  tlx_ranges = tlx_df %>% df2ranges(Rname, Junction, Junction, tlx_strand)

  rdc2tlx_df = as.data.frame(IRanges::mergeByOverlaps(rdc_ranges, tlx_ranges)) %>%
    dplyr::select(-dplyr::matches("_ranges\\.")) %>%
    dplyr::filter(!tlx_control)

  # Load Repli-SEQ data
  replication_df = readr::read_tsv("data/replication_reduced_subsets.tsv")
  replication_ranges = replication_df %>%
    dplyr::mutate(seqnames=replication_chrom, start=pmin(replication_start, replication_end), end=pmax(replication_start, replication_end)) %>%
    GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns=T)
  rdc2replication_df = as.data.frame(IRanges::mergeByOverlaps(rdc_ranges, replication_ranges))  %>%
    dplyr::select(-dplyr::matches("_ranges\\."))

  repliseq_df = readr::read_tsv("~/Workspace/Datasets/zhao_bmc_repliseq_2020/preprocessed/repliseq_NPC.tsv")
  repliseq_ranges = repliseq_df %>% df2ranges(repliseq_chrom, repliseq_start, repliseq_end)
  rdc2repliseq_df = as.data.frame(IRanges::mergeByOverlaps(rdc_ranges, repliseq_ranges))  %>%
    dplyr::select(-dplyr::matches("_ranges\\.")) %>%
    dplyr::group_by(rdc_cluster, rdc_cluster_display, repliseq_start) %>%
    dplyr::mutate(repliseq_value=repliseq_value-min(repliseq_value)+0.01) %>%
    dplyr::mutate(repliseq_value=repliseq_value^1.5) %>%
    dplyr::mutate(repliseq_value=repliseq_value/quantile(repliseq_value, 0.95, na.rm=T), repliseq_value=pmin(repliseq_value, 1)) %>%
    dplyr::ungroup()

  #
  # Plot specific RDC examples
  #
  rdc_samples_df = rdc2tlxcov_df %>%
    dplyr::arrange(as.numeric(gsub(".*([0-9]).*", "\\1", rdc_cluster))) %>%
    # dplyr::filter(grepl("MACS_002", rdc_cluster)) %>%
    dplyr::group_by(rdc_cluster, rdc_cluster_display) %>%
    dplyr::summarize(rdc_cluster_length=max(tlxcov_end)-min(tlxcov_start)) %>%
    dplyr::ungroup() %>%
    dplyr::left_join(readr::read_tsv("data/rdc_subsets.tsv"), by="rdc_cluster") %>%
    dplyr::mutate(rdc_repliseq_type=ifelse(is.na(rdc_repliseq_type), paste0("Other", cut(1:dplyr::n(), 5)), rdc_repliseq_type)) %>%
    dplyr::select(rdc_cluster, rdc_cluster_display, rdc_repliseq_type, rdc_rrs, rdc_cluster_length) %>%
    dplyr::filter(!is.na(rdc_repliseq_type)) %>%
    dplyr::group_by(rdc_repliseq_type) %>%
    dplyr::summarize(rdc_filter=paste(rdc_cluster, collapse="|"), rdc_examples_n=dplyr::n(), rdc_examples_length=sum(rdc_cluster_length)) %>%
    dplyr::mutate(rdc_examples_missing=max(rdc_examples_n) - rdc_examples_n, rdc_examples_mislength=max(rdc_examples_length) - rdc_examples_length)

  plist = lapply(split(rdc_samples_df, f=rdc_samples_df$rdc_repliseq_type), FUN=function(df) {
    dff<<-df
    # df = data.frame(rdc_filter="MACS_014", rdc_examples_mislength=0, rdc_examples_length=100)
    # df = rdc_samples_df %>% dplyr::slice(1)  %>% dplyr::mutate(rdc_filter="MACS_001")
    ggplot_rdc_breaks(rdc2rdc_df=rdc2rdc_df, rdc2tlxcov_df=rdc2tlxcov_df, rdc2tlx_df=rdc2tlx_df, rdc2genes_df=rdc2genes_df, rdc2repliseq_df=rdc2repliseq_df, rdc2offtargets_df=rdc2offtargets_df, rdc2replication_df=rdc2replication_df, rdc_filter=df$rdc_filter) +
      ggtitle(df$rdc_repliseq_type) +
      geom_segment(aes(x=0, xend=rdc_examples_mislength, y=-0.2, yend=-0.2), data=df %>% dplyr::mutate(rdc_cluster_display="X"), size=1, color="#FFFFFF") +
      facet_grid(.~rdc_cluster_display, scales="free_x", space="free_x") +
      scale_x_continuous(labels=scales::label_number(accuracy=0.1, scale=1e-6, suffix="Mb")) +
      ggpubr::theme_pubclean(base_size=5) +
      theme(axis.text.x=element_text(angle=45, hjust=1))
  })
  pdf("reports/tlx2repliseq_samples_baitchrom_macs6.pdf", width=5*8.27, height=1*11.6)
  g_legend = get_legend(plist[[1]] + theme(legend.box.margin=margin(0, 0, 0, 12)))
  g_plots = cowplot::plot_grid(plotlist=lapply(plist, function(p) p+theme(legend.position="none")), align="v", axis="t", ncol=1)
  cowplot::plot_grid(g_plots, g_legend, ncol=1, rel_heights=c(19, 1))
  dev.off()

  #
  # pdf("reports/tlx2repliseq_samples_baitchrom_macs6.pdf", width=4, height=4)
  # cowplot::plot_grid(g_legend, ncol=1)
  # dev.off()

  #
  # Calculate breaks in relationship to replication direction
  #

  rdc_ranges = rdc_df %>% df2ranges(rdc_chrom, rdc_region_start, rdc_region_end)
  tlx_ranges = GenomicRanges::makeGRangesFromDataFrame(tlx_df %>% dplyr::mutate(seqnames=Rname, start=Junction, end=Junction, strand=tlx_strand) %>% dplyr::select(-Strand), ignore.strand=F, keep.extra.columns=T)
  rdc2tlx_df = as.data.frame(IRanges::mergeByOverlaps(rdc_ranges, tlx_ranges)) %>%
    dplyr::select(-dplyr::matches("_ranges\\.")) %>%
    dplyr::filter(!tlx_control)
  rdc2tlx_ranges = rdc2tlx_df %>% df2ranges(Rname, Junction, Junction)

  replication_end_ranges = replication_df %>%
    dplyr::mutate(replication_start=ifelse(replication_strand=="+", replication_end-2e5, replication_start), replication_end=ifelse(replication_strand=="+", replication_end, replication_start+2e5)) %>%
    dplyr::mutate(start=pmin(replication_start, replication_end), end=pmax(replication_start, replication_end)) %>%
    df2ranges(replication_chrom, start, end)

  rdc2replication_breakscount_df = as.data.frame(IRanges::mergeByOverlaps(rdc2tlx_ranges, replication_end_ranges)) %>%
    dplyr::select(-dplyr::matches("_ranges\\.")) %>%
    dplyr::group_by(rdc_cluster, replication_chrom, replication_start, replication_end, tlx_strand, replication_strand) %>%
    dplyr::summarize(breaks_count=dplyr::n()) %>%
    dplyr::group_by(rdc_cluster, replication_chrom, replication_start, replication_end, replication_strand) %>%
    dplyr::mutate(breaks_prop=breaks_count/sum(breaks_count)) %>%
    dplyr::ungroup()

  rdc2replication_breakscount_largest_df = rdc2replication_breakscount_df %>%
    dplyr::arrange(dplyr::desc(breaks_prop)) %>%
    dplyr::group_by(rdc_cluster, tlx_strand) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup() %>%
    dplyr::select(rdc_cluster, replication_chrom, replication_start, replication_end, tlx_strand)


  rdc2replication_breakscount_largest_ranges = rdc2replication_breakscount_largest_df %>% df2ranges(replication_chrom, replication_start, replication_end)
  rdc2replication_breakscount_ggplot = as.data.frame(IRanges::mergeByOverlaps(rdc2replication_breakscount_largest_ranges, genes_ranges)) %>%
    dplyr::select(-dplyr::matches("_ranges\\.")) %>%
    dplyr::group_by(rdc_cluster, replication_chrom, replication_start, replication_end, tlx_strand) %>%
    dplyr::summarize(transcription_strand=paste(unique(gene_strand), collapse="")) %>%
    dplyr::ungroup() %>%
    dplyr::select(-tlx_strand, -replication_chrom) %>%
    dplyr::left_join(rdc2replication_breakscount_df, by=c("rdc_cluster", "replication_start", "replication_end")) %>%
    dplyr::mutate(breaks_count=tidyr::replace_na(breaks_count, 0))
    # dplyr::filter(grepl("RDC_083|RDC_093|RDC_103|RDC_104|RDC_101|RDC_001|RDC_032|RDC_067|RDC_089|RDC_090|RDC_099|RDC_028|RDC_038|RDC_050|RDC_054|RDC_056|RDC_076|RDC_081|RDC_020|RDC_033|RDC_049", rdc_cluster)) %>%
    # dplyr::filter(grepl("RDC_086|RDC_079|RDC_088|RDC_089|RDC_090|RDC_093|RDC_080|RDC_081|RDC_099|RDC_101|RDC_103|RDC_104", rdc_cluster))


  # pdf("reports/breaks2replication_direction.pdf", width=8.27, height=11.6)
  g = ggplot(rdc2replication_breakscount_ggplot) +
    geom_boxplot(aes(y=breaks_prop, x=replication_strand, fill=tlx_strand), outlier.shape=NA) +
    ggiraph::geom_point_interactive(aes(y=breaks_prop, x=replication_strand, fill=tlx_strand, tooltip=paste0(rdc_cluster)), pch = 21, position= position_jitterdodge()) +
    theme_bw(base_size = 9) +
    facet_wrap(~transcription_strand)
  g
  ggiraph::girafe(ggobj = g)
# dev.off()


  #
  # Shift
  #
  tlx_df = tlx_read_many(samples_df, threads=30)
  tlx_df = tlx_remove_rand_chromosomes(tlx_df)
  tlx_df = tlx_extract_bait(tlx_df, bait_size=19, bait_region=2e6)
  tlx_df = tlx_df %>% dplyr::inner_join(samples_df, by=c("tlx_sample"="sample"))
  # tlx_df = tlx_df %>% dplyr::filter(!tlx_is_bait_junction & tlx_is_bait_chromosome)

  rdc_filter = rdc_samples_df %>%
    dplyr::filter(grepl("External Ori|Internal Termination", rdc_repliseq_type)) %>%
    dplyr::summarize(rdc_filter=paste(rdc_filter, collapse="|")) %>%
    .$rdc_filter
  rdc_filter = "RDC_019|RDC_083|RDC_093|RDC_103|RDC_101|RDC_001|RDC_010|RDC_032|RDC_088|RDC_089|RDC_0|RDC_090|RDC_099|RDC_007|RDC_034|RDC_038|RDC_050|RDC_054|RDC_056|RDC_059|RDC_076|RDC_081|RDC_033|RDC_053|RDC_008|RDC_072|RDC_086|RDC_094|RDC_097|RDC_062|RDC_073|RDC_075|RDC_079|RDC_100"
  # rdc_filter = "RDC_088|RDC_089|RDC_090"
  rdc_filter = ".*"

  rdc_ranges = rdc_df %>% dplyr::filter(grepl(rdc_filter, rdc_cluster)) %>% df2ranges(rdc_chrom, rdc_start, rdc_end)
  tlx_ranges = tlx_df %>% df2ranges(Rname, Junction, Junction, tlx_strand)
  rdc2tlx_narrow_df = as.data.frame(IRanges::mergeByOverlaps(rdc_ranges, tlx_ranges)) %>%
    dplyr::select(-dplyr::matches("_ranges\\.")) %>%
    dplyr::filter(!tlx_control) %>%
    # dplyr::filter(grepl("chr8|chr17|chr6", tlx_bait_chrom) & Rname=="chr6") %>%
    # dplyr::filter(!grepl("SRR", tlx_sample)) %>%
    dplyr::group_by(tlx_sample, rdc_cluster, Rname) %>%
    dplyr::filter(dplyr::n()>=30) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(junction_position=ifelse(Junction>tlx_bait_start, "after bait", "before bait")) %>%
    dplyr::mutate(JunctionRel=(Junction-rdc_start)/(rdc_end-rdc_start))

  rdc2tlx_narrow_longdf = rdc2tlx_narrow_df %>%
    dplyr::group_by(tlx_sample, rdc_cluster, tlx_strand, tlx_bait_chrom, tlx_bait_strand, Rname, tlx_is_bait_chrom, junction_position) %>%
    dplyr::summarize(strand_breaks_n=dplyr::n()) %>%
    dplyr::mutate(tlx_strand_name=ifelse(tlx_strand=="+", "plus", "minus"))  %>%
    dplyr::mutate(strand_breaks_n=tidyr::replace_na(strand_breaks_n,0)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(rdc_cluster_number=as.numeric(factor(rdc_cluster)))
  rdc2tlx_narrow_widedf = rdc2tlx_narrow_longdf %>%
    reshape2::dcast(tlx_sample+rdc_cluster+tlx_bait_chrom+tlx_bait_strand+rdc_cluster_number+Rname+tlx_is_bait_chrom+junction_position ~ tlx_strand_name, value.var="strand_breaks_n")  %>%
    dplyr::mutate(minus=tidyr::replace_na(minus,0), plus=tidyr::replace_na(plus,0)) %>%
    dplyr::mutate(breaks_difference=plus-minus, breaks_prop=plus/minus) %>%
    dplyr::group_by(tlx_bait_chrom, tlx_bait_strand, tlx_is_bait_chrom) %>%
    dplyr::mutate(observations_n=length(unique(tlx_sample)), rdc_n=length(unique(rdc_cluster)), samples_n=length(unique(tlx_sample)), tlx_is_bait_chrom=ifelse(tlx_is_bait_chrom, "C-bait", "C-other")) %>%
    dplyr::ungroup()


  ggplot(rdc2tlx_narrow_widedf) +
    geom_jitter(aes(x=paste(rdc_cluster, junction_position), y=breaks_prop, color=tlx_bait_chrom, size=tlx_is_bait_chrom), width=0.2) +
    # geom_line(aes(x=rdc_cluster_number, y=breaks_prop, color=tlx_bait_chrom)) +
    labs(y="breaks(+)/breaks(-)", x="", title=paste0("Chr8+ bait junction with 3 clusters")) +
    scale_size_manual(values=c("C-bait"=5, "C-other"=2.5)) +
    facet_grid(Rname~., scales="free", space="free") +
    coord_flip()


  rdc2tlx_narrow_widedf %>%
    # dplyr::group_by(tlx_is_bait_chrom, tlx_bait_strand) %>%
    # dplyr::summarize(rdc_n=length(unique(rdc_cluster)), breaks_prop=mean(breaks_prop)) %>%
    dplyr::arrange(breaks_prop) %>%
    ggplot() +
      geom_boxplot(aes(x=paste0(tlx_is_bait_chrom, " (obs:", observations_n, " smpl:", samples_n, " rdc:", rdc_n, ")"), y=breaks_prop, fill="+"), width=0.2, size=0.1) +
      facet_grid(tlx_bait_chrom+tlx_bait_strand~., scales="free") +
      labs(y="breaks(+)/breaks(-)", x="", title=paste0("Chr8+ bait junction with 3 clusters")) +
      coord_flip()


  rdc2tlx_narrow_aov = rstatix::anova_test(data=rdc2tlx_narrow_ttest_widedf, dv=breaks_prop, wid=tlx_sample, within=rdc_cluster)
  rdc2tlx_narrow_aov_table = rstatix::get_anova_table(rdc2tlx_narrow_aov)


  pdf("reports/breaks_strand_balance.pdf", width=8.27, height=11.69, paper="a4")
  ggplot(rdc2tlx_narrow_widedf) +
    geom_boxplot(aes(x=-1, y=minus, fill=rdc_cluster), width=0.2) +
    geom_boxplot(aes(x=1, y=plus, fill=rdc_cluster), width=0.2) +
    geom_text(aes(x=-1, y=minus, label=tlx_sample), hjust=0.9, color="#000000") +
    geom_text(aes(x=1, y=plus, label=tlx_sample), hjust=0.2, color="#000000") +
    geom_segment(aes(x=-1, xend=1, y=minus, yend=plus, color=rdc_cluster)) +
    labs(y="Breaks count", x="Junction strand", title=paste0("Chr8+ bait junction with 3 clusters")) +
    scale_x_continuous(breaks=c(-1, 1), labels=c("-", "+")) +
    guides(fill="none") +
    theme_classic(base_size=16)


  ggplot(rdc2tlx_narrow_widedf) +
    geom_jitter(aes(x=rdc_cluster, y=breaks_prop, color=tlx_bait_chrom, size=ifelse(tlx_is_bait_chrom, "bait", "other")), width=0.2) +
    # geom_line(aes(x=rdc_cluster_number, y=breaks_prop, color=tlx_bait_chrom)) +
    labs(y="breaks(+)/breaks(-)", x="", title=paste0("Chr8+ bait junction with 3 clusters")) +
    scale_size_manual(values=c("bait"=5, "other"=2.5)) +
    facet_grid(tlx_bait_chrom~., scales="free", space="free") +
    coord_flip()


    scale_x_continuous(breaks=unique(rdc2tlx_narrow_widedf$rdc_cluster_number), labels=unique(rdc2tlx_narrow_widedf$rdc_cluster)) +
    theme_gray(base_size=16)
  dev.off()




  rdc2tlx_narrow_density_df = rdc2tlx_narrow_df %>%
    dplyr::group_by(tlx_sample, junction_position, rdc_cluster, tlx_strand, tlx_bait_strand) %>%
    dplyr::filter(dplyr::n()>=5) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(tlx_sample, junction_position, rdc_cluster, tlx_strand, tlx_bait_strand) %>%
    dplyr::do((function(z){
      zz<<-z
      z.density = density(z$JunctionRel)
      data.frame(JunctionRel=z.density$x, density=z.density$y*nrow(z))
      cbind(z[1, c("tlx_sample", "junction_position", "rdc_cluster", "tlx_strand", "tlx_bait_strand")], data.frame(JunctionRel=z.density$x, density=z.density$y*nrow(z)), JunctionRel.max=z.density$x[which.max(z.density$y)], JunctionRel.maxdensity=max(z.density$y)*nrow(z))
    })(.)) %>%
    dplyr::group_by(tlx_sample, rdc_cluster) %>%
    dplyr::mutate(JunctionRel.maxdensity=JunctionRel.maxdensity/max(density), density=density/max(density))

  rdc2tlx_narrow_peak_df = rdc2tlx_narrow_density_df %>% dplyr::distinct(tlx_sample, junction_position, rdc_cluster, tlx_strand, .keep_all=T)
  ggplot(rdc2tlx_narrow_density_df) +
    geom_ribbon(aes(x=JunctionRel, ymin=0, ymax=density, fill=tlx_strand), alpha=0.5) +
    # geom_point(aes(x=JunctionRel.max, y=JunctionRel.maxdensity), data=rdc2tlx_narrow_peak_df)  +
    facet_grid(tlx_sample~rdc_cluster+tlx_bait_strand, scales="free_y")

  ggplot(rdc2tlx_narrow_peak_df) +
    geom_boxplot(aes(x=junction_position, y=JunctionRel.maxdensity, fill=tlx_strand)) +
    facet_wrap(~rdc_cluster+junction_position, scales="free")

  rdc2tlx_narrow_peakdiff_df = rdc2tlx_narrow_peak_df %>%
    dplyr::group_by(tlx_sample, junction_position, rdc_cluster, tlx_bait_strand) %>%
    dplyr::summarise(JunctionRel.maxdensity_diff=JunctionRel.maxdensity[tlx_strand=="+"] - JunctionRel.maxdensity[tlx_strand=="-"])

  ggplot(rdc2tlx_narrow_peakdiff_df) +
    geom_boxplot(aes(x=junction_position, y=JunctionRel.maxdensity_diff, fill=paste0(rdc_cluster, "(", junction_position, ")")))
  #
  #
  #
  ks.test(x$JunctionRel[x$junction_position=="after bait"], x$JunctionRel[x$junction_position=="before bait"])



  strand_balance_df = rdc2tlx_narrow_df %>%
    dplyr::group_by(rdc_cluster, tlx_sample) %>%
    dplyr::mutate(sample_rdc_breaks_count=dplyr::n()) %>%
    dplyr::filter(sample_rdc_breaks_count>=20) %>%
    dplyr::group_by(tlx_strand, tlx_bait_strand) %>%
    dplyr::mutate(group_breaks_count=dplyr::n()) %>%
    dplyr::group_by(rdc_cluster) %>%
    dplyr::mutate(rdc_breaks_count=dplyr::n()) %>%
    dplyr::group_by(tlx_sample, rdc_chrom, tlx_strand, tlx_bait_strand, rdc_breaks_count, junction_position) %>%
    dplyr::summarize(sample_position_breaks_proportion=dplyr::n()/unique(group_breaks_count)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(sample_position_breaks_proportion=sample_position_breaks_proportion-min(sample_position_breaks_proportion[!is.infinite(sample_position_breaks_proportion)], na.rm=T), sample_position_breaks_proportion=sample_position_breaks_proportion/max(sample_position_breaks_proportion[!is.infinite(sample_position_breaks_proportion)], na.rm=T)) %>%
    dplyr::group_by(tlx_strand, tlx_bait_strand) %>%
    dplyr::mutate(samples_n=length(unique(tlx_sample))) %>%
    dplyr::group_by(rdc_chrom, tlx_strand, tlx_bait_strand, rdc_breaks_count, samples_n, junction_position) %>%
    dplyr::summarize(position_breaks_proportion=mean(sample_position_breaks_proportion)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(bait=paste0("bait", tlx_bait_strand), junction=paste0("junction", tlx_strand), chromosomes=paste0(bait, "/", junction, " (samples:", samples_n, ")"))
  dim(strand_balance_df)
  strand_balance_df$sample_position_breaks_proportion

  max(strand_balance_df$sample_position_breaks_proportion, na.rm=T)
  ggplot(strand_balance_df) +
    ggridges::geom_density_ridges(aes(x=position_breaks_proportion, y=chromosomes), bandwidth=0.02) +
    facet_grid(bait ~ junction_position, scales="free_y") +
    labs(y="", x="breaks proportion in RDC (averaged equaly from all samples)")


  # strand_densities_df = strand_balance_df %>%
  #   dplyr::mutate(chromosomes=paste0(tlx_bait_strand, tlx_strand)) %>%
  #   # dplyr::filter(repliseqTime_max>10) %>%
  #   # dplyr::mutate(rdc_strand=paste0("Transcription", rdc_strand)) %>%
  #   dplyr::group_by(tlx_strand, tlx_bait_strand, chromosomes) %>%
  #   dplyr::do((function(z){ zz<<-z; data.frame(loc=density(z$n)$x, dens=density(z$n)$y, bw=1e6) })(.)) %>%
  #   dplyr::mutate(dens=ifelse(tlx_strand=="+", dens, -dens))
  # ggplot(strand_densities_df, aes(dens, loc, fill=tlx_strand)) +
  #   geom_polygon() +
  #   labs(y='strand proportion', fill="Junction strand", title="Transcription dependant junctions strand proportion") +
  #   theme_minimal(base_size=20) +
  #   # coord_cartesian(ylim=c(0,100)) +
  #   facet_wrap(~chromosomes, scales="free") +
  #   theme(axis.title.x=element_blank())
}

analyze_examples = function() {
  rdc_df = readr::read_tsv("data/rdc_pnas_mm10.tsv") %>%
    dplyr::filter(!grepl(",|--", rdc_gene) & rdc_strand %in% c("+", "-"))
  rdc_ranges = rdc_df %>% df2ranges(rdc_chrom, rdc_start, rdc_end)

  repliseqTime_df = readr::read_tsv("~/Workspace/Emily/data/zhao_bmc_repliseq_2020/preprocessed/repliseqTime_merged.tsv") %>%
    dplyr::filter(repliseqTime_celltype=="esc/npc") %>%
    dplyr::select(repliseqTime_chrom, repliseqTime_start, repliseqTime_end, repliseqTime_avg)
  repliseqTime_ranges = repliseqTime_df %>% df2ranges(repliseqTime_chrom, repliseqTime_start, repliseqTime_end)
  rdc2repliseq_df = as.data.frame(IRanges::mergeByOverlaps(rdc_ranges, repliseqTime_ranges))  %>%
    dplyr::select(-dplyr::matches("_ranges\\.")) %>%
    dplyr::group_by(rdc_chrom, rdc_start, rdc_end, rdc_cluster, rdc_group, rdc_strand, rdc_gene) %>%
    dplyr::summarize(repliseqTime_max=max(repliseqTime_avg, na.rm=T))
  rdc2repliseq_ranges = rdc2repliseq_df %>% df2ranges(rdc_chrom, rdc_start, rdc_end)

  samples_df = readr::read_tsv("~/Workspace/HTGTS_data/All_samples.tsv", comment="#") %>%
    dplyr::filter(experiment %in% c("APH concentration", "Wei et al. PNAS 2018")) %>%
    dplyr::mutate(path=paste0("~/Workspace/HTGTS_data/", path))
  tlx_df = tlx_read_many(samples_df)
  tlx_df = tlx_remove_rand_chromosomes(tlx_df)
  tlx_df = tlx_mark_bait_chromosome(tlx_df)
  tlx_df = tlx_mark_bait_junctions(tlx_df, 4e6)
  tlx_dff = tlx_df %>% dplyr::filter(!tlx_is_bait_junction & tlx_is_bait_chromosome)
  tlx_ranges = tlx_dff %>% df2ranges(Rname, Rstart, Rend, tlx_strand)
  tlx2rdc_df = as.data.frame(IRanges::mergeByOverlaps(rdc2repliseq_ranges, tlx_ranges)) %>%
    dplyr::select(-dplyr::matches("_ranges\\."))

  strand_balance_df = tlx2rdc_df %>%
    dplyr::group_by(rdc_gene, tlx_sample) %>%
    dplyr::mutate(sample_rdc_n=dplyr::n()) %>%
    dplyr::filter(sample_rdc_n>0) %>%
    dplyr::group_by(rdc_gene) %>%
    dplyr::mutate(rdc_n=dplyr::n()) %>%
    # dplyr::mutate(rdc_n>1) %>%
    dplyr::group_by(tlx_sample, rdc_gene, rdc_chrom, rdc_strand, tlx_strand, rdc_n, repliseqTime_max) %>%
    dplyr::summarize(sample_n=dplyr::n()/sample_rdc_n[1]) %>%
    dplyr::group_by(rdc_gene, rdc_chrom, rdc_strand, tlx_strand, rdc_n, repliseqTime_max) %>%
    dplyr::summarize(n=mean(sample_n)) %>%
    dplyr::ungroup()
  strand_densities_df = strand_balance_df %>%
    # dplyr::filter(repliseqTime_max>10) %>%
    dplyr::mutate(rdc_strand=paste0("Transcription", rdc_strand)) %>%
    dplyr::group_by(rdc_strand, tlx_strand) %>%
    dplyr::do(data.frame(loc=density(.$n)$x, dens=density(.$n)$y, bw=1e6)) %>%
    dplyr::mutate(dens=ifelse(tlx_strand=="+", dens, -dens))
  ggplot(strand_densities_df, aes(dens, loc, fill=tlx_strand)) +
    geom_polygon() +
    labs(y='strand proportion', fill="Junction strand", title="Transcription dependant junctions strand proportion") +
    theme_minimal(base_size=20) +
    # coord_cartesian(ylim=c(0,100)) +
    facet_wrap(~rdc_strand) +
    theme(axis.title.x=element_blank())

  strand_balance_df = tlx2rdc_df %>%
    dplyr::group_by(tlx_sample, rdc_gene) %>%
    dplyr::mutate(rdc_n=n()) %>%
    # dplyr::mutate(rdc_n>1) %>%
    dplyr::filter(rdc_n>20) %>%
    dplyr::group_by(tlx_sample, rdc_gene, rdc_chrom, rdc_strand, tlx_strand, rdc_n, repliseqTime_max) %>%
    dplyr::summarize(n=n()/rdc_n[1]) %>%
    # dplyr::summarize(n=n()) %>%
    dplyr::ungroup() %>%
    dplyr::filter(tlx_strand=="+")


  strand_balance_df %>%
    dplyr::mutate(tlx_strand=ifelse(tlx_strand=="+", "plus", "minus")) %>%
    reshape2::dcast(... ~ tlx_strand, value.var="n") %>%
    ggplot() +
      ggridges::geom_density_ridges(aes(x=plus, fill=rdc_chrom, y=rdc_chrom))

  ggplot(strand_balance_df) +
    geom_boxplot(aes(x=rdc_strand, y=n, fill=tlx_strand), position="dodge") +
    facet_wrap(~rdc_chrom)
}
