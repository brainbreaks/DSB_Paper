library(introdataviz)
library(ggplot2)
source("00-utils.R")
devtools::load_all("breaktools/")

replication_fork_length = function() {
  dir.create("reports/09-replication_fork")

  #
  # Load Repliseq
  #
  repliseq_df = readr::read_tsv("data/repliseq_zhao_bmc2020/repliseq_ESC.tsv") %>% dplyr::mutate(Run="Zhao_mESC")

  #
  # Load genes
  #
  genes_df = gtf_read('genomes/mm10/annotation/mm10.ncbiRefSeq.gtf.gz')
  genes_reduced_df = genes_df %>%
    df2ranges(gene_chrom, gene_start, gene_end) %>%
    GenomicRanges::reduce() %>%
    as.data.frame() %>%
    dplyr::select(gene_chrom=seqnames, gene_start=start, gene_end=end) %>%
    dplyr::mutate(gene_length=gene_end-gene_start)

  #
  # Load GRO-seq data
  #
  groseq_df = data.frame()
  for(f in Sys.glob("data/groseq/bedgraph50k/*.bedgraph")) {
    g = readr::read_tsv(f, col_names=c("groseq_chrom", "groseq_start", "groseq_end", "groseq_score")) %>%
      dplyr::mutate(groseq_filename=basename(f)) %>%
      dplyr::mutate(groseq_strand=gsub(".*(\\+|\\-).bedgraph", "\\1", groseq_filename)) %>%
      dplyr::mutate(groseq_sample=gsub("_bin.*$", "", groseq_filename)) %>%
      dplyr::mutate(groseq_primary=dplyr::case_when(grepl("NXP010", groseq_filename, ignore.case=T)~"NXP010", grepl("NXP047", groseq_filename, ignore.case=T)~"NXP047", T~NA_character_))
    groseq_df = dplyr::bind_rows(groseq_df, g)
  }
  groseq_sumdf = groseq_df %>%
    dplyr::filter(groseq_end-groseq_start==50e3) %>%
    dplyr::filter(groseq_primary=="NXP047") %>%
    dplyr::group_by(groseq_chrom, groseq_start, groseq_end) %>%
    dplyr::summarise(groseq_score=sum(groseq_score))

  #
  # Repliseq data
  #
  forks_df = readr::read_tsv("data/replication_forks_NPC.tsv") %>%
    dplyr::mutate(fork_length=fork_end-fork_start)
  fork_collisions_df = readr::read_tsv("data/replication_collisions_NPC.tsv")


  #
  # Load RDC and test to filter out only the RDC that are significant in concentration samples
  #
  rdc_df = readr::read_tsv("data/rdc.tsv") %>%
    dplyr::filter(tlx_group=="APH-Inter" & rdc_subset=="Wei+DKFZ" & rdc_is_significant) %>%
    dplyr::select(-tlx_group, -rdc_subset)

  groseq_cutoff = rdc_df %>%
    df2ranges(rdc_chrom, (rdc_end+rdc_start)/2-5e5, (rdc_end+rdc_start)/2+1e5) %>%
    leftJoinByOverlaps(groseq_sumdf %>% df2ranges(groseq_chrom, groseq_start, groseq_end)) %>%
    dplyr::group_by(rdc_chrom, rdc_end, rdc_start) %>%
    dplyr::summarise(groseq_score=max(groseq_score))
  groseq_cutoff = quantile(groseq_cutoff$groseq_score, 0.005)


  subset_palette = c("RDC"="#00B9D3", "Gene>100kb"="#D34B00", "Gene<100kb"="#D39B00", "No gene"="#A2D300")
  rdc2replication_df = fork_collisions_df %>%
    df2ranges(collision_chrom, collision_iz_left, collision_iz_right) %>%
    leftJoinByOverlaps(genes_reduced_df %>% df2ranges(gene_chrom, gene_start, gene_end)) %>%
    dplyr::group_by(collision_chrom, collision_iz_left, collision_iz_right) %>%
    dplyr::summarise(gene_length=max(tidyr::replace_na(gene_length, 0), na.rm=T)) %>%
    df2ranges(collision_chrom, collision_iz_left, collision_iz_right) %>%
    leftJoinByOverlaps(rdc_df %>% df2ranges(rdc_chrom, rdc_extended_start, rdc_extended_end)) %>%
    dplyr::group_by(collision_chrom, collision_iz_left, collision_iz_right, gene_length) %>%
    dplyr::summarise(rdc_length=max(tidyr::replace_na(rdc_length, 0), na.rm=T), rdc_name=ifelse(any(!is.na(rdc_name)), rdc_name, NA_character_), rdc_group=ifelse(!is.na(rdc_group), na.omit(rdc_group)[1], NA_integer_)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(overlap_type=dplyr::case_when(rdc_length>0~"RDC", gene_length>100e3~"Gene>100kb", gene_length>0~"Gene<100kb", T~"No gene")) %>%
    df2ranges(collision_chrom, collision_iz_left, collision_iz_right) %>%
    leftJoinByOverlaps(groseq_sumdf %>% df2ranges(groseq_chrom, groseq_start, groseq_end)) %>%
    dplyr::mutate(overlap_type=factor(overlap_type, names(subset_palette))) %>%
    dplyr::mutate(tr=groseq_score>=groseq_cutoff) %>%
    # dplyr::filter(replication_chrom =="chr1", replication_start ==  20859897, replication_end == 21698599) %>%
    dplyr::group_by(collision_chrom, collision_iz_left, collision_iz_right, overlap_type, gene_length, rdc_group, rdc_name) %>%
    dplyr::do((function(z) {
      z = as.data.frame(z)
      zz<<-z
      if(any(!is.na(z$tr))) {
        is_any_transcribed = any(z$tr)
        consecutive_transcribed = max(unlist(rle(z$tr)$lengths[rle(z$tr)$values]))
        return(data.frame(
          is_any_transcribed=is_any_transcribed,
          consecutive_transcribed=consecutive_transcribed,
          sum_transcribed=sum(z$tr),
          length_transcribed=sum((z$groseq_end-z$groseq_start)[z$tr]),
          groseq_median = median(z$groseq_score[z$tr], na.rm=T),
          groseq_max = max(z$groseq_score[z$tr], na.rm=T),
          prop_transcribed=mean(z$tr),
          bins_number=nrow(z)) %>%
          dplyr::mutate(is_transcribed=prop_transcribed==1 || consecutive_transcribed>=4))
      } else {
        return(data.frame(is_any_transcribed=F, consecutive_transcribed=0, prop_transcribed=0, is_transcribed=F, sum_transcribed=0, bins_number=nrow(z)))
      }
    })(.)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(collision_forks_length=collision_iz_right-collision_iz_left, collision_forks_length_group=dplyr::case_when(collision_forks_length>=1e6~">1Mb", collision_forks_length>=5e5~">500kb", collision_forks_length>=1e5~">100kb", T~">0")) %>%
    dplyr::mutate(collision_forks_length.log10=log10(collision_forks_length))


  pdf("reports/09-replication_fork/replication_fork_collision-new.pdf", width=8.27, height=11.6)
  g0 = ggplot(forks_df) +
    geom_histogram(aes(y=nrow(forks_df)*100e3*..density.., x=fork_length), alpha=0.5, binwidth=100e3, position="identity", fill="#111111", color="#000000") +
    geom_vline(xintercept=median(forks_df$fork_length), color="#FF0000") +
    geom_text(x=median(forks_df$fork_length), y=1.5e-7, label=round(median(forks_df$fork_length), 2), vjust=1.1, hjust=-0.1, size=3, color="#FF0000", data=data.frame(1)) +
    theme_paper() +
    labs(x="Replication length", y="count")

  rdc2replication_sumdf = rdc2replication_df %>%
    dplyr::group_by(overlap_type) %>%
    dplyr::summarise(sum_transcribed=sum(sum_transcribed), bins_number=sum(bins_number), prop=sum_transcribed/bins_number)
  g1 = ggplot(rdc2replication_sumdf) +
    geom_bar(aes(x=overlap_type, y=prop, fill=overlap_type), stat="identity", alpha=0.5) +
    geom_text(aes(x=overlap_type, y=prop, label=round(prop, 2)), size=3, color="#FF0000", stat="identity", vjust=ifelse(max(rdc2replication_sumdf$prop)==rdc2replication_sumdf$prop, 1.2, -1)) +
    labs(x="", y="Average\ntranscription proportion", fill="Overlap") +
    guides(fill=guide_legend(title.position="top", title.hjust=0.5, ncol=2), color="none") +
    scale_fill_manual(values=subset_palette) +
    scale_y_continuous(expand=c(0, 0), breaks=scales::pretty_breaks(n=5)) +
    theme_paper(base_size=12) +
    theme(legend.position="bottom") +
    ggplot2::theme(axis.text.x=ggplot2::element_text(angle=45, hjust=1, vjust=1), legend.key=element_rect(color="#000000"))

  g2 = ggplot(rdc2replication_df, aes(x=overlap_type, y=prop_transcribed)) +
    geom_boxplot(aes(fill=overlap_type), outlier.shape = NA, outlier.alpha=0, alpha=0.5) +
    # geom_point(aes(group=overlap_type), position=position_jitter(height=0,width = 0.1, seed=2), shape=21, color="#00000033", fill="#FFFFFF33", size=2, data=rdc2replication_df %>% dplyr::filter(abs(prop_transcribed-0.5)==0.5)) +
    geom_dotplot(aes(group=overlap_type), binwidth=0.005, binaxis="y", stackdir="center", shape=21, color="#00000033", fill="#FFFFFF33", data=rdc2replication_df %>% dplyr::filter(abs(prop_transcribed-0.5)==0.5) %>% dplyr::mutate(prop_transcribed=jitter(prop_transcribed, 0.11))) +
    geom_dotplot(aes(group=overlap_type), binwidth=0.005, binaxis="y", stackdir="center", shape=21, color="#00000033", fill="#FFFFFF33", data=rdc2replication_df %>% dplyr::filter(abs(prop_transcribed-0.5)<0.5)) +
    geom_text(aes(x=overlap_type, y=prop_transcribed.values, label=round(prop_transcribed.values, 2)), size=3, stat="identity", color="#FF0000", vjust=-1, hjust=-1, data=rdc2replication_df %>% dplyr::group_by(overlap_type) %>% dplyr::summarise(prop_transcribed.quantiles=1:3*0.25, prop_transcribed.values=quantile(prop_transcribed, prop_transcribed.quantiles, na.rm=T))) +
    guides(fill=guide_legend(title.position="top", title.hjust=0.5, ncol=2), color="none") +
    scale_fill_manual(values=subset_palette) +
    labs(x="", y="Transcribed proportion", fill="Overlap") +
    theme_paper(base_size=12) +
    theme(legend.position="bottom") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle=45, hjust=1, vjust=1))

  g3 = ggplot(rdc2replication_df, aes(y=overlap_type, x=prop_transcribed, fill=overlap_type)) +
    ggridges::geom_density_ridges(aes(height=0.1*..density.., y=overlap_type, fill=overlap_type), alpha=0.5, binwidth=0.1, position="identity", stat="binline", draw_baseline=F) +
    # geom_point(aes(color=overlap_type), pch=1, position=position_jitterdodge(jitter.width=0.01, dodge.width=0.04, seed=2), fill="#000000", alpha=0.1) +
    geom_text(label="|", size=2, position=position_jitter(width=0.01, height=0, seed=2), fill="#000000", alpha=0.1) +
    geom_text(aes(label=label), size=3, color="#FF0000", data=rdc2replication_df %>% dplyr::group_by(overlap_type) %>% dplyr::summarise(prop_transcribed=median(prop_transcribed, na.rm=T), label=paste0("\n|\n", round(prop_transcribed, 2)))) +
    scale_color_manual(values=sapply(subset_palette, function(z) "#000000")) +
    scale_fill_manual(values=subset_palette) +
    scale_x_continuous(expand=c(0, 0), breaks=scales::pretty_breaks(n=5)) +
    labs(y="", x="Transcribed proportion", fill="Overlap") +
    guides(fill=guide_legend(title.position="top", title.hjust=0.5, ncol=2), color="none") +
    theme_paper(base_size=12) +
    theme(legend.position="bottom")

  g4 = ggplot(rdc2replication_df, aes(x=overlap_type, y=collision_forks_length.log10)) +
    geom_boxplot(aes(fill=overlap_type), outlier.shape = NA, outlier.alpha=0, alpha=0.5) +
    # geom_point(aes(group=overlap_type), position=position_jitter(height=0,width = 0.1, seed=2), shape=21, color="#00000033", fill="#FFFFFF33") +
    geom_dotplot(aes(group=overlap_type), binwidth=0.005, binaxis="y", stackdir="center", shape=21, color="#00000033", fill="#FFFFFF33", data=rdc2replication_df %>% dplyr::mutate(collision_forks_length.log10=jitter(collision_forks_length.log10, 10))) +
    geom_text(aes(x=overlap_type, y=collision_forks_length.log10.values, label=round(collision_forks_length.log10.values, 2)), size=3, stat="identity", color="#FF0000", vjust=-1, hjust=-1, data=rdc2replication_df %>% dplyr::group_by(overlap_type) %>% dplyr::summarise(collision_forks_length.log10.quantiles=1:3*0.25, collision_forks_length.log10.values=quantile(collision_forks_length.log10, collision_forks_length.log10.quantiles, na.rm=T))) +
    scale_fill_manual(values=subset_palette) +
    labs(x="", y="Fork length, log10", fill="Overlap") +
    guides(fill=guide_legend(title.position="top", title.hjust=0.5, ncol=2), color="none") +
    theme_paper(base_size=12) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle=45, hjust=1, vjust=1), legend.position="bottom")

  g5 = ggplot(rdc2replication_df, aes(y=overlap_type, x=collision_forks_length.log10)) +
    ggridges::geom_density_ridges(aes(height=0.1*..density.., y=overlap_type, fill=overlap_type), alpha=0.5, binwidth=0.1, position="identity", stat="binline", draw_baseline=F) +
    geom_text(label="|", size=2, position=position_jitter(width=0.05, height=0, seed=2), fill="#000000", alpha=0.1) +
    geom_text(aes(label=label), size=3, color="#FF0000", data=rdc2replication_df %>% dplyr::group_by(overlap_type) %>% dplyr::summarise(collision_forks_length.log10=median(collision_forks_length.log10, na.rm=T), label=paste0("\n|\n", round(collision_forks_length.log10, 2)))) +
    labs(x="Replication length, log10", y="", fill="Overlap") +
    # scale_y_continuous(labels = scales::percent) + +
    scale_x_continuous(expand=c(0, 0), breaks=scales::pretty_breaks(n=5)) +
    scale_fill_manual(values=subset_palette) +
    guides(fill=guide_legend(title.position="top", title.hjust=0.5, nrow=1), color="none") +
    theme_paper(base_size=12) +
    theme(legend.position="bottom")

  gridExtra::grid.arrange(
    gridExtra::arrangeGrob(
        grid::textGrob(paste0("Each point is one of ", nrow(rdc2replication_df)*2, " replication forks"), gp=grid::gpar(fontsize=12,font=3)),
        g0,
        cowplot::get_legend(g5),
      nrow=3, heights=c(1,4,1)),
    g1 + theme(legend.position="none"),
    g2 + theme(legend.position="none"),
    g3 + theme(legend.position="none"),
    g4 + theme(legend.position="none"),
    g5 + theme(legend.position="none"),
    ncol=2, padding=4)
  dev.off()
}