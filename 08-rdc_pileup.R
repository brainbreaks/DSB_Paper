library(readr)
library(dplyr)
library(ggplot2)
library(cowplot)
library(stringr)
devtools::load_all('breaktools/')
source("00-utils.R")


rdc_pileup = function()
{
  #
  # Load offtargets
  #
  offtargets_df = readr::read_tsv("data/offtargets_dkfz.tsv")

  #
  # Load genes
  #
  genes_df = gtf_read('genomes/mm10/annotation/mm10.ncbiRefSeq.gtf.gz') %>%
    dplyr::filter(gene_cluster_i<=3)

  #
  # Replication and replication termination site (annotated repli-seq data)
  #
  forks_df = readr::read_tsv("data/replication_forks_NPC.tsv") %>%
    dplyr::mutate(fork_strand=ifelse(fork_direction=="telomeric", "+", "-"))

  #
  # Load samples data
  #

  samples_df = tlx_read_paper_samples("data/htgts_samples.tsv", "data")  %>%
    dplyr::filter(!control) %>%
    dplyr::filter(
      grepl("(Ctnna2|Nrxn1) promoter/enhancer", experiment) |
      grepl("concentration", experiment) & concentration %in% 0.4 |
      grepl("Wei|Tena", experiment)
    ) %>% dplyr::mutate(group="Inter")

  #
  # Load TLX
  #
  tlx_all_df = tlx_read_many(samples_df, threads=16) %>%
    tlx_extract_bait(bait_size=19, bait_region=12e6) %>%
    tlx_calc_copynumber(bowtie2_index="genomes/mm10/mm10", max_hits=100, threads=16) %>%
    tlx_mark_offtargets(offtargets_df, offtarget_region=1e5, bait_region=1e4)
  libfactors_df = tlx_all_df %>% tlx_libsizes()

  tlx_viz_df = tlx_all_df %>%
    tlx_remove_rand_chromosomes() %>%
    dplyr::filter(!tlx_duplicated & !tlx_is_bait_junction & !tlx_is_offtarget) %>%
    dplyr::group_by(Rname) %>%
    dplyr::mutate(dbscan_cluster=dbscan::dbscan(matrix(Junction), minPts=5, eps=30)$cluster) %>%
    dplyr::filter(dbscan_cluster==0) %>%
    dplyr::ungroup()

  #
  # Load GRO-seq data
  #
  groseq_df = data.frame()
  for(f in Sys.glob("data/groseq/bedgraph5k/*.bedgraph")) {
    g = readr::read_tsv(f, col_names=c("groseq_chrom", "groseq_start", "groseq_end", "groseq_score")) %>%
      dplyr::mutate(groseq_filename=basename(f)) %>%
      dplyr::mutate(groseq_strand=gsub(".*(\\+|\\-).bedgraph", "\\1", groseq_filename)) %>%
      dplyr::mutate(groseq_sample=gsub("_bin.*$", "", groseq_filename)) %>%
      dplyr::mutate(groseq_primary=dplyr::case_when(grepl("NXP010", groseq_filename, ignore.case=T)~"NXP010", grepl("NXP047", groseq_filename, ignore.case=T)~"NXP047", T~NA_character_))
    groseq_df = dplyr::bind_rows(groseq_df, g)
  }
  groseq_sumdf = groseq_df %>%
    dplyr::filter(groseq_end-groseq_start==5e3) %>%
    # dplyr::filter(groseq_primary=="NXP047") %>%
    dplyr::group_by(groseq_chrom, groseq_start, groseq_end, groseq_strand) %>%
    dplyr::summarise(groseq_score=sum(groseq_score))

  #
  # TLX coverage
  #
  tlxcov_viz_df = tlx_viz_df %>%
    tlx_coverage(group="group", extsize=5e4, exttype="symmetrical", libfactors_df=libfactors_df, ignore.strand=F)

  tlxcov_viz_df %>%
    dplyr::group_by(tlx_group, tlx_strand) %>%
    dplyr::do((function(z) {
      z %>%
        dplyr::select(tlxcov_chrom, tlxcov_start, tlxcov_end, tlxcov_pileup) %>%
        readr::write_tsv(paste0(z$tlx_group[1], "_", z$tlx_strand[1], ".bedgraph"), col_names=F)
    })(.))


  plot_viewport = 7e5
  plot_center_strategy = "collision"

  #
  # Find middle of RDC
  #
  if(plot_center_strategy=="between peaks") {
    # Visually detect center between telomeric and centromeric peaks
    collision_df = readr::read_tsv("data/rdc.tsv") %>%
      # dplyr::filter(grepl("RDC-chr12-53.7", rdc_name)) %>%
      dplyr::filter(tlx_group=="APH-Inter" & rdc_subset=="Wei+DKFZ" & rdc_is_significant) %>%
      dplyr::mutate(rdc_region_start=rdc_extended_start-rdc_extended_length/4, rdc_region_end=rdc_extended_end+rdc_extended_length/4) %>%
      df2ranges(rdc_chrom, rdc_region_start, rdc_region_end) %>%
      innerJoinByOverlaps(tlx_viz_df %>% dplyr::filter(tlx_group=="Inter") %>% dplyr::select(-tlx_group) %>% df2ranges(Rname, Junction, Junction)) %>%
      dplyr::arrange(rdc_name, Rname, Junction) %>%
      dplyr::group_by(collision_chrom=rdc_chrom, rdc_name) %>%
      dplyr::summarise(
        collision_tz=Junction,
        sum_minus=sum(tlx_strand=="-"),
        sum_plus=sum(tlx_strand=="+"),
        count_minus=cumsum(tlx_strand=="-"),
        count_plus=sum(tlx_strand=="+")-cumsum(tlx_strand=="+"),
        # proportion=abs(log2(count_minus/count_plus)),
        # proportion=abs(log2((count_plus/sum_minus)/(count_plus/sum_plus))),
        proportion=count_minus/sum_minus*count_plus/sum_plus,
        empty=""
      ) %>%
      # dplyr::summarise(proportion=abs(log2((cumsum(tlx_strand=="-")/sum(tlx_strand=="-"))/((sum(tlx_strand=="+")-cumsum(tlx_strand=="+"))/sum(tlx_strand=="+"))))) %>%
      dplyr::filter(!is.infinite(proportion)) %>%
      dplyr::group_by(rdc_name) %>%
      dplyr::arrange(dplyr::desc(proportion)) %>%
      dplyr::slice(1) %>%
      dplyr::ungroup() %>%
      dplyr::select(collision_chrom, collision_tz)
  } else {
    if(plot_center_strategy=="collision") {
      # Collissions as predicted by tzNN
      collision_df = readr::read_tsv("data/replication_collisions_NPC.tsv")
    } else {
      stop("Unsupported value of 'plot_center_strategy'")
    }
  }


  #
  # Load RDC and add the most middle termination site position
  #
  rdc_df = readr::read_tsv("data/rdc.tsv") %>%
    dplyr::filter(tlx_group=="APH-Inter" & rdc_subset=="Wei+DKFZ" & rdc_is_significant) %>%
    df2ranges(rdc_chrom, rdc_extended_start, rdc_extended_end)  %>%
    innerJoinByOverlaps(collision_df %>% df2ranges(collision_chrom, collision_tz, collision_tz)) %>%
    dplyr::mutate(rdc_region_start=collision_tz-plot_viewport, rdc_region_end=collision_tz+plot_viewport) %>%
    dplyr::arrange(dplyr::desc(pmax(abs(rdc_extended_start-collision_tz), abs(rdc_extended_end-collision_tz)))) %>%
    dplyr::distinct(tlx_group, rdc_subset, rdc_name, .keep_all=T) %>%
    dplyr::select(dplyr::matches("rdc"), collision_tz, -rdc_subset)

  rdc2tlxcov_df = rdc_df %>%
    df2ranges(rdc_chrom, rdc_region_start, rdc_region_end) %>%
    innerJoinByOverlaps(tlxcov_viz_df %>% df2ranges(tlxcov_chrom, tlxcov_start, tlxcov_end))
  rdc_filter = rdc_df %>%
    df2ranges(rdc_chrom, rdc_start, rdc_end) %>%
    innerJoinByOverlaps(tlx_viz_df %>% df2ranges(Rname, Junction, Junction)) %>%
    dplyr::filter(!is.na(rdc_gene_strand) & !is.na(collision_tz) & !is.na(rdc_group) & rdc_group==1 & tlx_group=="Inter" & grepl("Twin", rdc_peak_shape)) %>%
    dplyr::group_by(rdc_name, tlx_group) %>%
    dplyr::mutate() %>%
    dplyr::group_by(rdc_name, tlx_group) %>%
    dplyr::summarise(junctions_count=dplyr::n()) %>%
    dplyr::filter(junctions_count > 100) %>%
    dplyr::distinct(tlx_group, rdc_name)

  #
  # Prepare ggplot2 data
  #
  rdc2tlxcov_ggplot = rdc2tlxcov_df %>%
    dplyr::inner_join(rdc_filter, by=c("tlx_group", "rdc_name")) %>%
    dplyr::mutate(tlxcov_rel_start=tlxcov_start-collision_tz, tlxcov_rel_end=tlxcov_end-collision_tz) %>%
    dplyr::mutate(rdc_name_i=as.numeric(as.factor(rdc_name))) %>%
    dplyr::group_by(rdc_name) %>%
    dplyr::mutate(tlxcov_rel_pileup=tlxcov_pileup/max(tlxcov_pileup, na.rm=T)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(tlxcov_start=tlxcov_rel_start, tlxcov_end=tlxcov_rel_end, tlxcov_pileup=tlxcov_rel_pileup)

  rdc2replication_ggplot = rdc_df %>%
    dplyr::inner_join(rdc_filter %>% dplyr::distinct(rdc_name), by="rdc_name") %>%
    df2ranges(rdc_chrom, rdc_region_start, rdc_region_end) %>%
    innerJoinByOverlaps(forks_df %>% df2ranges(fork_chrom, fork_start, fork_end)) %>%
    dplyr::mutate(replication_start=pmin(rdc_region_end, pmax(rdc_region_start, fork_iz)), replication_end=pmin(rdc_region_end, pmax(rdc_region_start, fork_tz))) %>%
    dplyr::filter(fork_end>=rdc_region_start & fork_start<=rdc_region_end) %>%
    dplyr::mutate(replication_rel_start=replication_start-collision_tz, replication_rel_end=replication_end-collision_tz)

  rdc2genes_ggplot = rdc_df %>%
    dplyr::inner_join(rdc_filter %>% dplyr::distinct(rdc_name), by="rdc_name") %>%
    df2ranges(rdc_chrom, rdc_region_start, rdc_region_end) %>%
    innerJoinByOverlaps(genes_df %>% df2ranges(gene_chrom, gene_start, gene_end)) %>%
    dplyr::mutate(gene_start=pmax(rdc_region_start, gene_start), gene_end=pmin(rdc_region_end, gene_end)) %>%
    dplyr::mutate(gene_rel_start=gene_start-collision_tz, gene_rel_end=gene_end-collision_tz) %>%
    dplyr::mutate(segment_start=ifelse(gene_strand=="+", gene_rel_start, gene_rel_end), segment_end=ifelse(gene_strand=="+", gene_rel_end, gene_rel_start))

  rdc2groseq_ggplot = rdc_df %>%
    dplyr::inner_join(rdc_filter %>% dplyr::distinct(rdc_name), by="rdc_name") %>%
    df2ranges(rdc_chrom, rdc_region_start, rdc_region_end) %>%
    innerJoinByOverlaps(groseq_sumdf %>% df2ranges(groseq_chrom, groseq_start, groseq_end)) %>%
    dplyr::mutate(groseq_rel_start=groseq_start-collision_tz, groseq_rel_end=groseq_end-collision_tz) %>%
    dplyr::group_by(rdc_chrom, rdc_region_start, rdc_region_end) %>%
    dplyr::mutate(groseq_score=groseq_score/max(groseq_score)) %>%
    dplyr::ungroup()

  plist = rdc2tlxcov_ggplot %>%
    # dplyr::filter(rdc_name=="RDC_166" & tlx_group=="Inter") %>%
    dplyr::group_split(tlx_group) %>%
    lapply(function(rdc2tlxcov_ggplot.s) {
      rdc2tlxcov_ggplot.s <<- rdc2tlxcov_ggplot.s
      key = rdc2tlxcov_ggplot.s %>% dplyr::distinct(tlx_group)
      rdc_filter.s = rdc2tlxcov_ggplot.s %>% dplyr::distinct(rdc_name)
      rdc2genes_ggplot.s = rdc2genes_ggplot %>%
        dplyr::inner_join(rdc_filter.s, by=c("rdc_name")) %>%
        dplyr::arrange(dplyr::desc(gene_length)) %>%
        dplyr::distinct(rdc_name, .keep_all=T)
      rdc2genes_ggplot.s[rdc2genes_ggplot.s$rdc_gene_strand=="-",c("segment_start", "segment_end")] = -rdc2genes_ggplot.s[rdc2genes_ggplot.s$rdc_gene_strand=="-",c("segment_start", "segment_end")]

      rdc2replication_ggplot.s = rdc2replication_ggplot %>%
        dplyr::inner_join(rdc_filter.s, by="rdc_name")
      rdc2replication_ggplot.s[rdc2replication_ggplot.s$rdc_gene_strand=="-",c("replication_rel_start", "replication_rel_end")] = -rdc2replication_ggplot.s[rdc2replication_ggplot.s$rdc_gene_strand=="-",c("replication_rel_start", "replication_rel_end")]

      rdc2tlxcov_ggplot.s = rdc2tlxcov_ggplot.s %>%
        # dplyr::arrange(rdc_bootstrap_pvalue) %>%
        dplyr::arrange(dplyr::desc(rdc_extended_length)) %>%
        dplyr::mutate(rdc_name=factor(rdc_name, unique(rdc_name)))
      rdc2tlxcov_ggplot.s[rdc2tlxcov_ggplot.s$rdc_gene_strand=="-",c("tlxcov_start", "tlxcov_end")] = -rdc2tlxcov_ggplot.s[rdc2tlxcov_ggplot.s$rdc_gene_strand=="-",c("tlxcov_start", "tlxcov_end")]

      rdc2groseq_ggplot.s = rdc2groseq_ggplot %>%
        dplyr::inner_join(rdc_filter.s, by="rdc_name")
      rdc2groseq_ggplot.s[rdc2groseq_ggplot.s$rdc_gene_strand=="-",c("groseq_rel_start", "groseq_rel_end")] = -rdc2groseq_ggplot.s[rdc2groseq_ggplot.s$rdc_gene_strand=="-",c("groseq_rel_start", "groseq_rel_end")]

      strand_colors = c("+"="#FF0000", "-"="#0000FF")
      y_dist = 1.5
      g = ggplot() +
        geom_vline(xintercept=0, size=0.3, alpha=0.4) +
        # geom_bar(aes(x=groseq_rel_start, y=-groseq_score, fill=groseq_strand), data=rdc2groseq_ggplot.s, stat="identity", alpha=0.5) +
        geom_tlxcov(rdc2tlxcov_ggplot.s, alpha=0.4) +
        geom_segment(aes(x=segment_start, xend=segment_end, y=-0.1*y_dist*gene_cluster_i, yend=-0.1*y_dist*gene_cluster_i), size=0.3, arrow=grid::arrow(length=unit(4,"pt")), data=rdc2genes_ggplot.s) +
        # geom_segment(aes(x=segment_start, xend=segment_start, y=-0.1*gene_cluster_i-0.05, yend=-0.1*gene_cluster_i+0.05), size=0.3, data=rdc2genes_ggplot.s) +
        geom_text(aes(x=(segment_start+segment_end)/2, y=0.2*(gene_cluster_i+1)+0.05, label=gene_id), data=rdc2genes_ggplot.s %>% dplyr::distinct(rdc_name, .keep_all=T), size=3, vjust=1, hjust=0, color="#3F2B2D", fontface="bold") +
        # geom_text(aes(x=(gene_rel_end+gene_rel_start)/2, y=0.2*(gene_cluster_i+1)+0.05, label=gene_id), data=rdc2genes_ggplot.s %>% dplyr::distinct(rdc_name, .keep_all=T), size=3, vjust=1, hjust=0, color="#000000") +
        geom_segment(aes(x=replication_rel_start, xend=replication_rel_end, y=-0.2*y_dist+y, yend=-0.2*y_dist+y, color=fork_strand), size=0.1, arrow=grid::arrow(length=unit(4,"pt")), data=rdc2replication_ggplot.s %>% dplyr::mutate(y=ifelse(fork_strand=="+", 0.01, -0.01))) +
        facet_grid(rdc_name~., scales="free_y") +
        labs(x="Relative distance from RTS (100kb)", y="", color="Replication direction", fill="Translocation orientation", title=paste0(names(key), ": ", key, collapse="    ")) +
        scale_x_continuous(labels=scales::label_number(accuracy=1, scale=1e-5, suffix="")) +
        scale_fill_manual(values=strand_colors) +
        scale_color_manual(values=strand_colors) +
        coord_cartesian(xlim=c(-plot_viewport, plot_viewport)) +
        theme_classic(base_size=12) +
        theme(
          strip.background=element_blank(),
          strip.text.y=element_text(size=10, angle=0),
          legend.position="right",
          axis.line.y = element_blank(),
          axis.text.y=ggplot2::element_blank(),
          axis.ticks.y=ggplot2::element_blank(),
          axis.title.y=ggplot2::element_blank(),
          legend.key.size=unit(1.0, 'cm'),
          panel.spacing = unit(0, "lines")
        )

      # if(s<length(plist_strands)) {
      #   g = g + theme(axis.title.x=ggplot2::element_blank())
      # }
      g
  })

  pdf("reports/06-rdc_pileup/tlxcov_relative.pdf", width=8.27, height=11.69*2)
  plist
  dev.off()
}