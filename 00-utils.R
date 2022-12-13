tlx_read_paper_samples = function(annotation_path, data_path)
{
  tlx_read_samples(annotation_path, data_path) %>%
    dplyr::filter(
        organism=="mouse" & celltype=="NPC" & sample!="VI035" & (
          experiment=="APH concentration" |
          grepl("Nrxn1|Ctnna2 promoter/enhancer", experiment) & grepl(".*\\((NXP047||NXP010|22|22/37|22/5|47/5|18/4|38/3)\\)", group) |
          grepl("Wei|Tena", experiment)
        )
   )
}


tlx_offtarget_libfactor = function(tlx_df, offtargets_df)
{
  #
  # Off-target based normalization calculation
  #
  offtargets_best_df = offtargets_df %>%
    dplyr::arrange(offtarget_strand_pvalue) %>%
    dplyr::group_by(offtarget_bait_name) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup() %>%
    dplyr::select(offtarget_bait_name, offtarget_chrom, offtarget_start, offtarget_end, offtarget_end, offtarget_strand_pvalue)
  tlx_offtargets_df = tlx_df %>%
    dplyr::filter(tlx_is_offtarget) %>%
    df2ranges(Rname, Junction, Junction) %>%
    innerJoinByOverlaps(offtargets_best_df %>% df2ranges(offtarget_chrom, (offtarget_start+offtarget_end)/2-5e3, (offtarget_start+offtarget_end)/2+5e3)) %>%
    dplyr::filter(offtarget_bait_name==bait_name)
  libfactors_centration_df = tlx_offtargets_df %>%
    tlx_libsizes() %>%
    tlx_libfactors_between(min(library_size)/library_size)

  list(libfactors=libfactors_centration_df, offtargets=offtargets_best_df)
}

ggplot_rdc_breaks = function(rdc2rdc_df, rdc2tlxcov_df, rdc2tlx_df, rdc2genes_df, rdc2repliseq_df, rdc2replication_df=NULL, rdc2offtargets_df=NULL, rdc_filter=".*")
{
  rdc2rdc_dff = rdc2rdc_df %>% dplyr::filter(grepl(rdc_filter, rdc_name))
  rdc2tlxcov_dff = rdc2tlxcov_df %>% dplyr::filter(grepl(rdc_filter, rdc_name))
  rdc2tlx_dff = rdc2tlx_df %>% dplyr::filter(grepl(rdc_filter, rdc_name))
  rdc2genes_dff = rdc2genes_df %>% dplyr::filter(grepl(rdc_filter, rdc_name))
  rdc2repliseq_dff = rdc2repliseq_df %>% dplyr::filter(grepl(rdc_filter, rdc_name))

  rdc2genes_dff = rdc2genes_dff %>% dplyr::mutate(segment_start=ifelse(gene_strand=="+", gene_start, gene_end), segment_end=ifelse(gene_strand=="+", gene_end, gene_start))
  rdc2genes_dff.long = rdc2genes_dff %>% dplyr::filter(gene_length>=2e5)
  rdc2genes_dff.short = rdc2genes_dff %>% dplyr::filter(gene_length<2e5)

  g = ggplot()
  if(!is.null(rdc2replication_df)) {
    rdc2replication_dff = rdc2replication_df %>% dplyr::filter(grepl(rdc_filter, rdc_name))
    if(nrow(rdc2replication_dff)>0) {
      g = g +
        geom_segment(aes(x=fork_iz, xend=fork_tz, y=-19+y, yend=-19+y, color=fork_strand), size=0.1, arrow=grid::arrow(length=unit(1,"pt")), data=rdc2replication_dff %>% dplyr::mutate(fork_strand=ifelse(fork_direction=="telomeric", "+", "-"), y=ifelse(fork_direction=="telomeric", 0.2, -0.2)))

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
    rdc2offtargets_dff = rdc2offtargets_df %>% dplyr::filter(grepl(rdc_filter, rdc_name))
    if(nrow(rdc2offtargets_dff)>0) {
      g = g +
        geom_rect(aes(ymin=-25.5, ymax=-20.5, xmin=offtarget_start-1e5, xmax=offtarget_end+1e5), fill="#666666", data=rdc2offtargets_dff, size=0.1) +
        geom_rect(aes(ymin=-Inf, ymax=Inf,    xmin=offtarget_start,      xmax=offtarget_end), alpha=0.3, fill="#FF0000", data=rdc2offtargets_dff, size=0.1)
    }
  }

  g = g +
    geom_segment(aes(x=segment_start, xend=segment_end, y=-19-gene_cluster_i*4-0.2, yend=-19-gene_cluster_i*4-0.2), size=0.3, arrow=grid::arrow(length=unit(0.8,"pt")), color="#CCCCCC", data=rdc2genes_dff.short) +
    geom_segment(aes(x=segment_start, xend=segment_end, y=-19-gene_cluster_i*4-0.2, yend=-19-gene_cluster_i*4-0.2), size=0.3, arrow=grid::arrow(length=unit(0.8,"pt")), color=ifelse(rdc2genes_dff.long$gene_strand=="+", "#a14f4f", "#4f87a1"), data=rdc2genes_dff.long) +
    ggrepel::geom_text_repel(aes(x=(gene_end+gene_start)/2, y=-19-gene_cluster_i*4-2, label=gene_id), data=rdc2genes_dff.long, size=1.5, color="#FF6500", direction="x") +
    geom_hline(yintercept=-0.2, size=0.1) +
    facet_wrap(~rdc_name, scales="free_x") +
    guides(fill=guide_legend(nrow=1, byrow=TRUE)) +
    ggpubr::theme_pubclean(base_size=5) +
    theme(legend.position="bottom", axis.text.x=element_text(angle=45, hjust=1)) +
    scale_fill_manual(values=c("-"="#A6CEE3", "+"="#FB9A99", "cluster"="#00FF00")) +
    scale_color_manual(values=c("-"="#E31A1C", "+"="#1F78B4", "cluster"="#00FF00")) +
    scale_y_continuous(breaks=c(-3, c(1, 16), 27.5)-19, labels=c("gene", "early", "late", "density")) +
    scale_x_continuous(labels=scales::label_number(accuracy=0.1, scale=1e-6, suffix="Mb")) +
    labs(x="", y="")
  g
}

theme_paper = function(base_size=12) {
  theme_bw(base_size=base_size) +
  theme(
    strip.background=element_blank(),
    strip.text=element_text(size=base_size*1.5),
    axis.line = element_line(colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank())
}

theme_x_factors = function(size=NULL) {
  theme(axis.title.x=element_blank(), axis.text.x=element_text(angle=45, vjust=1, hjust=1, size=size))
}

theme_x_blank = function(size=NULL) {
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
}
