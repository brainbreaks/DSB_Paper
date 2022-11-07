tlx_offtarget_libfactor = function(tlx_df, offtargets_df)
{
  #
  # Off-target based normalization calculation
  #
  offtargets_best_ranges = offtargets_df %>%
    dplyr::arrange(offtarget_strand_pvalue) %>%
    dplyr::group_by(offtarget_bait_name) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup() %>%
    dplyr::select(offtarget_bait_name, offtarget_chrom, offtarget_start, offtarget_end, offtarget_end, offtarget_strand_pvalue) %>%
    df2ranges(offtarget_chrom, (offtarget_start+offtarget_end)/2-5e3, (offtarget_start+offtarget_end)/2+5e3)
  tlx_offtargets_df = tlx_df %>%
    dplyr::filter(tlx_is_offtarget) %>%
    df2ranges(Rname, Junction, Junction) %>%
    innerJoinByOverlaps(offtargets_best_ranges) %>%
    dplyr::filter(offtarget_bait_name==bait_name)
  libfactors_centration_df = tlx_offtargets_df %>%
    tlx_libsizes() %>%
    tlx_libfactors_between(min(library_size)/library_size)

  list(libfactors=libfactors_centration_df, offtargets=offtargets_best_df)
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
