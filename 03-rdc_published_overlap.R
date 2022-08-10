library(readr)
library(dplyr)
library(ggplot2)
library(grid)
devtools::load_all("breaktools/")

rdc_published_overlap = function()
{
  rdc_df = readr::read_tsv("data/rdc.tsv")
  rdc_combined_df = rdc_df %>%
    dplyr::filter(grepl("-(Inter|Intra)", tlx_group)) %>% dplyr::mutate(tlx_group=gsub("Inter|Intra", "Combined", tlx_group)) %>%
    dplyr::bind_rows(rdc_df) %>%
    dplyr::group_by(tlx_group, rdc_subset) %>%
    dplyr::mutate(rdc_name=paste0("RDC_", stringr::str_pad((0:(dplyr::n()))[-1], 3, pad="0")), rdc_length=rdc_end-rdc_start, rdc_extended_length=rdc_extended_end-rdc_extended_start) %>%
    dplyr::ungroup()

  #
  # Load public RDC dataset
  #
  pubrdc_raw_df = readr::read_tsv("data/pubrdc.tsv") %>%
    dplyr::filter(pubrdc_source %in% c("Wei2018", "Wei2018_DMSO") & pubrdc_celline %in% c("NPC", "NSPC")) %>%
    # dplyr::filter(pubrdc_rrs>=7) %>%
    tidyr::separate_rows(rdcpub_bait_chrom, rdcpub_bait_chrom, sep=", ?") %>%
    dplyr::group_by(pubrdc_chrom, pubrdc_start, pubrdc_end) %>%
    dplyr::mutate(pubrdc_bait_count=length(unique(rdcpub_bait_chrom))) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(tlx_group=dplyr::case_when(
      pubrdc_source=="Wei2018_DMSO" & rdcpub_bait_chrom==pubrdc_chrom ~ "DMSO-Intra",
      pubrdc_source=="Wei2018_DMSO" & rdcpub_bait_chrom!=pubrdc_chrom ~ "DMSO-Inter",
      rdcpub_bait_chrom==pubrdc_chrom ~ "APH-Intra",
      rdcpub_bait_chrom!=pubrdc_chrom ~ "APH-Inter")) %>%
    dplyr::group_by(pubrdc_source, tlx_group, pubrdc_chrom, pubrdc_start, pubrdc_end, pubrdc_bait_count) %>%
    dplyr::mutate(pubrdc_groupbait_count=dplyr::n()) %>%
    dplyr::ungroup() %>%
    dplyr::filter(!is.na(pubrdc_start))
  pubrdc_df = dplyr::bind_rows(
    pubrdc_raw_df,
    pubrdc_raw_df %>% dplyr::mutate(tlx_group=gsub("APH", "DMSO", tlx_group)),
    pubrdc_raw_df %>% dplyr::mutate(tlx_group=dplyr::case_when(pubrdc_source=="Wei2018_DMSO" ~ "DMSO",          pubrdc_source %in% c("Tena2020", "Wei2018") ~ "APH")),
    pubrdc_raw_df %>% dplyr::mutate(tlx_group=dplyr::case_when(pubrdc_source=="Wei2018_DMSO" ~ "DMSO-Combined", pubrdc_source %in% c("Tena2020", "Wei2018") ~ "APH-Combined"))) %>%
    dplyr::select(tlx_group, pubrdc_source, pubrdc_celline, pubrdc_chrom, pubrdc_start, pubrdc_end)

  #
  # Combine public and our RDC datasets
  #
  allrdc_df = dplyr::bind_rows(
    pubrdc_df %>%
      dplyr::mutate(rdc_chrom=pubrdc_chrom, rdc_start=pubrdc_start, rdc_end=pubrdc_end, rdc_extended_start=pubrdc_start, rdc_extended_end=pubrdc_end, rdc_is_significant=T) %>%
      dplyr::mutate(rdc_source=paste0(pubrdc_source, "-", pubrdc_celline), rdc_significant_sumarea=1e7, rdc_bootstrap_pvalue=min(rdc_combined_df$rdc_bootstrap_pvalue), rdc_bootstrap_zscore=max(rdc_combined_df$rdc_bootstrap_zscore), rdc_bootstrap_fc=max(rdc_combined_df$rdc_bootstrap_fc), rdc_maxscore=max(rdc_combined_df$rdc_maxscore), rdc_sumscore=max(rdc_combined_df$rdc_sumscore), rdc_length=max(rdc_combined_df$rdc_length), rdc_extended_length=max(rdc_combined_df$rdc_extended_length), rdc_maxscore.adjusted=max(rdc_combined_df$rdc_maxscore.adjusted)) %>%
      tidyr::crossing(data.frame(rdc_subset=unique(rdc_combined_df$rdc_subset))),
    rdc_combined_df %>%
        tidyr::crossing(rdc_source="DKFZ") %>%
        dplyr::select(tlx_group, rdc_chrom, rdc_extended_start, rdc_extended_end, rdc_start, rdc_end, rdc_source, rdc_subset, rdc_is_significant, rdc_bootstrap_zscore, rdc_bootstrap_pvalue, rdc_bootstrap_fc, rdc_significant_sumarea, rdc_maxscore, rdc_maxscore.adjusted, rdc_sumscore, rdc_length, rdc_extended_length)) %>%
    dplyr::mutate(rdc_bait_chrom=as.character(ifelse(grepl("Intra", tlx_group), as.character(rdc_chrom), "Other chromosome")))

  #
  # Write debugging information from final list of RDC (includding published
  #
  if(F)
  {
    allrdc_reduced_df %>%
      dplyr::select(rdc_reduced_chrom, rdc_reduced_start, rdc_reduced_end) %>%
      readr::write_tsv("reports/detect_rdc/allrdc_reduced.bed", col_names=F)

    allrdc_df %>%
      dplyr::group_by(rdc_source, tlx_group) %>%
      dplyr::do((function(df){
        df %>%
          dplyr::select(rdc_chrom, rdc_start, rdc_end) %>%
          readr::write_tsv(paste0("reports/detect_rdc/allrdc_-", df$rdc_source[1], "-", df$tlx_group[1], ".bed"), col_names=F)
      })(.))
  }


  allrdc_ranges = allrdc_df %>% df2ranges(rdc_chrom, rdc_extended_start, rdc_extended_end)
  allrdc_reduced_df = allrdc_ranges %>%
    GenomicRanges::reduce(min.gapwidth=200e3) %>%
    as.data.frame() %>%
    dplyr::mutate(rdc_reduced_id=paste0(seqnames, ":", start, "-", end)) %>%
    dplyr::select(rdc_reduced_chrom=seqnames, rdc_reduced_start=start, rdc_reduced_end=end, rdc_reduced_id)
  overlaps_df = allrdc_reduced_df %>%
    df2ranges(rdc_reduced_chrom, rdc_reduced_start, rdc_reduced_end) %>%
    innerJoinByOverlaps(allrdc_ranges)

  if(F) {
    overlaps_venn_mean = overlaps_df %>%
      dplyr::filter(grepl("APH", tlx_group)) %>%
      # dplyr::filter(rdc_length>=100e3 & rdc_bait_chrom %in% c("chr5", "chr6", "chr8", "chr17", "Other chromosome")) %>%
      dplyr::filter(rdc_significant_sumarea>=100e3 & rdc_bait_chrom %in% c("chr5", "chr6", "chr8", "Other chromosome")) %>%
      dplyr::mutate(th=-log10(rdc_bootstrap_pvalue)) %>%
      tidyr::crossing(cur_cutoff=seq(2, 100, length.out=100)) %>%
      dplyr::group_by(tlx_group, rdc_subset, cur_cutoff) %>%
      dplyr::do((function(odf) {
        ooo <<- odf
        is_aph = grepl("APH", odf$tlx_group[1])
        tlx_group_short = gsub("(APH|DMSO)-?", "", odf$tlx_group[1])
        tlx_group_short = ifelse(tlx_group_short=="", "All", tlx_group_short)

        if(all(odf$th<odf$cur_cutoff)) {
          d = data.frame(tlx_group_short=tlx_group_short, condition=ifelse(is_aph, "APH", "DMSO"), fpr=0, tpr=0)
        } else {
          x = odf %>%
            dplyr::mutate(tlx_group_short=tlx_group_short) %>%
            dplyr::filter(is_aph & rdc_source=="Wei2018-NPC" | !is_aph & rdc_source=="Wei2018_DMSO-NPC" | rdc_source=="DKFZ") %>%
            dplyr::filter(th>=cur_cutoff | rdc_source!="DKFZ") %>%
            reshape2::dcast(rdc_reduced_id+tlx_group_short~rdc_source, value.var="cur_cutoff", fun.aggregate=function(y) length(y)) %>%
            tibble::column_to_rownames("rdc_reduced_id")
          if("DKFZ" %in% colnames(x)) {
            d = data.frame(
              tlx_group_short=tlx_group_short,
              condition=ifelse(is_aph, "APH", "DMSO"),
              fpr=sum(x[["DKFZ"]]>0 & x[,grepl("Wei", colnames(x))]==0)/sum(x[["DKFZ"]]>0),
              tpr=sum(x[["DKFZ"]]>0 & x[,grepl("Wei", colnames(x))]>0)/sum(x[,grepl("Wei", colnames(x))]>0)
            )
          } else {
            d = data.frame(tlx_group_short=tlx_group_short, condition=ifelse(is_aph, "APH", "DMSO"), fpr=0, tpr=0)
          }
        }

        d
      })(.))
      ggplot(overlaps_venn_mean %>% dplyr::filter(cur_cutoff<=5e6)) +
        geom_line(aes(x=cur_cutoff, y=fpr, color=tlx_group_short, linetype="fpr")) +
        geom_line(aes(x=cur_cutoff, y=tpr, color=tlx_group_short, linetype="tpr")) +
        geom_vline(xintercept=0.6) +
        # coord_cartesian(xlim=c(0, 1e6), ylim=c(0, 1)) +
        facet_grid(condition~rdc_subset, scales="free") +
        labs(title="pvalue, gap50, minlen50")
  }


  pdf("reports/rdc_compare_with_published.pdf", width=8.27, height=8.27)
  overlaps_venn1 = overlaps_df %>%
    dplyr::filter(rdc_subset=="Wei+DKFZ") %>%
    dplyr::filter(rdc_is_significant) %>%
    # dplyr::filter(rdc_significant_sumarea>=100e3 & rdc_bootstrap_pvalue<=0.01) %>%
    # dplyr::filter(rdc_maxscore>=2 & rdc_extended_length>=300e3 & rdc_significant_sumarea>=100e3) %>%
    dplyr::filter(rdc_bait_chrom %in% c("chr5", "chr6", "chr8", "Other chromosome")) %>%
    dplyr::filter(grepl("APH-", tlx_group) & grepl("^(Wei2018-NPC|DKFZ|Wei)$", rdc_source) | grepl("DMSO-", tlx_group) & grepl("^(Wei2018_DMSO-NPC|Wei2018-NPC|DKFZ|Wei)$", rdc_source))  %>%
    dplyr::group_split(tlx_group) %>%
    lapply(FUN=function(odf) {
      ooo <<- odf
      odf_filter = odf %>% dplyr::mutate(rdc_source=factor(rdc_source))
      odf_list = split(odf_filter$rdc_reduced_id, odf_filter$rdc_source, drop=F)
      odf_list = odf_list[sort(names(odf_list))]
      if(length(odf_list)<2) {
        return(grid::grob())
      }
      ggvenn::ggvenn(odf_list, fill_color=RColorBrewer::brewer.pal(8, "Pastel2")[1:length(odf_list)], stroke_size=0.5, set_name_size=2, show_percentage=F) +
        ggtitle(gsub("Intra", "Intra (chr5, chr6, chr8)", odf$tlx_group[1])) +
        theme(plot.margin=unit(rep(8, 4), "pt"))
    })
  gridExtra::grid.arrange(grobs=c(overlaps_venn1), ncol=3, top=grid::textGrob("Overlap between newly calculated\nand published RDCs (Wei+DKFZ)", gp=grid::gpar(fontsize=20,font=3)), padding=grid::unit(100, "pt"))

  overlaps_venn2 = overlaps_df %>%
    dplyr::filter(rdc_is_significant) %>%
    # dplyr::filter(rdc_bait_chrom %in% c(bait_chromosomes, "Other chromosome")) %>%
    dplyr::filter(grepl("APH-Inter|DMSO-Inter|APH-Intra|DMSO-Intra", tlx_group) & rdc_source=="DKFZ" & rdc_subset=="Wei+DKFZ") %>%
    dplyr::mutate(tlx_subgroup=paste0(gsub("(APH|DMSO)-", "", tlx_group), "---", gsub("-(Inter|Intra)", "", tlx_group))) %>%
    tidyr::separate_rows(tlx_subgroup, sep="---") %>%
    dplyr::filter(tlx_subgroup %in% c("DMSO", "APH") & rdc_chrom %in% c("chr5", "chr6", "chr8") | tlx_subgroup %in% c("Intra") & rdc_chrom %in% c("chr5", "chr6", "chr8") | tlx_subgroup %in% c("Inter", "Intra")) %>%
    dplyr::mutate(tlx_subgroup=ifelse(tlx_subgroup %in% c("DMSO", "APH", "Intra"), paste(tlx_subgroup, " (chr5, chr6, chr8)"), tlx_subgroup)) %>%
    dplyr::group_split(tlx_subgroup) %>%
    lapply(FUN=function(odf) {
      ooo <<- odf
      odf_filter = odf %>% dplyr::mutate(tlx_group=factor(tlx_group))
      odf_list = split(odf_filter$rdc_reduced_id, odf_filter$tlx_group, drop=F)
      odf_list = odf_list[sort(names(odf_list))]
      ggvenn::ggvenn(odf_list, fill_color=RColorBrewer::brewer.pal(8, "Pastel2")[1:length(odf_list)], stroke_size=0.5, set_name_size=2, show_percentage=F) +
        ggtitle(odf$tlx_subgroup[1]) +
        theme(plot.margin=unit(rep(8, 4), "pt"))
    })

  gridExtra::grid.arrange(grobs=c(overlaps_venn2), ncol=2, top=grid::textGrob("Overlap between DMSO\nand APH RDCs (Wei+DKFZ)", gp=gpar(fontsize=20,font=3)))
  dev.off()
}