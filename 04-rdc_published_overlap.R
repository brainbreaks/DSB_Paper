library(readr)
library(dplyr)
library(ggplot2)
library(grid)
devtools::load_all("breaktools/")
source("00-utils.R")

rdc_published_overlap = function()
{
  dir.create("reports/03-rdc_published_overlap", recursive=T, showWarnings=F)

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
        dplyr::mutate(rdc_source=paste0(pubrdc_source, "-", pubrdc_celline), rdc_significant_length=1e7, rdc_bootstrap_pvalue=min(rdc_combined_df$rdc_bootstrap_pvalue), rdc_bootstrap_fc=max(rdc_combined_df$rdc_bootstrap_fc), rdc_significant_maxscore=max(rdc_combined_df$rdc_significant_maxscore), rdc_significant_length=max(rdc_combined_df$rdc_significant_length), rdc_length=max(rdc_combined_df$rdc_length), rdc_extended_length=max(rdc_combined_df$rdc_extended_length)) %>%
        tidyr::crossing(data.frame(rdc_subset=unique(rdc_combined_df$rdc_subset))),
    rdc_combined_df %>%
        tidyr::crossing(rdc_source="DKFZ") %>%
        dplyr::select(tlx_group, rdc_chrom, rdc_extended_start, rdc_extended_end, rdc_start, rdc_end, rdc_source, rdc_subset, rdc_is_significant, rdc_bootstrap_pvalue, rdc_bootstrap_fc, rdc_significant_length, rdc_significant_maxscore, rdc_significant_length, rdc_length, rdc_extended_length)) %>%
        dplyr::mutate(rdc_bait_chrom=as.character(ifelse(grepl("Intra", tlx_group), as.character(rdc_chrom), "Other chromosome")))

  #
  # Join RDC that are close to each other (<200kb)
  #
  allrdc_ranges = allrdc_df %>% df2ranges(rdc_chrom, rdc_extended_start, rdc_extended_end)
  allrdc_reduced_df = allrdc_ranges %>%
    GenomicRanges::reduce(min.gapwidth=200e3) %>%
    as.data.frame() %>%
    dplyr::mutate(rdc_reduced_id=paste0(seqnames, ":", start, "-", end)) %>%
    dplyr::select(rdc_reduced_chrom=seqnames, rdc_reduced_start=start, rdc_reduced_end=end, rdc_reduced_id)
  overlaps_df = allrdc_reduced_df %>%
    df2ranges(rdc_reduced_chrom, rdc_reduced_start, rdc_reduced_end) %>%
    innerJoinByOverlaps(allrdc_ranges)

  pdf("reports/03-rdc_published_overlap/wei2018-venn-new.pdf", width=8.27, height=8.27)
  overlaps_venn1 = overlaps_df %>%
    dplyr::filter(rdc_subset %in% c("Wei+DKFZ", "DKFZ")) %>%
    dplyr::filter(rdc_is_significant) %>%
    dplyr::filter(rdc_bait_chrom %in% c("chr5", "chr6", "chr8", "chr17", "Other chromosome")) %>%
    dplyr::filter(grepl("APH-", tlx_group) & grepl("^(Wei2018-NPC|DKFZ|Wei)$", rdc_source) | grepl("DMSO-", tlx_group) & grepl("^(Wei2018_DMSO-NPC|Wei2018-NPC|DKFZ|Wei)$", rdc_source))  %>%
    dplyr::group_split(tlx_group, rdc_subset) %>%
    lapply(FUN=function(odf) {
      ooo <<- odf
      odf_filter = odf %>% dplyr::mutate(rdc_source=factor(rdc_source))
      odf_list = split(odf_filter$rdc_reduced_id, odf_filter$rdc_source, drop=F)
      odf_list = odf_list[sort(names(odf_list))]
      if(length(odf_list)<2) {
        return(grid::grob())
      }

      g = ggvenn::ggvenn(odf_list, fill_color=RColorBrewer::brewer.pal(8, "Pastel2")[1:length(odf_list)], stroke_size=0.5, set_name_size=2, show_percentage=F) +
        ggtitle(gsub("Intra", "Intra (chr5, chr6, chr8, chr17)", odf$tlx_group[1])) +
        theme(plot.margin=unit(rep(8, 4), "pt"))
      g$tlx_group = odf$tlx_group[[1]]
      g$rdc_subset = odf$rdc_subset[[1]]
      g
    })

  gridExtra::grid.arrange(grobs=overlaps_venn1[sapply(overlaps_venn1, function(g) g$rdc_subset)=="Wei+DKFZ"], ncol=3, top=grid::textGrob("Overlap between newly calculated (Wei+DKFZ)\nand published RDCs", gp=grid::gpar(fontsize=20,font=3)), padding=grid::unit(100, "pt"))
  gridExtra::grid.arrange(grobs=overlaps_venn1[sapply(overlaps_venn1, function(g) g$rdc_subset)=="DKFZ"], ncol=3, top=grid::textGrob("Overlap between newly calculated (DKFZ)\nand published RDCs", gp=grid::gpar(fontsize=20,font=3)), padding=grid::unit(100, "pt"))

  overlaps_venn2 = overlaps_df %>%
    dplyr::filter(rdc_source=="DKFZ" & rdc_is_significant) %>%
    dplyr::filter(grepl("APH-Inter|APH-Intra|APH-Combined", tlx_group) & rdc_subset %in% c("Wei+DKFZ", "DKFZ", "Wei")) %>%
    dplyr::group_split(tlx_group) %>%
    lapply(FUN=function(odf) {
      ooo <<- odf
      odf_filter = odf %>% dplyr::mutate(rdc_subset=factor(rdc_subset))
      odf_list = split(odf_filter$rdc_reduced_id, odf_filter$rdc_subset, drop=F)
      odf_list = odf_list[sort(names(odf_list))]
      ggvenn::ggvenn(odf_list, fill_color=RColorBrewer::brewer.pal(8, "Pastel2")[1:length(odf_list)], stroke_size=0.5, set_name_size=2, show_percentage=F) +
        ggtitle(odf$tlx_group[1]) +
        theme(plot.margin=unit(rep(8, 4), "pt"))
    })
  gridExtra::grid.arrange(grobs=overlaps_venn2, ncol=2, top=grid::textGrob("Overlap between newly calculated RDC using different datasets", gp=grid::gpar(fontsize=20,font=3)), padding=grid::unit(100, "pt"))



  overlaps_venn2 = overlaps_df %>%
    dplyr::filter(rdc_is_significant) %>%
    dplyr::filter(grepl("APH-Inter|DMSO-Inter|APH-Intra|DMSO-Intra", tlx_group) & rdc_source=="DKFZ" & rdc_subset=="Wei+DKFZ") %>%
    dplyr::mutate(tlx_subgroup=paste0(gsub("(APH|DMSO)-", "", tlx_group), "---", gsub("-(Inter|Intra)", "", tlx_group))) %>%
    tidyr::separate_rows(tlx_subgroup, sep="---") %>%
    dplyr::filter(tlx_subgroup %in% c("DMSO", "APH") & rdc_chrom %in% c("chr5", "chr6", "chr8", "chr17") | tlx_subgroup %in% c("Intra") & rdc_chrom %in% c("chr5", "chr6", "chr8", "chr17") | tlx_subgroup %in% c("Inter", "Intra")) %>%
    dplyr::mutate(tlx_subgroup=ifelse(tlx_subgroup %in% c("DMSO", "APH", "Intra"), paste(tlx_subgroup, " (chr5, chr6, chr8, chr17)"), tlx_subgroup)) %>%
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