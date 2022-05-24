library(readr)
library(dplyr)
library(ggplot2)
library(cowplot)
library(stringr)
devtools::load_all('~/Workspace/breaktools/')


main = function()
{
  color_scheme = c(
    "CGNPs;Atoh1-EGFP"="#666666",
    "DMSO (ESC-NPC;NXP010)"="#A6CEE3",
    "ESC-NPC;NXP010"="#1F78B4",
    "DMSO (ESC-NPC;NXP010;18/4)"="#C2BCCA",
    "ESC-NPC;NXP010;18/4"="#605881",
    "DMSO (ESC-NPC;NXP010;38/3)"="#DEABB1",
    "ESC-NPC;NXP010;38/3"="#A1394E",
    "DMSO (ESC-NPC;NXP010;71/3)"="#FB9A99",
    "ESC-NPC;NXP010;71/3"="#E31A1C"
  )

  #
  # Read TLX files
  #
  samples_df = tlx_read_samples("~/Workspace/Datasets/HTGTS/samples/All_samples.tsv", "~/Workspace/Datasets/HTGTS") %>%
    dplyr::filter(grepl("promoter/enhancer", experiment)) %>%
    dplyr::mutate(group=dplyr::case_when(
      group=="Perental cell (NXP010)"~"NXP010",
      group=="Perental cell (NXP047)"~"NXP047",
      grepl("Ctnna2", experiment) & grepl("Allelic deletion", group) ~ "fndr (47/5 + 71/3)",
      grepl("Ctnna2", experiment) & grepl(".*promoter deletion", group) ~ "fndr-prom (18/4)",
      grepl("Csmd1", experiment) & grepl("Allelic deletion", group) ~ "fndr (73/16 + 73/7)",
      grepl("Csmd1", experiment) & grepl(".*promoter deletion", group) ~ "fndr-prom (13/12 + 37/8)",
      grepl("Nrxn1", experiment) & grepl("Allelic deletion", group) ~ "fndr (22)",
      grepl("Nrxn1", experiment) & grepl(".*promoter deletion", group) ~ "fndr-prom (22/1 + 22/39 + 22/37 +22/5)"
    )) %>%
    dplyr::mutate(group_short=gsub("(fndr|prom) .*$", "\\1", group)) %>%
    dplyr::mutate(sample_number=gsub("^[^0-9]+0*(\\d+).*", "\\1", basename(path))) %>%
    dplyr::mutate(treatment=ifelse(control, "DMSO", "APH"))

  tlx_all_df = tlx_read_many(samples_df)
  tlx_df = tlx_all_df %>%
    tlx_extract_bait(bait_size=19, bait_region=12e6) %>%
    tlx_remove_rand_chromosomes() %>%
    tlx_mark_dust() %>%
    tlx_calc_copynumber(bowtie2_index="~/Workspace/genomes/mm10/mm10", max_hits=100, threads=24) %>%
    dplyr::filter(tlx_copynumber==1 & !tlx_duplicated) %>%
    dplyr::group_by(run, tlx_sample) %>%
    dplyr::filter(dplyr::n()>3000) %>%
    dplyr::ungroup()

  #
  # Normalization
  #
  # effective_size = 1.87e9
  # sgRNA_length = 20
  # extsize = 1e5
  # maxgap = extsize*2
  # exttype = "symetrical"
  # threshold_qvalue = 1e-2
  # threshold_pileup = 1
  # slocal = 1e7
  # llocal = 1e7
  # bait_region=6e6
  params = macs2_params(extsize=100e3, exttype="symmetrical")
  libfactors_df = tlx_libfactors(tlx_df, normalize_within="group", normalize_between="none", normalization_target="min")

  roi_df = readr::read_tsv("data/promoter-enhancer_deletion_roi.tsv")
  roi_ranges = roi_df %>% df2ranges(roi_chrom, roi_start, roi_end)
  tlx_ranges = tlx_df %>% df2ranges(Rname, Junction, Junction)
  tlx2roi_df = innerJoinByOverlaps(roi_ranges, tlx_ranges) %>%
    dplyr::group_by(experiment, group, group_short, treatment, tlx_sample, roi_gene) %>%
    dplyr::summarize(breaks_count=n()) %>%
    dplyr::inner_join(libfactors_df, by="tlx_sample") %>%
    dplyr::mutate(breaks_norm_count=breaks_count*library_factor) %>%
    dplyr::group_by(experiment, group, roi_gene) %>%
    dplyr::mutate(breaks_norm=breaks_norm_count/max(breaks_norm_count)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(group_short=factor(group_short, unique(samples_df$group_short))) %>%
    dplyr::arrange(group_short)

  pdf("reports/difference_in_breaks.pdf", width=11.69, height=8.27, paper="a4r")
  set.seed(123)
  p1 = ggplot(tlx2roi_df %>% dplyr::filter(grepl("Ccser1|Ctnna2|Grid2", roi_gene)), aes(x=group_short, y=breaks_norm_count, fill=treatment, group=paste(group_short, treatment))) +
    geom_boxplot(outlier.shape=NA) +
    geom_point(color="#000000", size=4, position=position_jitterdodge(jitter.width=0.1), show.legend=F) +
    facet_wrap(~experiment, scales="free")
  set.seed(123)
  p1 = p1 +
    # geom_text(aes(label=sample_number), color="#FFFFFF", size=2.5, position=position_jitterdodge(jitter.width=0.1), show.legend=F) +
    geom_point(color="#000000", size=2.5, position=position_jitterdodge(jitter.width=0), show.legend=F) +
    # geom_point(aes(color=batch), size=2, position=position_jitterdodge(jitter.width=0), show.legend=T) +
    facet_grid(experiment~roi_gene, scales="free") +
    labs(x="", y="Junctions", fill="Treatment") +
    ggpubr::theme_pubclean(base_size=12) +
    theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1))
    # ggpubr::stat_compare_means(aes(label = ..p.format..), method="t.test", hide.ns=T, bracket.size=1)
  p1

  p2 = ggplot(tlx2roi_df, aes(x=group_short, y=breaks_norm_count, fill=treatment, group=paste(group_short, treatment))) +
    geom_boxplot(outlier.shape=NA) +
    geom_point(color="#000000", size=2.5, position=position_jitterdodge(jitter.width=0.1), show.legend=F) +
    # geom_point(color="#000000", size=2.5, position=position_jitterdodge(jitter.width=0), show.legend=F) +
    # geom_point(aes(color=batch), size=2, position=position_jitterdodge(jitter.width=0), show.legend=T) +
    facet_wrap(~roi_gene, scales="free") +
    labs(x="", y="Junctions", fill="Treatment") +
    ggpubr::theme_pubclean(base_size=12) +
    theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1))
    # ggpubr::stat_compare_means(aes(label = ..p.format..), method="t.test", hide.ns=T, bracket.size=1)
  p2

  tlx2roi_pulled_df = tlx2roi_df %>%
    dplyr::mutate(roi_gene=ifelse(grepl("Ctnna2", roi_gene), "Ctnna2", "Other genes pulled")) %>%
    dplyr::group_by(roi_gene, tlx_sample, group_short, treatment) %>%
    dplyr::summarize(breaks_norm_count=sum(breaks_norm_count))
  p3 = ggplot(tlx2roi_pulled_df, aes(x=group_short, y=breaks_norm_count, fill=treatment, group=paste(group_short, treatment))) +
    geom_boxplot(outlier.shape=NA) +
    geom_point(color="#000000", size=2.5, position=position_jitterdodge(jitter.width=0), show.legend=F) +    # coord_flip() +
    geom_point(size=2, position=position_jitterdodge(jitter.width=0), show.legend=F) +    # coord_flip() +
    facet_wrap(~roi_gene, scales="free") +
    ggpubr::stat_compare_means(aes(label = ..p.format..), method="t.test", hide.ns=T, bracket.size=1)
  p3

  # gridExtra::grid.arrange(p1, p2, p3, ncol=2)
  dev.off()
}