library(readr)
library(dplyr)
library(ggplot2)
library(cowplot)
library(stringr)
devtools::load_all('~/Workspace/breaktools/')


main = function()
{
  #
  # Read TLX files
  #
  samples_df = tlx_read_samples("~/Workspace/Datasets/HTGTS/samples/All_samples.tsv", "~/Workspace/Datasets/HTGTS") %>%
    dplyr::filter(grepl("(Ctnna2|Nrxn1) promoter/enhancer", experiment)) %>%
    dplyr::filter(grepl("NXP010|NXP047", group) | grepl("\\((47/5|18/4|38/3|22|22/5|22/37)\\)", group)) %>%
    dplyr::mutate(group=dplyr::case_when(
      control~"DMSO",
      grepl("NXP010|NXP047", group)~"WT",
      grepl("Allelic deletion", group) ~ gsub("Allelic deletion", "fndr", group),
      grepl("Allelic\\+promoter deletion", group) ~ gsub("Allelic\\+promoter deletion", "fndr-prom", group)
    )) %>%
    # dplyr::mutate(group=dplyr::case_when(
    #   grepl("NXP010|NXP047", group)~"WT",
    #   grepl("Ctnna2", experiment) & grepl("Allelic deletion", group) ~ "fndr (47/5 + 71/3)",
    #   grepl("Ctnna2", experiment) & grepl(".*promoter deletion", group) ~ "fndr-prom (18/4)",
    #   grepl("Csmd1", experiment) & grepl("Allelic deletion", group) ~ "fndr (73/16 + 73/7)",
    #   grepl("Csmd1", experiment) & grepl(".*promoter deletion", group) ~ "fndr-prom (13/12 + 37/8)",
    #   grepl("Nrxn1", experiment) & grepl("Allelic deletion", group) ~ "fndr (22)",
    #   grepl("Nrxn1", experiment) & grepl(".*promoter deletion", group) ~ "fndr-prom (22/1 + 22/39 + 22/37 +22/5)"
    # )) %>%
    dplyr::mutate(group_short=gsub("(fndr|prom) .*$", "\\1", group)) %>%
    dplyr::mutate(sample_number=gsub("^[^0-9]+0*(\\d+).*", "\\1", basename(path))) %>%
    dplyr::mutate(treatment=ifelse(control, "DMSO", "APH"))

  tlx_all_df = tlx_read_many(samples_df) %>% tlx_extract_bait(bait_size=19, bait_region=12e6)
  tlx_df = tlx_all_df %>%
    tlx_remove_rand_chromosomes() %>%
    tlx_calc_copynumber(bowtie2_index="~/Workspace/genomes/mm10/mm10", max_hits=100, threads=24) %>%
    dplyr::filter(tlx_copynumber==1 & !tlx_duplicated) %>%
    dplyr::group_by(run, tlx_sample) %>%
    dplyr::filter(dplyr::n()>1000) %>%
    dplyr::ungroup()

  #
  # Normalization
  #
  tlx_clean_df = tlx_df
    # dplyr::filter(tlx_is_bait_chrom & !tlx_is_bait_junction)
  libfactors_df = tlx_clean_df %>%
    dplyr::mutate(tlx_group=experiment) %>%
    tlx_libfactors(normalize_within="group", normalize_between="none", normalization_target="max") %>%
    dplyr::filter(library_size>=5000)

  #
  # Export bedgraph
  #
  if(F) {
    params = macs2_params(extsize=5e4, exttype="symmetrical")
    tlxcov_df = tlx_df %>%
      dplyr::filter(tlx_is_bait_chrom & !tlx_is_bait_junction) %>%
      dplyr::mutate(tlx_group=group_short) %>%
      tlx_coverage(group="group", extsize=params$extsize, exttype=params$exttype, libfactors_df=libfactors_df, ignore.strand=T)
    tlxcov_df %>%
      tlxcov_write_bedgraph(path="reports/promoter-enhancer_deletion/bedgraph", group="group")
  }

  #
  # Plots
  #
  roi_df = readr::read_tsv("data/promoter-enhancer_deletion_roi.tsv")
  missing_roi_df = tlx_df %>%
    dplyr::distinct(experiment, group, group_short, treatment, tlx_sample, tlx_control) %>%
    dplyr::inner_join(roi_df, by=c("experiment"="roi_experiment")) %>%
    dplyr::select(tlx_sample=tlx_sample, roi_gene=roi_gene, missing_experiment=experiment, missing_group=group, missing_group_short=group_short, missing_treatment=treatment, missing_sample=tlx_sample, missing_control=tlx_control, missing_gene=roi_gene, missing_is_effected=roi_is_effected)
  roi_ranges = roi_df %>% df2ranges(roi_chrom, roi_start, roi_end)
  tlx_ranges = tlx_df %>% df2ranges(Rname, Junction, Junction)
  tlx2roi_df = innerJoinByOverlaps(roi_ranges, tlx_ranges) %>%
    dplyr::filter(roi_experiment==experiment) %>%
    dplyr::group_by(experiment, group, group_short, treatment, tlx_sample, tlx_control, roi_gene, roi_start, roi_end) %>%
    dplyr::summarize(breaks_count=n()) %>%
    dplyr::ungroup() %>%
    dplyr::right_join(missing_roi_df, by=c("tlx_sample", "roi_gene")) %>%
    dplyr::mutate(experiment=missing_experiment, group=missing_group, group_short=missing_group_short, treatment=missing_treatment, tlx_sample=missing_sample, tlx_control=missing_control, roi_gene=missing_gene, roi_is_effected=missing_is_effected, breaks_count=tidyr::replace_na(breaks_count, 0)) %>%
    dplyr::select(-dplyr::matches("^missing_")) %>%
    dplyr::group_by(tlx_sample) %>%
    dplyr::mutate(roi_breaks_count=sum(breaks_count[!roi_is_effected])) %>%
    dplyr::inner_join(libfactors_df %>% dplyr::select(tlx_sample, library_factor), by="tlx_sample") %>%
    dplyr::mutate(breaks_norm_count=(breaks_count)*library_factor) %>%
    dplyr::mutate(group_short=ifelse(tlx_control, "DMSO", as.character(group_short))) %>%
    dplyr::group_by(experiment, roi_gene) %>%
    dplyr::mutate(breaks_norm_mbr=1e6*breaks_norm_count/(roi_end-roi_start)) %>%
    dplyr::mutate(breaks_norm_rel=breaks_norm_count/ifelse(max(breaks_norm_count)==0, 1, max(breaks_norm_count))) %>%
    dplyr::ungroup()  %>%
    dplyr::mutate(group_short=factor(group_short, c("WT", "fndr", "fndr-prom", "DMSO"))) %>%
    dplyr::arrange(group_short) %>%
    dplyr::mutate(group=factor(group, unique(group))) %>%
    dplyr::mutate(roi_gene = factor(roi_gene, rev(roi_df %>% dplyr::arrange(roi_is_effected) %>% .$roi_gene)))



  # pdf("reports/promoter_enhancer_boxplots.pdf", width=8.27, height=11.69, paper="a4")
  tlx2roi_pulled_df = tlx2roi_df %>%
    dplyr::mutate(roi_gene=ifelse(roi_is_effected, "Affected gene", "Other genes pulled")) %>%
    dplyr::group_by(experiment, roi_gene, tlx_sample, group, group_short, treatment) %>%
    dplyr::summarize(breaks_norm_count=sum(breaks_norm_count), breaks_norm_rel=mean(breaks_norm_rel), breaks_norm_mbr=mean(breaks_norm_mbr)) %>%
    dplyr::group_by(experiment, roi_gene) %>%
    dplyr::mutate(breaks_norm_rel=breaks_norm_rel/max(breaks_norm_rel)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(facet=paste0(roi_gene, " (", experiment, ")"), facet=factor(facet, unique(facet)))

  set.seed(123)
  tlx2roi_pulled_stat = tlx2roi_pulled_df %>%
    dplyr::group_by(facet, roi_gene, group) %>%
    dplyr::mutate(n=dplyr::n()) %>%
    dplyr::group_by(facet, roi_gene) %>%
    dplyr::filter(n>=2) %>%
    dplyr::group_by(facet) %>%
    rstatix::t_test(breaks_norm_rel ~ group) %>%
    rstatix::remove_ns(col="p", signif.cutoff=0.05) %>%
    rstatix::add_xy_position(scales="free") %>%
    dplyr::group_by(facet) %>%
    dplyr::mutate(y.position=1+cumsum(rep(0.1, dplyr::n()))) %>%
    dplyr::ungroup()

  fill_pal = dplyr::case_when(
    levels(tlx2roi_pulled_df$group)=="DMSO"~"#999999",
    levels(tlx2roi_pulled_df$group)=="WT"~"#E69F00",
    grepl("fndr-prom", levels(tlx2roi_pulled_df$group))~"#009E73",
    grepl("fndr", levels(tlx2roi_pulled_df$group))~"#56B4E9") %>%
    setNames(., levels(tlx2roi_pulled_df$group))

  ggplot(tlx2roi_pulled_df, aes(x=group, y=breaks_norm_rel)) +
    geom_boxplot(aes(fill=group), outlier.shape=NA, show.legend=F) +
    geom_point(color="#000000", size=2.5, position="jitter", show.legend=F) +
    facet_wrap(~facet, ncol=2, scales="free") +
    ggprism::add_pvalue(tlx2roi_pulled_stat, label="p", tip.length = 0.005) +
    ggpubr::theme_pubclean(base_size=12) +
    labs(y="Relative number of translocations", title="Promoter/enhancer deletion breaks (pulled together)") +
    scale_fill_manual(values=fill_pal) +
      scale_y_continuous(limits=c(0, NA)) +
    theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1)) +
    theme(legend.position="bottom", axis.title.x=element_blank(), axis.ticks.x=element_blank())

  set.seed(123)
  pList = lapply(split(tlx2roi_df, tlx2roi_df$experiment), function(df) {
    # df = tlx2roi_df %>% dplyr::filter(grepl("Nrxn1", experiment))
    highlight_df = df %>% dplyr::filter(roi_is_effected) %>% dplyr::distinct(experiment, roi_gene)
    ggplot(df %>% dplyr::mutate(group=droplevels(group))) +
      geom_boxplot(aes(x=roi_gene, fill=group, y=breaks_norm_mbr), outlier.shape=NA) +
      geom_point(aes(x=roi_gene, fill=group, y=breaks_norm_mbr), color="#000000", size=3, position=position_jitterdodge(jitter.width=0.1), show.legend=F, alpha=0.2) +
      # geom_text(aes(x=roi_gene, fill=group, y=breaks_norm_rel, label=tlx_sample), color="#000000", size=3, position=position_jitterdodge(jitter.width=0.1), show.legend=F, alpha=0.6) +
      geom_rect(aes(xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf), alpha=0.05, fill="#FF0000", color="#00000000", data=highlight_df) +
      facet_wrap(~roi_gene, scales="free", nrow=1) +
      scale_fill_manual(values=fill_pal[as.character(unique(df$group))]) +
      labs(y="Number of translocation per 1Mb", title=df$experiment[1]) +
      scale_y_continuous(limits=c(0, NA)) +
      ggpubr::theme_pubclean(base_size=12) +
      theme(legend.position="bottom", axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
  })
  cowplot::plot_grid(plotlist=pList, ncol=1, align='v')
  dev.off()
}