library(readr)
library(dplyr)
library(ggplot2)
library(stringr)
devtools::load_all('breaktools/')


main = function()
{
  #
  # Read TLX files
  #
  samples_df = tlx_read_samples("data/htgts_samples.tsv", "~/Workspace/Datasets/HTGTS") %>%
    dplyr::filter(grepl("(Ctnna2|Nrxn1) promoter/enhancer", experiment)) %>%
    dplyr::filter(grepl("NXP010|NXP047", group) | grepl("\\((47/5|18/4|38/3|22|22/5|22/37|73/16|13/12|37/8)\\)", group)) %>%
    dplyr::mutate(group=dplyr::case_when(
      control~"DMSO",
      grepl("NXP010|NXP047", group)~"WT",
      grepl("Allelic deletion", group) ~ gsub("Allelic deletion", "fndr", group),
      grepl("Allelic\\+promoter deletion", group) ~ gsub("Allelic\\+promoter deletion", "fndr-prom", group)
    )) %>%
    dplyr::mutate(group_short=gsub("(fndr|prom) .*$", "\\1", group)) %>%
    dplyr::mutate(sample_number=gsub("^[^0-9]+0*(\\d+).*", "\\1", basename(path))) %>%
    dplyr::mutate(treatment=ifelse(control, "DMSO", "APH"))

  tlx_all_df = tlx_read_many(samples_df, threads=10) %>%
    tlx_extract_bait(bait_size=19, bait_region=12e6) %>%
    tlx_calc_copynumber(bowtie2_index="genomes/mm10/mm10", max_hits=100, threads=24)
  tlx_df = tlx_all_df %>%
    tlx_remove_rand_chromosomes()%>%
    dplyr::filter(tlx_copynumber==1 & !tlx_duplicated & Rname!="chrY") %>%
    dplyr::group_by(run, tlx_sample) %>%
    dplyr::filter(dplyr::n()>5000) %>%
    dplyr::ungroup()

  libfactors_df = tlx_df %>%
    tlx_libsizes() %>%
    tlx_libfactors_within(min(library_size)/library_size)

  rdc_df = readr::read_tsv("data/rdc.tsv") %>%
    dplyr::filter(tlx_group=="APH-Inter" & rdc_subset=="Wei+DKFZ" & rdc_is_significant) %>%
    dplyr::select(-tlx_group, -rdc_subset)

  rdc2tlx_df = rdc_df %>%
    df2ranges(rdc_chrom, rdc_extended_start, rdc_extended_end) %>%
    innerJoinByOverlaps(df2ranges(tlx_df, Rname, Junction, Junction)) %>%
    dplyr::mutate(tlx_translocation=ifelse(tlx_is_bait_chrom, "Intra", "Inter")) %>%
    dplyr::mutate(rdc_is_effected=!is.na(rdc_gene) & stringr::str_detect(experiment, rdc_gene)) %>%
    dplyr::group_by(rdc_is_effected, experiment, group, group_short, treatment, tlx_sample, tlx_control, tlx_translocation, rdc_name, rdc_gene, rdc_start, rdc_end, rdc_extended_start, rdc_extended_end) %>%
    dplyr::summarize(breaks_count=dplyr::n()) %>%
    dplyr::ungroup() %>%
    dplyr::inner_join(libfactors_df %>% dplyr::select(tlx_sample, library_factor), by="tlx_sample") %>%
    dplyr::mutate(breaks_norm_count=breaks_count*library_factor) %>%
    dplyr::mutate(group_short=ifelse(tlx_control, "DMSO", as.character(group_short))) %>%
    dplyr::group_by(experiment, rdc_name) %>%
    dplyr::mutate(breaks_norm_mbr=1e6*breaks_norm_count/(rdc_end-rdc_start), breaks_norm_rel=breaks_norm_count/ifelse(max(breaks_norm_count)==0, 1, max(breaks_norm_count))) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(group_short=factor(group_short, c("WT", "fndr", "fndr-prom", "DMSO"))) %>%
    dplyr::arrange(group_short) %>%
    dplyr::mutate(group=factor(group, unique(group))) %>%
    dplyr::group_by(experiment, rdc_name, group) %>%
    dplyr::mutate(breaks_norm_mbr.group_mean=mean(breaks_norm_mbr, na.rm=T)) %>%
    dplyr::group_by(experiment, rdc_name) %>%
    dplyr::filter(max(breaks_norm_mbr.group_mean)>=5 | any(rdc_is_effected)) %>%
    dplyr::ungroup()

  tlx2roi_pulled_df = rdc2tlx_df %>%
    dplyr::mutate(roi_gene=ifelse(rdc_is_effected, "Affected gene", "Other genes pulled")) %>%
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

  fill_pal = with(tlx2roi_pulled_df, dplyr::case_when(
    levels(group)=="DMSO"~"#999999", levels(group)=="WT"~"#E69F00", grepl("fndr-prom", levels(group))~"#009E73", grepl("fndr", levels(group))~"#56B4E9") %>%
    setNames(., levels(group)))


  pdf("reports/05-promoter_enhancer_deletion-boxplots_pulled.pdf", width=8.27, height=11.69, paper="a4")
  ggplot(tlx2roi_pulled_df, aes(x=group, y=breaks_norm_rel)) +
    geom_boxplot(aes(fill=group), outlier.shape=NA, show.legend=F) +
    geom_point(color="#000000", size=2.5, position=position_jitter(width=0.2, height=0, seed=2), show.legend=F) +
    facet_wrap(~facet, ncol=2, scales="free") +
    ggprism::add_pvalue(tlx2roi_pulled_stat, label="p", tip.length = 0.005) +
    ggpubr::theme_pubclean(base_size=14) +
    labs(y="Relative number of translocations", title="Promoter/enhancer deletion breaks (pulled together)") +
    scale_fill_manual(values=fill_pal) +
    scale_y_continuous(limits=c(0, NA)) +
    theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1)) +
    theme(legend.position="bottom", axis.title.x=element_blank(), axis.ticks.x=element_blank())
  dev.off()

  set.seed(123)
  pdf("reports/05-promoter_enhancer_deletion-boxplots.pdf", width=11.69, height=8.27, paper="a4r")
  rdc2tlx_df %>%
    dplyr::filter(!is.na(rdc_gene)) %>%
    dplyr::group_split(experiment) %>%
    lapply(FUN=function(df) {
      dff<<-df
      highlight_df = df %>% dplyr::filter(rdc_is_effected) %>% dplyr::distinct(experiment, rdc_name, rdc_gene)
      ggplot(df %>% dplyr::mutate(group=droplevels(group))) +
        geom_boxplot(aes(x=rdc_gene, fill=group, y=breaks_norm_mbr), outlier.shape=NA) +
        geom_point(aes(x=rdc_gene, fill=group, y=breaks_norm_mbr), color="#000000", size=3, position=position_jitterdodge(jitter.width=0.1), show.legend=F, alpha=0.2) +
        geom_rect(aes(xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf), alpha=0.05, fill="#FF0000", color="#00000000", data=highlight_df) +
        facet_wrap(~paste0(gsub("RDC_", "", rdc_name), " (", rdc_gene, ")"), scales="free") +
        scale_fill_manual(values=fill_pal[as.character(unique(df$group))]) +
        labs(y="Number of translocation per 1Mb", title=df$experiment[1]) +
        scale_y_continuous(limits=c(0, NA)) +
        ggpubr::theme_pubclean(base_size=12) +
        theme(legend.position="bottom", axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), legend.key.size=unit(1.2, 'cm'), strip.text.x = element_text(size=10))
    })
  dev.off()
}