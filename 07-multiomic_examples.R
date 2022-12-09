Sys.setenv(TZ='GMT')
library(readr)
library(dplyr)
library(ggplot2)
library(cowplot)
library(stringr)
devtools::load_all("breaktools/")
source("00-utils.R")

multiomics_examples = function()
{
  dir.create("reports/07-multiomic_examples", recursive=T, showWarnings=F)
  params = macs2_params(extsize=1e5, exttype="symmetrical", minqvalue=0.01, effective_size=1.87e9, maxgap=2e5, minlen=1e5)

  #
  # Load offtargets
  #
  offtargets_df = readr::read_tsv("data/offtargets_dkfz.tsv")

  #
  # Load gene annotations
  #
  genes_df = gtf_read('genomes/mm10/annotation/mm10.ncbiRefSeq.gtf.gz')

  #
  # Load RDC
  #
  rdc_df = readr::read_tsv("data/rdc.tsv") %>%
    dplyr::filter(tlx_group %in% c("APH-Inter", "APH-Intra") & rdc_subset=="Wei+DKFZ" & rdc_is_significant) %>%
    dplyr::mutate(rdc_region_start=rdc_extended_start-1e6, rdc_region_end=rdc_extended_end+1e6, rdc_tlx_group=tlx_group) %>%
    dplyr::select(dplyr::matches("rdc_")) %>%
    dplyr::mutate(rdc_no=as.numeric(gsub("RDC-chr[^-]+-", "", rdc_name)), rdc_name_display=paste0(rdc_chrom, ":", rdc_name))

  #
  # Load TLX
  #
  samples_df = tlx_read_paper_samples("data/htgts_samples.tsv", "data") %>%
    dplyr::filter(!control & (
      grepl("(Ctnna2|Nrxn1) promoter/enhancer", experiment) |
      grepl("concentration$", experiment) & concentration==0.4 |
      grepl("Wei|Tena", experiment)))


  tlx_all_df = tlx_read_many(samples_df, threads=16) %>%
    tlx_extract_bait(bait_size=19, bait_region=12e6) %>%
    tlx_mark_dust() %>%
    tlx_mark_offtargets(offtargets_df, offtarget_region=1e5, bait_region=1e4) %>%
    tlx_calc_copynumber(bowtie2_index="genomes/mm10/mm10", max_hits=100, threads=24) %>%
    dplyr::mutate(tlx_group=ifelse(tlx_is_bait_chrom, "Intra", "Inter"))

  #
  # No normalziation here
  #
  libfactors_df = tlx_all_df %>% tlx_libsizes()

  #
  # Filter data that will be displayed
  #
  tlx_clean_df = tlx_all_df %>%
    tlx_remove_rand_chromosomes() %>%
    dplyr::filter(tlx_copynumber==1 & !tlx_duplicated) %>%
    dplyr::filter(!tlx_control & !tlx_is_offtarget)

  #
  # Calculate TLX coverage
  #
  tlxcov_clean_strand_df = tlx_clean_df %>%
    tlx_coverage(group="group", exttype=params$exttype, extsize=params$extsize, libfactors_df=libfactors_df, ignore.strand=F)
  tlxcov_clean_all_df = tlx_clean_df %>%
    tlx_coverage(group="group", exttype=params$exttype, extsize=params$extsize, libfactors_df=libfactors_df, ignore.strand=T)

  #
  # Write TLX coverage
  #
  if(F) {
    tlxcov_write_bedgraph(tlxcov_clean_all_df, path=paste0("reports/bedgaph-", extsize), group="group", ignore.strand=T)
    tlxcov_write_bedgraph(tlxcov_clean_all_df, path=paste0("reports/bedgaph-", extsize), group="group", ignore.strand=F)
  }

  #
  # Expand RDC displayed range to include other neighbouring RDC
  #
  rdc_wide_ranges = rdc_df %>%
    dplyr::rename(rdc_region_chrom="rdc_chrom") %>%
    dplyr::select(-rdc_extended_start, -rdc_extended_end) %>%
    df2ranges(rdc_region_chrom, rdc_region_start, rdc_region_end)
  rdc_narrow_ranges = rdc_df %>%
    dplyr::select(rdc_tlx_group, rdc_chrom, rdc_extended_start, rdc_extended_end) %>%
    df2ranges(rdc_chrom, rdc_extended_start, rdc_extended_end)
  rdc2rdc_df = innerJoinByOverlaps(rdc_wide_ranges, rdc_narrow_ranges) %>%
    filter(rdc_tlx_group.x==rdc_tlx_group.y) %>%
    dplyr::group_by(rdc_tlx_group=rdc_tlx_group.x, rdc_name, rdc_name_display, rdc_no) %>%
    dplyr::mutate(rdc_region_start=min(c(rdc_start, rdc_region_start)), rdc_region_end=max(c(rdc_region_end, rdc_end))) %>%
    dplyr::ungroup()

  rdc2tlxcov_clean_strand_df = tlxcov_clean_strand_df %>%
    df2ranges(tlxcov_chrom, tlxcov_start, tlxcov_end) %>%
    innerJoinByOverlaps(rdc_df %>% df2ranges(rdc_chrom, rdc_region_start, rdc_region_end)) %>%
    dplyr::filter(tlxcov_start>=rdc_region_start & tlxcov_end<=rdc_region_end) %>%
    dplyr::group_by(rdc_tlx_group, tlx_group, rdc_name, rdc_name_display, rdc_no) %>%
    dplyr::mutate(tlxcov_pileup=tlxcov_pileup/max(tlxcov_pileup[tlx_strand %in% c("+", "-")])) %>%
    dplyr::mutate(tlxcov_pileup=pmin(tlxcov_pileup, 1)) %>%
    dplyr::ungroup()

  #
  # Search for offtargets
  #
  rdc2offtargets_df = rdc_df %>%
    df2ranges(rdc_chrom, rdc_region_start, rdc_region_end) %>%
    innerJoinByOverlaps(offtargets_df %>% df2ranges(offtarget_chrom, offtarget_start, offtarget_end)) %>%
    dplyr::filter(rdc_chrom==offtarget_bait_chrom)

  rdc2genes_df = rdc_df %>%
    df2ranges(rdc_chrom, rdc_region_start, rdc_region_end) %>%
    innerJoinByOverlaps(genes_df %>% df2ranges(gene_chrom, gene_start, gene_end)) %>%
    dplyr::select(-dplyr::matches("_ranges\\.")) %>%
    dplyr::filter(gene_cluster_i<=3) %>%
    dplyr::mutate(gene_start=pmax(rdc_region_start, gene_start), gene_end=pmin(rdc_region_end, gene_end))

  rdc2tlx_df = tlx_clean_df %>%
    df2ranges(Rname, Junction, Junction) %>%
    innerJoinByOverlaps(rdc_df %>% df2ranges(rdc_chrom, rdc_region_start, rdc_region_end)) %>%
    dplyr::filter(!tlx_control)

  # Load Repli-SEQ data
  #replication_df = readr::read_tsv("data/replication_reduced_subsets.tsv")
  forks_df = readr::read_tsv("data/replication_forks_NPC.tsv") %>%
    dplyr::mutate(fork_strand=ifelse(fork_direction=="telomeric", "+", "-"))


  rdc2replication_df = forks_df %>%
    df2ranges(fork_chrom, fork_start, fork_end) %>%
    innerJoinByOverlaps(rdc_df %>% df2ranges(rdc_chrom, rdc_region_start, rdc_region_end)) %>%
    dplyr::mutate(
      fork_start=pmax(fork_start, rdc_region_start),
      fork_end=pmin(fork_end, rdc_region_end),
      fork_iz=ifelse(fork_direction=="telomeric", fork_start, fork_end),
      fork_tz=ifelse(fork_direction=="telomeric", fork_end, fork_start))

  repliseq_df = readr::read_tsv("data/repliseq_zhao_bmc2020/repliseq_NPC.tsv")
  rdc2repliseq_df = rdc_df %>%
    df2ranges(rdc_chrom, rdc_region_start, rdc_region_end) %>%
    innerJoinByOverlaps(repliseq_df %>% df2ranges(repliseq_chrom, repliseq_start, repliseq_end))  %>%
    dplyr::select(-dplyr::matches("_ranges\\.")) %>%
    dplyr::group_by(rdc_tlx_group, rdc_name, rdc_name_display, repliseq_start) %>%
    dplyr::mutate(repliseq_value=repliseq_value-min(repliseq_value)+0.01) %>%
    dplyr::mutate(repliseq_value=repliseq_value^1.5) %>%
    dplyr::mutate(repliseq_value=repliseq_value/quantile(repliseq_value, 0.95, na.rm=T), repliseq_value=pmin(repliseq_value, 1)) %>%
    dplyr::ungroup()

  #
  # Plot specific RDC examples
  #
  rdc_samples_df = rdc2tlxcov_clean_strand_df %>%
    dplyr::arrange(rdc_no) %>%
    dplyr::group_by(rdc_tlx_group, rdc_name) %>%
    dplyr::summarize(rdc_length=max(tlxcov_end)-min(tlxcov_start), .groups="keep") %>%
    dplyr::ungroup() %>%
    dplyr::mutate(rdc_no=match(rdc_name, unique(rdc_name)), rdc_repliseq_type=paste0("Other ", cut(1:dplyr::n(), 8))) %>%
    dplyr::select(rdc_tlx_group, rdc_name, rdc_repliseq_type, rdc_length) %>%
    dplyr::group_by(rdc_tlx_group, rdc_repliseq_type) %>%
    dplyr::summarize(rdc_filter=paste(gsub("\\.", "\\\\\\.", rdc_name), collapse="|"), rdc_examples_n=dplyr::n(), rdc_examples_length=sum(rdc_length), .groups="keep") %>%
    dplyr::mutate(rdc_examples_missing=max(rdc_examples_n) - rdc_examples_n, rdc_examples_mislength=max(rdc_examples_length) - rdc_examples_length)

  pdf("reports/07-multiomic_examples/multiomics_examples-new.pdf", width=5*8.27, height=2*11.6)
    for(tg in rdc2tlx_df %>% dplyr::distinct(rdc_tlx_group, tlx_group) %>% split(f=1:nrow(.))) {
      rdc_samples_df.tg = rdc_samples_df %>% dplyr::filter(rdc_tlx_group==tg$rdc_tlx_group)
      plist = lapply(split(rdc_samples_df.tg, f=rdc_samples_df.tg$rdc_repliseq_type), FUN=function(df) {
        dff<<-df
        ggplot_rdc_breaks(
          rdc2rdc_df=rdc2rdc_df %>% dplyr::filter(rdc_tlx_group==tg$rdc_tlx_group),
          rdc2tlxcov_df=rdc2tlxcov_clean_strand_df %>% dplyr::inner_join(tg, by=c("rdc_tlx_group","tlx_group")),
          rdc2tlx_df=rdc2tlx_df %>% dplyr::inner_join(tg, by=c("rdc_tlx_group", "tlx_group")),
          rdc2genes_df=rdc2genes_df %>% dplyr::filter(rdc_tlx_group==tg$rdc_tlx_group),
          rdc2repliseq_df=rdc2repliseq_df %>% dplyr::filter(rdc_tlx_group==tg$rdc_tlx_group),
          # rdc2offtargets_df=rdc2offtargets_df,
          rdc2replication_df=rdc2replication_df %>% dplyr::filter(rdc_tlx_group==tg$rdc_tlx_group),
          rdc_filter=df$rdc_filter) +
            ggtitle(df$rdc_repliseq_type) +
            geom_segment(aes(x=0, xend=rdc_examples_mislength, y=-0.2, yend=-0.2), data=df %>% dplyr::mutate(rdc_name="X", rdc_name_display="X"), size=1, color="#FFFFFF") +
            facet_grid(.~rdc_name_display, scales="free_x", space="free_x") +
            scale_x_continuous(labels=scales::label_number(accuracy=0.1, scale=1e-6, suffix="Mb")) +
            ggpubr::theme_pubclean(base_size=6) +
            theme(axis.text.x=element_text(angle=45, hjust=1), axis.text=element_text(size=8))
      })
      g_title = cowplot::ggdraw() + draw_label(paste0(tg$tlx_group, "-chromosomal translocations examples (RDC source: ", tg$rdc_tlx_group, ")"), fontface='bold', x=0, hjust=0)
      g_legend = get_legend(plist[[1]] + theme(legend.box.margin=margin(0, 0, 0, 12)))
      g_plots = cowplot::plot_grid(plotlist=lapply(plist, function(p) p+theme(legend.position="none")), align="v", axis="t", ncol=1)
      print(cowplot::plot_grid(g_title, g_plots, g_legend, ncol=1, rel_heights=c(1, 19, 1)))
    }
  dev.off()
}