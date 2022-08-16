library(readr)
library(dplyr)
library(ggplot2)
library(stringr)
devtools::load_all('breaktools/')

APH_concentration = function()
{
  debug = F
  dir.create("reports/04-concentration", recursive=T, showWarnings=F)
  group_palette = c("APH 0.2 uM 96h"="#C6DBEF", "APH 0.3 uM 96h"="#6DAACE", "APH 0.4 uM 96h"="#317BA5", "APH 0.6 uM 96h"="#335E9D", "DMSO"="#CCCCCC")
  translocation_palette = c("Inter"="#00D1FF", "Intra"="#E40006")

  #
  # Load offtargets
  #
  offtargets_df = readr::read_tsv("data/offtargets_dkfz.tsv")

  #
  # Repliseq data
  #
  replication_df = readr::read_tsv("data/replication_reduced_subsets.tsv")
  replication_ranges = replication_df %>%
    df2ranges(replication_chrom, pmin(replication_start, replication_end), pmax(replication_start, replication_end))

  #
  # Read TLX
  #
  samples_df = tlx_read_samples(annotation_path="data/htgts_samples.tsv", samples_path="~/Workspace/Datasets/HTGTS") %>%
    dplyr::filter(grepl("APH concentration", experiment))

  tlx_all_df = tlx_read_many(samples_df, threads=24) %>%
    tlx_extract_bait(bait_size=19, bait_region=12e6) %>%
    tlx_calc_copynumber(bowtie2_index="genomes/mm10/mm10", max_hits=100, threads=24) %>%
    tlx_mark_offtargets(offtargets_df, offtarget_region=1e5, bait_region=1e4)

  tlx_concentration_df = tlx_all_df %>%
    tlx_remove_rand_chromosomes() %>%
    dplyr::filter(tlx_copynumber==1 & !tlx_duplicated) %>%
    dplyr::group_by(tlx_sample) %>%
    dplyr::filter(dplyr::n()>5000) %>%
    dplyr::ungroup() %>%
    dplyr::filter(!tlx_is_offtarget) %>%
    dplyr::mutate(tlx_group=paste0(tlx_group, " (", ifelse(tlx_is_bait_chrom, "Intra", "Inter"), ")"))

  params_concentration = macs2_params(extsize=50e3, exttype="symmetrical")
  libfactors_centration_df = tlx_concentration_df %>% tlx_libfactors(normalize_within="none", normalize_between="group", normalization_between_target="samplecount_min")
  # tlx_concentration_df = tlx_clean_df %>%
  #   dplyr::filter(!tlx_is_offtarget & tlx_is_bait_chrom)
  tlxcov_concentration_df = tlx_concentration_df %>%
    tlx_coverage(group="group", extsize=params_concentration$extsize, exttype=params_concentration$exttype, libfactors_df=libfactors_centration_df, ignore.strand=F)
  tlx_concentration_ranges = tlx_concentration_df %>% df2ranges(Rname, Junction, Junction)
  tlxcov_concentration_ranges = tlxcov_concentration_df %>% df2ranges(tlxcov_chrom, tlxcov_start, tlxcov_end)


  #
  # Load RDC and test to filter out only the RDC that are significant in concentration samples
  #
  rdc_df = readr::read_tsv("data/rdc.tsv") %>%
    dplyr::filter(tlx_group=="APH-Inter" & rdc_subset=="Wei+DKFZ" & rdc_is_significant) %>%
    dplyr::select(-tlx_group, -rdc_subset)
  rdc2tlxcov_df = rdc_df %>%
    df2ranges(rdc_chrom, rdc_extended_start, rdc_extended_end) %>%
    innerJoinByOverlaps(tlxcov_concentration_ranges) %>%
    dplyr::mutate(tlx_translocation=gsub(".*(Inter|Intra).*", "\\1", tlx_group), tlx_concentration=gsub("(.*96h|DMSO) .*", "\\1", tlx_group))


  #
  # Bootstrap this specific dataset to further fileter significant RDC
  #
  rdc_bootstrap_df = rdc_df %>%
    tidyr::crossing(rdc2tlxcov_df %>% dplyr::distinct(tlx_group, tlx_translocation, tlx_concentration)) %>%
    dplyr::group_by(tlx_group, tlx_translocation, tlx_concentration) %>%
    dplyr::do((function(z){
      zz<<-z
      bootstrap_data_overlaps(
        evaluated_ranges=z %>% df2ranges(rdc_chrom, rdc_extended_start, rdc_extended_end),
        data_ranges=tlx_concentration_ranges[tlx_concentration_ranges$tlx_group==z$tlx_group[1]],
        genome_tiles_step=10000, genome_tiles_width=50000, n_samples=1000)
    })(.)) %>%
    dplyr::ungroup()
  rdc_bootstrap_sumdf = rdc_df %>%
    dplyr::inner_join(rdc_bootstrap_df, by=c("rdc_chrom"="bootstrap_chrom", "rdc_extended_start"="bootstrap_start", "rdc_extended_end"="bootstrap_end")) %>%
    # tidyr::crossing(tlxcov_concentration_df %>% dplyr::distinct(tlx_group)) %>%
    dplyr::group_by(tlx_group, tlx_translocation, tlx_concentration, rdc_name) %>%
    dplyr::filter(!all(bootstrap_data_count==0)) %>%
    dplyr::do((function(z){
      sg = z$bootstrap_data_count[z$bootstrap_type=="signal"]
      fit_gamma = fitdistrplus::fitdist(z$bootstrap_data_count[z$bootstrap_type=="background"], distr="nbinom", method="qme", probs=c(0.1, 0.9))
      fc = sg/fit_gamma$estimate["mu"]
      pvalue = pnbinom(sg, mu=fit_gamma$estimate["mu"], size=fit_gamma$estimate["size"], lower.tail=F)
      data.frame(concentration_bootstrap_pvalue=pvalue, concentration_bootstrap_fc=fc)
    })(.)) %>%
    dplyr::ungroup()

  #
  # Describe pileup to select pileups with enough data
  #
  rdc_peaks_sumdf = rdc2tlxcov_df %>%
    dplyr::group_by(tlx_group, tlx_translocation, tlx_concentration, rdc_name) %>%
    dplyr::do((function(z){
      z %>%
        dplyr::filter(tlxcov_pileup>=4 & tlxcov_pileup>=max(0.5*tlxcov_pileup)) %>%
        df2ranges(tlxcov_chrom, tlxcov_start, tlxcov_end, tlx_strand) %>%
        GenomicRanges::reduce(min.gapwidth=50e3) %>%
        as.data.frame()
    })(.)) %>%
    dplyr::group_by(tlx_group, tlx_translocation, tlx_concentration, rdc_name) %>%
    dplyr::summarise(peaks_telomeric_count=sum(strand=="+"), peaks_centromeric_count=sum(strand=="-"), peaks_telomeric_length=sum((end-start)[strand=="+"]), peaks_centromeric_length=sum((end-start)[strand=="-"]))

  rdc_concentration_df = rdc_df %>%
    dplyr::select(-dplyr::matches("^(rdc_bootstrap_|peaks_)")) %>%
    inner_join(rdc_bootstrap_sumdf %>% dplyr::select(tlx_group, rdc_name, concentration_bootstrap_pvalue), by="rdc_name") %>%
    inner_join(rdc_peaks_sumdf, by=c("tlx_group", "rdc_name")) %>%
    dplyr::mutate(concentration_is_significant=concentration_bootstrap_pvalue<=0.01) %>%
    dplyr::mutate(concentration_has_single_peak=peaks_telomeric_count==1 & peaks_centromeric_count==1 & pmin(peaks_telomeric_length, peaks_centromeric_length)>=50e3)
  rdc2tlxcov_df = rdc_df %>%
    df2ranges(rdc_chrom, rdc_extended_start-250e3, rdc_extended_end+250e3) %>%
    innerJoinByOverlaps(tlxcov_concentration_ranges) %>%
    dplyr::mutate(tlx_translocation=gsub(".*(Inter|Intra).*", "\\1", tlx_group), tlx_concentration=gsub("(.*96h|DMSO) .*", "\\1", tlx_group))
  rdc2tlxcov_concentration_df = rdc2tlxcov_df %>%
    dplyr::inner_join(rdc_concentration_df %>% dplyr::select(tlx_group, rdc_name, dplyr::starts_with("concentration_")), by=c("tlx_group", "rdc_name")) %>%
    dplyr::filter(concentration_is_significant)

  #
  # Cross-correlation ()
  #
  ccs_df = rdc2tlxcov_concentration_df %>%
    dplyr::filter(tlxcov_pileup>1) %>%
    dplyr::group_by(tlx_group, tlx_concentration, tlx_translocation, rdc_name, rdc_gene_strand, rdc_chrom, rdc_extended_start, rdc_extended_end, concentration_is_significant) %>%
    dplyr::do((function(y){
      yy<<-y
      r = tlx_strand_crosscorrelation(tlxcov_df=y, step=1000, negative_lag=T, negative_correlation=F)
      r$pileup_sense = max(y$tlxcov_pileup[y$tlx_strand=="+"], na.rm=T)
      r$pileup_anti = max(y$tlxcov_pileup[y$tlx_strand=="-"], na.rm=T)
      r
    })(.)) %>%
    dplyr::ungroup()
  ccs_good_df = ccs_df %>%
    dplyr::filter(!is.na(crosscorrelation_lag) & crosscorrelation_cor>0.1 & crosscorrelation_cor_pvalue<=0.01 & pmin(pileup_anti, pileup_sense)>=5)

  #
  # Simple shift in mean
  #
  tlx_shift_df = rdc_df %>%
    df2ranges(rdc_chrom, rdc_start, rdc_end) %>%
    innerJoinByOverlaps(tlx_concentration_ranges) %>%
    dplyr::mutate(tlx_translocation=gsub(".*(Inter|Intra).*", "\\1", tlx_group), tlx_concentration=gsub("(.*96h|DMSO) .*", "\\1", tlx_group)) %>%
    dplyr::arrange(concentration) %>%
    dplyr::mutate(tlx_concentration=factor(tlx_concentration, unique(tlx_concentration))) %>%
    dplyr::ungroup()%>%
    dplyr::group_by(tlx_group, tlx_concentration, tlx_translocation, rdc_name, rdc_gene_strand) %>%
    dplyr::summarize(
      n_sense=sum(tlx_strand=="+"),
      n_anti=sum(tlx_strand=="-"),
      tlx_strand_shift=mean(Junction[tlx_strand=="+"])-mean(Junction[tlx_strand=="-"]),
      tlx_strand_relshift=tlx_strand_shift/(max(rdc_end)-min(rdc_start)),
      tlx_fc=n_sense/(n_sense+n_anti),
      tlx_fc=dplyr::case_when(is.finite(tlx_fc)~tlx_fc, T~NA_real_),
      tlx_relfc=tlx_fc/(max(rdc_end)-min(rdc_start))
    )

  tlx_proportion_df = rdc_df %>%
    df2ranges(rdc_chrom, rdc_start, rdc_end) %>%
    innerJoinByOverlaps(tlx_concentration_ranges) %>%
    dplyr::mutate(tlx_translocation=gsub(".*(Inter|Intra).*", "\\1", tlx_group), tlx_concentration=gsub("(.*96h|DMSO) .*", "\\1", tlx_group)) %>%
    dplyr::arrange(concentration) %>%
    dplyr::mutate(tlx_concentration=factor(tlx_concentration, unique(tlx_concentration))) %>%
    df2ranges(rdc_chrom, Junction, Junction) %>%
    innerJoinByOverlaps(replication_ranges) %>%
    dplyr::mutate(collision=dplyr::case_when(rdc_gene_strand==replication_strand~"co-directional", T~"head-on"))  %>%
    dplyr::ungroup()%>%
    dplyr::group_by(tlx_group, tlx_concentration, tlx_translocation, rdc_name, rdc_gene_strand, collision) %>%
    dplyr::summarize(
      n_sense=sum(tlx_strand=="+"),
      n_anti=sum(tlx_strand=="-"),
      tlx_fc=n_sense/(n_sense+n_anti),
      tlx_relfc=tlx_fc/(max(rdc_end)-min(rdc_start))
    ) %>%
    dplyr::filter(pmin(n_sense, n_anti)>=1)

  pdf("reports/APH_concentration12.pdf", width=11.69, height=8.27, paper="a4r")
  #
  # 1. Plot shift of junctions peak at each APH each concentration
  #
  ggplot(tlx_shift_df %>% dplyr::filter(is.finite(tlx_strand_shift)), aes(y=tlx_strand_relshift, x=tlx_concentration)) +
    geom_hline(yintercept=0, linetype="dashed") +
    geom_boxplot(aes(fill=tlx_concentration), outlier.shape = NA, outlier.alpha = 0) +
    geom_point(aes(size=n_sense+n_anti), position=position_jitter(width=0.2, height=0, seed=2), alpha=0.3) +
    # geom_text(aes(label=gsub("RDC_", "", rdc_name)), color="#FF0000", position=position_jitter(width=0.2, height=0, seed=2)) +
    labs(y="Relative distance between centromeric(+) and telomeric(-)\n junctions average locations", size="No. of breaks", fill="APH concentration") +
    scale_fill_manual(values=group_palette) +
    theme_bw(base_size=12) +
    facet_wrap(~tlx_translocation, scales="free") +
    scale_y_continuous(labels=scales::percent, n.breaks=20) +
    theme(axis.title.x=element_blank(), axis.text.x=element_text(angle=45, vjust=1, hjust=1, size=10), legend.key.size=unit(1.2, 'cm'))

  #
  # 2. Plot proportion of junctions in each concentration
  #
  tlx_shift_df %>%
    dplyr::filter(!is.na(rdc_gene_strand)) %>%
    reshape2::melt(measure.vars=c("n_sense", "n_anti")) %>%
    dplyr::group_by(tlx_translocation, rdc_name) %>%
    dplyr::mutate(junctions_count=(value)/sum(value)) %>%
    dplyr::mutate(variable=ifelse(variable=="n_sense", "+", "-")) %>%
    ggplot(aes(y=junctions_count, x=tlx_concentration)) +
      ggpattern::geom_boxplot_pattern(aes(fill=tlx_concentration, pattern=variable), outlier.shape=NA, outlier.alpha=0, pattern_fill="#000000", pattern_color="#00000000", pattern_spacing=0.01) +
      geom_point(aes(fill=paste0(tlx_concentration, variable)), position=position_jitterdodge(jitter.width=0.2, jitter.height=0, seed=2), alpha=0.3) +
      labs(y="Proportion of RDC breaks per APH concentration (each point is RDC)", size="No. of breaks", fill="APH concentration", pattern="Junction orientation") +
      scale_fill_manual(values=group_palette, guide=guide_legend(override.aes=list(pattern="none"))) +
      ggpattern::scale_pattern_manual(values=c("-"='none', "+"='stripe'), guide=guide_legend(override.aes=list(fill=group_palette[3], colour="#00000000"))) +
      theme_bw(base_size=12) +
      facet_grid(~tlx_translocation) +
      scale_y_continuous(labels=scales::percent, n.breaks=20) +
      theme(axis.title.x=element_blank(), axis.text.x=element_text(angle=45, vjust=1, hjust=1, size=10), legend.key.size=unit(1.2, 'cm'))

  #
  # 3. Plot proportion of junctions in each concentration
  #
  ggplot(tlx_proportion_df %>% dplyr::filter(!is.na(rdc_gene_strand)), aes(y=tlx_fc, x=tlx_concentration)) +
    geom_hline(yintercept=0.5, linetype="dashed") +
    ggpattern::geom_boxplot_pattern(aes(fill=tlx_concentration, pattern=rdc_gene_strand), outlier.shape=NA, outlier.alpha=0, pattern_fill="#000000", pattern_color="#00000000", pattern_spacing=0.01) +
    geom_point(aes(x=tlx_concentration, size=n_sense+n_anti, group=rdc_gene_strand), position=position_jitterdodge(jitter.width=0.1, jitter.height=0), shape=1, alpha=0.8) +
    labs(pattern="Gene strand", y="Proportion of right-moving fork breaks, cent/(cent+tel)\nMost genes will have a head-on area and co-directional area", size="No. of breaks", fill="APH concentration") +
    ggpattern::scale_pattern_manual(values=c("-"='none', "+"='stripe'), guide=guide_legend(override.aes=list(fill=group_palette[3], colour="#00000000"))) +
    scale_fill_manual(values=group_palette, guide=guide_legend(override.aes=list(pattern="none"))) +
    facet_grid(collision~tlx_translocation) +
    scale_y_continuous(labels=scales::percent, n.breaks=10) +
    theme_bw(base_size=12) +
    theme(axis.title.x=element_blank(), axis.text.x=element_text(angle=45, vjust=1, hjust=1, size=10), legend.key.size=unit(1.2, 'cm'))

  #
  # 3. Cross-correlation
  #
  ggplot(ccs_good_df, aes(y=crosscorrelation_lag, x=tlx_concentration)) +
    geom_hline(yintercept=0, linetype="dashed") +
    geom_boxplot(aes(fill=tlx_concentration), outlier.shape = NA, outlier.alpha=0) +
    geom_point(aes(x=tlx_concentration, size=pileup_sense+pileup_anti), color="#FFFFFF", shape=16, position=position_jitter(width=0.2, height=0, seed=2), alpha=0.9) +
    geom_point(aes(x=tlx_concentration, size=pileup_sense+pileup_anti),  shape=21, position=position_jitter(width=0.2, height=0, seed=2), alpha=0.9) +
    labs(y="Relative distance between centromeric(+) and telomeric(-)\nstrand breaks peaks (cross-correlation)", size="Peak height", fill="APH concentration") +
    scale_fill_manual(values=group_palette) +
    scale_color_manual(values=translocation_palette) +
    facet_wrap(~tlx_translocation, scales="free") +
    theme_bw(base_size=12) +
    theme(axis.title.x=element_blank(), axis.text.x=element_text(angle=45, vjust=1, hjust=1, size=10), legend.key.size=unit(1.2, 'cm'))
  dev.off()

  #
  # Export debuging info (offtarget)
  #
  if(debug) {
    tlxcov_concentration_df %>% tlxcov_write_bedgraph(path="reports/04-concentration/bedgraph", group="group", ignore.strand=F)
  }
}