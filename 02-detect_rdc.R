library(readr)
library(dplyr)
library(ggplot2)
library(stringr)
library(grid)
devtools::load_all("breaktools/")
source("00-utils.R")



detect_rdc = function()
{
  debug=F
  dir.create("reports/02-detect_rdc", recursive=T, showWarnings=F)

  #
  # Load offtargets
  #
  offtargets_df = readr::read_tsv("data/offtargets_dkfz.tsv")

  #
  # Load gene annotations
  #
  genes_df = gtf_read('genomes/mm10/annotation/mm10.ncbiRefSeq.gtf.gz')
  genes_ranges = genes_df %>%
    dplyr::filter(gene_length>=1e5) %>%
    df2ranges(gene_chrom, gene_start, gene_end)

  #
  # Load samples data
  #
  samples_df = tlx_read_samples(annotation_path="data/htgts_samples.tsv", samples_path="data") %>%
    dplyr::filter(tlx_exists & celltype=="NPC" & organism=="mouse" & sample!="VI035" & (
      grepl("(Csmd1|Ctnna2|Nrxn1) promoter/enhancer", experiment) |
      grepl("concentration", experiment) & concentration %in% c(0, 0.4) |
      grepl("Wei", experiment))
    )

  #
  # Load TLX
  #
  tlx_all_df = tlx_read_many(samples_df, threads=16) %>%
    tlx_extract_bait(bait_size=19, bait_region=12e6) %>%
    tlx_calc_copynumber(bowtie2_index="genomes/mm10/mm10", max_hits=100, threads=8) %>%
    tlx_mark_offtargets(offtargets_df, offtarget_region=1e5, bait_region=1e4)
  libfactors_df = tlx_all_df %>% tlx_libsizes()

  #
  # Select junctions suitable for analysis
  #
  tlx_rdc_all_df = tlx_all_df %>%
    tlx_remove_rand_chromosomes() %>%
    dplyr::filter(tlx_copynumber==1) %>%
    dplyr::filter(!tlx_duplicated & !tlx_is_bait_junction & !tlx_is_offtarget & (!grepl("prom", group) | !tlx_is_bait_chrom)) %>%
    dplyr::mutate(tlx_group=dplyr::case_when(
      !tlx_control & tlx_is_bait_chrom ~ "APH-Intra",
      !tlx_control & !tlx_is_bait_chrom ~ "APH-Inter",
      tlx_control & tlx_is_bait_chrom ~ "DMSO-Intra",
      tlx_control & !tlx_is_bait_chrom ~ "DMSO-Inter",
    ))
  tlx_rdc_all_df = dplyr::bind_rows(
    tlx_rdc_all_df %>% dplyr::mutate(tlx_group=paste0(tlx_group, " (Wei+DKFZ)")),
    tlx_rdc_all_df %>% dplyr::mutate(tlx_group=paste0(tlx_group, " (", ifelse(grepl("Wei", experiment), "Wei", "DKFZ"), ")"))
  ) %>% dplyr::mutate(tlx_control=F)
  tlx_rdc_df = tlx_rdc_all_df

  pdf("reports/02-detect_rdc/qualified_reads.pdf", width=11.69, height=8.27, paper="a4r")
  tlx_all_sumdf = tlx_all_df %>%
    dplyr::mutate(total_reads=dplyr::n()) %>%
    dplyr::mutate(qualified_read=ifelse(grepl("^chr[0-9X]+$", Rname) & tlx_copynumber==1 & !tlx_duplicated & !tlx_is_bait_junction & !tlx_is_offtarget, "Qualified read", "Discarded read"), translocation=ifelse(tlx_is_bait_chrom, "Intra", "Inter")) %>%
    dplyr::group_by(qualified_read, translocation, total_reads) %>%
    dplyr::summarize(group_reads=dplyr::n(), group_prop=group_reads/total_reads[1])
  ggplot(tlx_all_sumdf) +
    geom_bar(aes(x=paste(translocation, qualified_read), y=group_prop, fill=paste(translocation, qualified_read)), stat="identity") +
    scale_y_continuous(labels = scales::label_percent())
  dev.off()

  #
  # Evaluate extsize
  #
  if(F) {
    extsize_df = data.frame(extsize=c(seq(100, 500, 100), seq(1000, 9000, 1000), seq(10000, 90000, 10000), seq(1e5, 3e5, 5e5))) %>%
      dplyr::rowwise() %>%
      dplyr::do((function(z){
        zz<<-z
        params_extsize = macs2_params(extsize=z$extsize, exttype="symmetrical", llocal=1e7, minqvalue=0.01, effective_size=1.87e9, maxgap=2e5, minlen=2e5, baseline=2)
        tlxcov_rdc_df = tlx_rdc_df %>%
          dplyr::filter(!tlx_is_bait_junction & (!grepl("prom", group) | !tlx_is_bait_chrom)) %>%
          tlx_coverage(group="group", extsize=params_extsize$extsize, exttype=params_extsize$exttype, libfactors_df=libfactors_df, ignore.strand=T)
        macs_rdc = tlxcov_macs2(tlxcov_rdc_df, group="group", params=params_extsize)
        macs_rdc$islands %>%
          tidyr::crossing(as.data.frame(z))
      })(.)) %>%
      dplyr::ungroup()
    extsize_sumdf = extsize_df %>%
      dplyr::group_by(extsize, tlx_group) %>%
      dplyr::summarize(count=dplyr::n(), width=mean(island_end-island_start))
    pdf("reports/extsize_selection.pdf", width=11.69, height=8.27, paper="a4r")
    ggplot(extsize_sumdf) +
      geom_line(aes(x=extsize, y=count)) +
      geom_vline(xintercept=5e4, color="#FF0000") +
      geom_smooth(aes(x=extsize, y=count)) +
      labs(y="RDC count") +
      facet_wrap(~tlx_group, scales="free")
    ggplot(extsize_sumdf) +
      geom_line(aes(x=extsize, y=width)) +
      geom_vline(xintercept=5e4, color="#FF0000") +
      geom_smooth(aes(x=extsize, y=width)) +
      labs(y="RDC width") +
      facet_wrap(~tlx_group, scales="free")
    dev.off()
  }
  #
  # Detect RDC. Split coverage into telomeric, centromeric and combined and detect islands for each of the 3 subsets
  #
  rdc_power_df = data.frame()
  for(sample_frac in seq(1, 1, 0.04)) {
    print(sample_frac)
    tlx_rdc_df = tlx_rdc_all_df %>%
      # dplyr::filter(tlx_group %in% c("APH-Inter (DKFZ)")) %>%
      dplyr::group_by(tlx_group) %>%
      dplyr::sample_frac(sample_frac) %>%
      dplyr::ungroup()

    params_rdc = macs2_params(extsize=5e4, exttype="symmetrical", llocal=1e7, minpvalue=0.01, maxgap=100e3, minlen=100e3, seedlen=50e3, seedgap=10e3, baseline=2)
    tlxcov_rdc_df = tlx_rdc_df %>% tlx_coverage(group="group", extsize=params_rdc$extsize, exttype=params_rdc$exttype, libfactors_df=libfactors_df, ignore.strand=T, recalculate_duplicate_samples=F)
    tlxcov_rdc_strand_df = tlx_rdc_df %>% tlx_coverage(group="group", extsize=params_rdc$extsize, exttype=params_rdc$exttype, libfactors_df=libfactors_df, ignore.strand=F, recalculate_duplicate_samples=F)
    tlxcov_rdc_combined_df = dplyr::bind_rows(tlxcov_rdc_strand_df, tlxcov_rdc_df)
    bgmodel_combned_df = tlxcov_rdc_combined_df %>%
      dplyr::group_by(tlx_group) %>%
      dplyr::do((function(z) {
        z.ranges = df2ranges(z, tlxcov_chrom, tlxcov_start, tlxcov_end, tlx_strand)
        z.mask = coverage_find_empty_intervals(coverage_ranges=z.ranges, coverage_column="tlxcov_pileup", minlen=1e3, maxcoverage=0)
        z.bgmodel_df = macs2_coverage_bgmodel(coverage_ranges=z.ranges, distr="nbinom", coverage_column="tlxcov_pileup", mask_ranges=z.mask, debug_plots=F)
        cbind(z[1, "tlx_group", drop=F], z.bgmodel_df)
      })(.))
    macs_combined_rdc = tlxcov_macs2(tlxcov_df=tlxcov_rdc_combined_df, group="group", params=params_rdc, bgmodel_df=bgmodel_combned_df, debug_plots=F, extended_islands=T, extended_islands_dist=1e6, extended_islands_significance=0.1)

    #
    # Find biggest qvalue among all strands for each interval
    #
    qvalues_corrected_df = macs_combined_rdc$qvalues %>%
      dplyr::group_by(tlx_group) %>%
      dplyr::do((function(tlx_g) {
        tlx_g %>%
          df2ranges(qvalue_chrom, qvalue_start, qvalue_end, qvalue_strand) %>%
          coverage_merge_strands(aggregate_fun=max, score_column="qvalue_score") %>%
          as.data.frame() %>%
          dplyr::select(qvalue_chrom=seqnames, qvalue_start=start, qvalue_end=end, qvalue_strand=strand, qvalue_score=score)
      })(.)) %>%
      dplyr::ungroup()

    #
    # Reduce strand-specific and non-specific extended peaks to a single set of detected RDC (stratified by tlx_group)
    #
    rdc_maxgap = 200e3
    rdc_minlen = 100e3
    islands_combined_reduced_df = macs_combined_rdc$islands %>%
      dplyr::group_by(tlx_group) %>%
      dplyr::do(GenomicRanges::reduce(df2ranges(., island_chrom, island_extended_start, island_extended_end), min.gapwidth=rdc_maxgap) %>% as.data.frame()) %>%
      dplyr::ungroup()  %>%
      dplyr::select(island_combined_group=tlx_group, island_combined_chrom=seqnames, island_combined_extended_start=start, island_combined_extended_end=end) %>%
      dplyr::group_by(island_combined_group) %>%
      dplyr::mutate(island_combined_name=paste0("ISLAND_", stringr::str_pad((0:(dplyr::n()))[-1], 3, pad="0"))) %>%
      dplyr::ungroup()


    #
    # Find overlap between stranded and non-stranded peak detection and reduce
    #
    rdc_genes_df = islands_combined_reduced_df %>%
      df2ranges(island_combined_chrom, island_combined_extended_start, island_combined_extended_end) %>%
      leftJoinByOverlaps(macs_combined_rdc$islands %>% df2ranges(island_chrom, island_extended_start, island_extended_end)) %>%
      dplyr::filter(island_combined_group==tlx_group) %>%
      dplyr::group_by(tlx_group, rdc_chrom=island_combined_chrom, rdc_extended_start=island_combined_extended_start, rdc_extended_end=island_combined_extended_end) %>%
      dplyr::summarize(
        rdc_start=min(island_start),
        rdc_end=max(island_end),
        rdc_significant_orientation = dplyr::case_when(
          any(island_strand=="+") & any(island_strand=="-") ~ "bidirectional",
          any(island_strand=="+") ~ "telomeric",
          any(island_strand=="-") ~ "centromeric",
          any(island_strand=="*") ~ "combined",
          T ~ "none"
        ),
        rdc_significant_orientation_count=paste0(sum(island_strand=="+"), "+", sum(island_strand=="-"), "-", sum(island_strand=="*"), "*")) %>%
      df2ranges(rdc_chrom, rdc_extended_start, rdc_extended_end) %>%
      leftJoinByOverlaps(qvalues_corrected_df %>% dplyr::rename(qvalue_tlx_group="tlx_group") %>% df2ranges(qvalue_chrom, qvalue_start, qvalue_end)) %>%
      dplyr::filter(tlx_group==qvalue_tlx_group) %>%
      dplyr::mutate(qvalue_start=pmax(qvalue_start, rdc_extended_start), qvalue_end=pmin(qvalue_end, rdc_extended_end)) %>%
      dplyr::group_by(tlx_group, rdc_chrom, rdc_start, rdc_end) %>%
      dplyr::do((function(z){
        z %>%
          dplyr::filter(qvalue_score>=-log10(params_rdc$minsignif)) %>%
          dplyr::summarize(rdc_significant_maxscore=max(qvalue_score, na.rm=T), rdc_significant_length=sum(qvalue_end-qvalue_start, na.rm=T)-params_rdc$extsize) %>%
          dplyr::bind_cols(z %>% dplyr::select(dplyr::starts_with("rdc_")) %>% dplyr::slice(1))
      })(.)) %>%
      dplyr::filter(rdc_significant_length>10e3) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(rdc_subset=gsub(".*\\((.*)\\)", "\\1", tlx_group), tlx_group=gsub(" ?\\(.*", "", tlx_group)) %>%
      df2ranges(rdc_chrom, rdc_extended_start, rdc_extended_end) %>%
      leftJoinByOverlaps(genes_ranges) %>%
      dplyr::arrange(dplyr::desc(gene_length)) %>%
      dplyr::distinct(tlx_group, rdc_subset, rdc_chrom, rdc_extended_start, rdc_extended_end, .keep_all=T) %>%
      dplyr::mutate(rdc_gene=gene_id, rdc_gene_strand=gene_strand, rdc_gene_overlap=pmin(rdc_extended_end, gene_end)-pmax(rdc_extended_start, gene_start)) %>%
      dplyr::group_by(tlx_group, rdc_subset) %>%
      dplyr::mutate(rdc_name=paste0("RDC-", rdc_chrom, "-", sprintf((rdc_end+rdc_start)/2e6, fmt='%#.1f')), rdc_length=rdc_end-rdc_start, rdc_extended_length=rdc_extended_end-rdc_extended_start) %>%
      dplyr::ungroup() %>%
      dplyr::select(tlx_group, dplyr::starts_with("rdc_"))

    #
    # Bootstrap background to find exact probability and signal-to-noise fold change
    #
    tlx_rdc_ranges = tlx_rdc_df %>% df2ranges(Rname, Junction, Junction)
    rdc_bootstrap_df = rdc_genes_df %>%
      dplyr::group_by(tlx_group, rdc_subset) %>%
      dplyr::do((function(z){
        bootstrap_data_overlaps(
          evaluated_ranges=z %>% df2ranges(rdc_chrom, rdc_extended_start, rdc_extended_end),
          data_ranges=tlx_rdc_ranges[tlx_rdc_ranges$tlx_group==paste0(z$tlx_group[1], " (", z$rdc_subset[1], ")")],
          genome_tiles_step=10000,
          genome_tiles_width=50000,
          n_samples=1000)
      })(.)) %>%
      dplyr::ungroup()

    rdc_bootstrap_sumdf = rdc_genes_df %>%
      dplyr::inner_join(rdc_bootstrap_df %>% dplyr::select(tlx_group, rdc_subset, dplyr::matches("bootstrap_")), by=c("tlx_group", "rdc_subset", "rdc_chrom"="bootstrap_chrom", "rdc_extended_start"="bootstrap_start", "rdc_extended_end"="bootstrap_end")) %>%
      dplyr::group_by(tlx_group, rdc_subset, rdc_name) %>%
      dplyr::do((function(z){
        sg = z$bootstrap_data_count[z$bootstrap_type=="signal"]
        fit_gamma = fitdistrplus::fitdist(z$bootstrap_data_count[z$bootstrap_type=="background"], distr="nbinom", method="qme", probs=c(0.1, 0.9))
        fc = sg/fit_gamma$estimate["mu"]
        pvalue = pnbinom(sg, mu=fit_gamma$estimate["mu"], size=fit_gamma$estimate["size"], lower.tail=F)
        data.frame(rdc_bootstrap_pvalue=pvalue, rdc_bootstrap_fc=fc)
      })(.)) %>%
      dplyr::ungroup()

    rdc_df = rdc_genes_df %>%
      dplyr::inner_join(rdc_bootstrap_sumdf, by=c("tlx_group", "rdc_subset", "rdc_name")) %>%
      dplyr::mutate(rdc_is_significant=rdc_significant_length>=100e3 & rdc_bootstrap_pvalue<=0.01) %>%
      dplyr::arrange(tlx_group, rdc_chrom, rdc_start)

    rdc_sample_df = rdc_df %>%
      dplyr::mutate(tlx_orig_group=paste0(tlx_group, " (", rdc_subset, ")")) %>%
      dplyr::inner_join(tlx_rdc_all_df %>% dplyr::group_by(tlx_group) %>% dplyr::summarise(tlx_group_n=dplyr::n()), by=c("tlx_orig_group"="tlx_group")) %>%
      dplyr::mutate(sample_n=sample_frac*tlx_group_n, sample_frac=sample_frac, extsize=params_rdc$extsize)
    rdc_power_df = rbind(rdc_power_df, rdc_sample_df)
  }
  #readr::write_tsv(rdc_subset_df, file="data/rdc_power.tsv")


  #
  # Add manual RDC group annotation
  #
  rdc_groups_df = readr::read_tsv("data/rdc_groups.tsv")
  rdc_df = rdc_df %>%
    df2ranges(rdc_chrom, rdc_extended_start, rdc_extended_end) %>%
    leftJoinByOverlaps(rdc_groups_df %>%
      dplyr::select(group_chrom=rdc_chrom, group_translocation=tlx_translocation, group_start=rdc_extended_start, group_end=rdc_extended_end, rdc_group, rdc_peak_shape, rdc_is_biphasic, rdc_DRIP_peaks, rdc_forks, rdc_sizesize) %>%
      df2ranges(group_chrom, group_start, group_end)) %>%
    group_by_at(colnames(rdc_df)) %>%
    # dplyr::filter(rdc_name=="RDC-chr14-122.4") %>%
    # group_by(rdc_name) %>%
    dplyr::summarise(
      rdc_group=dplyr::case_when(
        any(group_translocation==tlx_group & !is.na(rdc_group)) ~ rdc_group[group_translocation==tlx_group & !is.na(rdc_group)][1],
        any(!is.na(rdc_group)) ~ na.omit(rdc_group)[1]),
      rdc_peak_shape=dplyr::case_when(
        any(group_translocation==tlx_group & !is.na(rdc_peak_shape)) ~ rdc_peak_shape[group_translocation==tlx_group & !is.na(rdc_peak_shape)][1],
        any(!is.na(rdc_peak_shape)) ~ na.omit(rdc_peak_shape)[1]),
      rdc_is_biphasic=dplyr::case_when(
        any(group_translocation==tlx_group & !is.na(rdc_is_biphasic)) ~ rdc_is_biphasic[group_translocation==tlx_group & !is.na(rdc_is_biphasic)][1]=="1",
        any(!is.na(rdc_is_biphasic)) ~ na.omit(rdc_is_biphasic)[1]=="1"),
      rdc_DRIP_peaks=dplyr::case_when(
        any(group_translocation==tlx_group & !is.na(rdc_DRIP_peaks)) ~ rdc_DRIP_peaks[group_translocation==tlx_group & !is.na(rdc_DRIP_peaks)][1],
        any(!is.na(rdc_DRIP_peaks)) ~ na.omit(rdc_DRIP_peaks)[1]),
      rdc_forks=dplyr::case_when(
        any(group_translocation==tlx_group & !is.na(rdc_forks)) ~ rdc_forks[group_translocation==tlx_group & !is.na(rdc_forks)][1],
        any(!is.na(rdc_forks)) ~ na.omit(rdc_forks)[1]),
      rdc_sizesize=dplyr::case_when(
        any(group_translocation==tlx_group & !is.na(rdc_sizesize)) ~ rdc_sizesize[group_translocation==tlx_group & !is.na(rdc_sizesize)][1],
        any(!is.na(rdc_sizesize)) ~ na.omit(rdc_sizesize)[1]))

  #
  # Export RDC
  #
  readr::write_tsv(rdc_df, file="data/rdc.tsv")
  rdc_df %>%
    dplyr::filter(rdc_is_significant) %>%
    dplyr::group_by(tlx_group, rdc_subset) %>%
    dplyr::do((function(df){
      dff<<-df
      df %>%
        dplyr::mutate(rdc_display_name=paste0(rdc_gene, " (", rdc_name, ")"), rdc_score=-log10(rdc_bootstrap_pvalue)) %>%
        dplyr::select(rdc_chrom, rdc_extended_start, rdc_extended_end, rdc_display_name, rdc_score, rdc_gene_strand, rdc_start, rdc_end) %>%
        readr::write_tsv(paste0("reports/02-detect_rdc/rdc-", df$tlx_group[1], "_", df$rdc_subset[1], ".bed"), col_names=F)
    })(.))

  #
  # Write debugging information from RDC calling
  #
  if(debug)
  {
    tlx_rdc_df %>% tlx_write_bed(path="reports/02-detect_rdc/bed", group="group", ignore.strand=T)
    tlxcov_rdc_combined_df %>% tlxcov_write_bedgraph(path="reports/02-detect_rdc/bedgraph", group="group", ignore.strand=F)

    macs_combined_rdc$islands %>%
      dplyr::group_by(tlx_group, island_strand) %>%
      dplyr::do((function(df){
        df %>%
          dplyr::select(island_chrom, island_extended_start, island_extended_end, island_name, island_score, island_strand, island_start, island_end) %>%
          readr::write_tsv(paste0("reports/02-detect_rdc/islands-", df$tlx_group[1], "_", df$island_strand[1], ".bed"), col_names=F)
      })(.))

    islands_combined_reduced_df %>%
      dplyr::group_by(island_combined_group) %>%
      dplyr::do((function(df){
        df %>%
          dplyr::mutate(l=island_combined_extended_end-island_combined_extended_start, thickStart=island_combined_extended_start, thickEnd=island_combined_extended_end) %>%
          dplyr::mutate(i=tidyr::replace_na(findInterval(l, seq(rdc_minlen, max(l), length.out=1000))+1, 1), score=1, strand="*", rgb=dplyr::case_when(l<=rdc_minlen~"#BBBBFF", T~colorRampPalette(c("#0000B2", "#FA00FF"))(1000)[i])) %>%
          dplyr::select(island_combined_chrom, island_combined_extended_start, island_combined_extended_end, island_combined_name, score, strand, thickStart, thickEnd, rgb) %>%
          readr::write_tsv(paste0("reports/02-detect_rdc/islands-reduced-", df$island_combined_group[1], ".bed"), col_names=F)
      })(.))

    rdc_df %>%
      dplyr::filter(rdc_is_significant) %>%
      dplyr::group_by(tlx_group, rdc_subset) %>%
      dplyr::do((function(df){
        dff<<-df
        df %>%
          dplyr::mutate(rdc_strand="*") %>%
          dplyr::select(rdc_chrom, rdc_extended_start, rdc_extended_end, rdc_name, rdc_significant_length, rdc_strand, rdc_start, rdc_end) %>%
          readr::write_tsv(paste0("reports/02-detect_rdc/rdc-", df$tlx_group[1], " (", df$rdc_subset[1], ").bed"), col_names=F)
      })(.))

    macs_combined_rdc$qvalues %>%
      dplyr::group_by(tlx_group, qvalue_strand) %>%
      dplyr::do((function(df){
        writeLines('track color="255,102,102" altColor="255,0,0"', con=paste0("reports/02-detect_rdc/qvalues-", df$tlx_group[1], "_", df$qvalue_strand[1], ".bedgraph"))
        df %>%
          dplyr::select(qvalue_chrom, qvalue_start, qvalue_end, qvalue_score) %>%
          readr::write_tsv(paste0("reports/02-detect_rdc/qvalues-", df$tlx_group[1], "_", df$qvalue_strand[1], ".bedgraph"), col_names=F, append=T)
        writeLines('track color="255,153,51" altColor="255,0,0"', con=paste0("reports/02-detect_rdc/qvalue-", df$tlx_group[1], ".bedgraph"))
        df %>%
          dplyr::select(qvalue_chrom, qvalue_start, qvalue_end, bgmodel_signal) %>%
          readr::write_tsv(paste0("reports/02-detect_rdc/baseline-", df$tlx_group[1], "_", df$qvalue_strand[1], ".bedgraph"), col_names=F, append=T)
      })(.))
  }
}

estimate_power = function()
{
  dir.create("reports/02-detect_rdc", recursive=T, showWarnings=F)

  rdc_power_df = readr::read_tsv(file="data/rdc_power.tsv")
  rdc_power_sumdf = rdc_power_df %>%
    dplyr::filter(rdc_is_significant) %>%
    df2ranges(rdc_chrom, rdc_extended_start, rdc_extended_end) %>%
    leftJoinByOverlaps(
      rdc_subset_df %>%
        dplyr::filter(extsize==50e3 & sample_frac==1 & rdc_is_significant) %>%
        dplyr::select(ref_chrom=rdc_chrom, ref_extended_start=rdc_extended_start, ref_extended_end=rdc_extended_end) %>%
        df2ranges(ref_chrom, ref_extended_start, ref_extended_end)
    ) %>%
    dplyr::group_by(tlx_group, extsize, sample_frac, sample_n) %>%
    dplyr::summarize(n_significant=dplyr::n(), n_overlap=sum(!is.na(ref_chrom)))

  pdf("reports/02-detect_rdc/rdc_power_estimation.pdf", width=11.69, height=8.27, paper="a4r")
  ggplot(rdc_power_sumdf) +
    # geom_smooth(aes(x=sample_n, y=n_overlap, color=factor(as.integer(extsize))), se=F) +
    geom_line(aes(x=sample_n, y=n_overlap, color=factor(as.integer(extsize))), se=F) +
    labs(x="Number of reads", y="Number of predicted RDC\noverlapping with final RDC list", color="Ext. size") +
    scale_x_continuous(labels=scales::label_number(accuracy=1, scale=1e-3, suffix="K"))
  dev.off()
}