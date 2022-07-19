setwd("~/Workspace/DSB_Paper")

library(readr)
library(dplyr)
library(ggplot2)
library(cowplot)
library(stringr)
library(UpSetR)
library(grid)
devtools::load_all("~/Workspace/breaktools/")



detect_rdc = function()
{
  debug=F
  dir.create("reports/detect_rdc", recursive=T, showWarnings=F)

  #
  # Load offtargets
  #
  offtargets_df = readr::read_tsv("data/offtargets_dkfz.tsv")
  offtargets_ranges = offtargets_df %>%
    df2ranges(offtarget_chrom, offtarget_start, offtarget_end)



  #
  # Load gene annotations
  #
  genes_cache = "tmp/genes.rda"
  if(file.exists(genes_cache)) {
    load(genes_cache)
  } else {}
    genes_df = gtf_read('~/Workspace/genomes/mm10/annotation/mm10.ncbiRefSeq.gtf.gz')
    genes_ranges = genes_df %>%
      dplyr::filter(gene_length>=1e5) %>%
      df2ranges(gene_chrom, gene_start, gene_end)
    save(genes_df, genes_ranges, file=genes_cache)
  }
  #
  # Load RDC
  #
  samples_df = tlx_read_samples("~/Workspace/Datasets/HTGTS/samples/All_samples.tsv", "~/Workspace/Datasets/HTGTS") %>%
    dplyr::filter(tlx_exists & celltype=="NPC" & organism=="mouse" & sample!="VI035" & (
      grepl("(Csmd1|Ctnna2|Nrxn1) promoter/enhancer", experiment) |
      grepl("concentration", experiment) & concentration==0.4 |
      grepl("Wei", experiment))
    )

  tlx_all_df = tlx_read_many(samples_df, threads=10) %>%
    tlx_extract_bait(bait_size=19, bait_region=12e6) %>%
    tlx_calc_copynumber(bowtie2_index="~/Workspace/genomes/mm10/mm10", max_hits=100, threads=24) %>%
    tlx_mark_offtargets(offtargets_df, offtarget_region=1e5, bait_region=1e4)
  libfactors_df = tlx_all_df %>%
    tlx_libsizes()

  #
  # Select junctions suitable for analysis
  #
  tlx_rdc_df = tlx_all_df %>%
    tlx_remove_rand_chromosomes() %>%
    dplyr::filter(tlx_copynumber==1 & !tlx_duplicated & !tlx_is_bait_junction & !tlx_is_offtarget & (!grepl("prom", group) | !tlx_is_bait_chrom)) %>%
    dplyr::mutate(tlx_group=dplyr::case_when(
      !tlx_control & tlx_is_bait_chrom ~ "APH-Intra",
      !tlx_control & !tlx_is_bait_chrom ~ "APH-Inter",
      tlx_control & tlx_is_bait_chrom ~ "DMSO-Intra",
      tlx_control & !tlx_is_bait_chrom ~ "DMSO-Inter",
    ))
  tlx_rdc_df = dplyr::bind_rows(
    tlx_rdc_df %>% dplyr::filter(!tlx_control) %>% dplyr::mutate(tlx_group="APH (Wei+DKFZ)"),
    tlx_rdc_df %>% dplyr::filter(tlx_control) %>% dplyr::mutate(tlx_group="DMSO (Wei+DKFZ)"),
    tlx_rdc_df %>% dplyr::mutate(tlx_group=paste0(tlx_group, " (Wei+DKFZ)")),
    tlx_rdc_df %>% dplyr::mutate(tlx_group=paste0(tlx_group, " (", ifelse(grepl("Wei", experiment), "Wei", "DKFZ"), ")"))
  ) %>% dplyr::mutate(tlx_control=F)


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


  # save(tlx_rdc_df, tlxcov_rdc_df, params_rdc, genes_df, libfactors_df, libfactors_df, offtargets_df, file="backup.rda")
  # load("backup.rda")

  #
  # Detect RDC
  #
  devtools::load_all("~/Workspace/breaktools/")
  # params_rdc = macs2_params(extsize=5e4, exttype="symmetrical", llocal=1e7, minpvalue=0.01, effective_size=1.87e9, maxgap=10e3, minlen=100e3, baseline=2)
  params_rdc = macs2_params(extsize=5e4, exttype="symmetrical", llocal=1e7, minpvalue=0.01, effective_size=1.87e9, maxgap=100e3, minlen=200e3, seedlen=25e3, seedgap=10e3, baseline=2)
  tlxcov_rdc_df = tlx_rdc_df %>%
    tlx_coverage(group="group", extsize=params_rdc$extsize, exttype=params_rdc$exttype, libfactors_df=libfactors_df, ignore.strand=T, recalculate_duplicate_samples=F)
  tlxcov_rdc_strand_df = tlx_rdc_df %>%
    tlx_coverage(group="group", extsize=params_rdc$extsize, exttype=params_rdc$exttype, libfactors_df=libfactors_df, ignore.strand=F, recalculate_duplicate_samples=F)
  macs_rdc = tlxcov_rdc_df %>% tlxcov_macs2(group="group", params=params_rdc, debug_plots=F)
  macs_strand_rdc = tlxcov_rdc_strand_df %>% tlxcov_macs2(group="group", params=params_rdc, debug_plots=F)

  #
  # Write debugging information from RDC calling
  #
  if(debug)
  {
    tlxcov_rdc_strand_df %>% tlxcov_write_bedgraph(path="reports/detect_rdc/bedgraph", group="group")
    tlxcov_rdc_df %>% tlxcov_write_bedgraph(path="reports/detect_rdc/bedgraph", group="group")
    macs_rdc$islands %>%
      dplyr::group_by(tlx_group) %>%
      dplyr::do((function(df){
        df %>%
          dplyr::mutate(score=1) %>%
          dplyr::select(island_chrom, island_start, island_end, island_name, score, island_strand) %>%
          readr::write_tsv(paste0("reports/detect_rdc/islands-", df$tlx_group[1], ".bed"), col_names=F)
      })(.))
    macs_strand_rdc$islands %>%
      dplyr::group_by(tlx_group) %>%
      dplyr::do((function(df){
        df %>%
          dplyr::mutate(score=1) %>%
          dplyr::select(island_chrom, island_start, island_end, island_name, score, island_strand) %>%
          readr::write_tsv(paste0("reports/detect_rdc/islands-strand-", df$tlx_group[1], ".bed"), col_names=F)
      })(.))


    macs_rdc$qvalues %>%
      dplyr::group_by(tlx_group) %>%
      dplyr::do((function(df){
        writeLines('track color="255,102,102" altColor="255,0,0"', con=paste0("reports/detect_rdc/qvalues-", df$tlx_group[1], ".bedgraph"))
        df %>%
          dplyr::select(qvalue_chrom, qvalue_start, qvalue_end, qvalue_score) %>%
          readr::write_tsv(paste0("reports/detect_rdc/qvalues-", df$tlx_group[1], ".bedgraph"), col_names=F, append=T)
        writeLines('track color="255,153,51" altColor="255,0,0"', con=paste0("reports/detect_rdc/qvalue-", df$tlx_group[1], ".bedgraph"))
        df %>%
          dplyr::select(qvalue_chrom, qvalue_start, qvalue_end, bgmodel_signal) %>%
          readr::write_tsv(paste0("reports/detect_rdc/baseline-", df$tlx_group[1], ".bedgraph"), col_names=F, append=T)
      })(.))
    macs_strand_rdc$qvalues %>%
      dplyr::mutate(qvalue_strand_name=ifelse(qvalue_strand=="+", "plus", "minus")) %>%
      dplyr::group_by(tlx_group, qvalue_strand_name) %>%
      dplyr::do((function(df){
        # df = macs_strand_rdc$qvalues %>%
        #   dplyr::filter(tlx_group=="APH-Intra (Wei+DKFZ)" & qvalue_strand=="+")

        writeLines('track color="255,102,102" altColor="255,0,0"', con=paste0("reports/detect_rdc/qvalues-", df$tlx_group[1], "-", df$qvalue_strand_name[1], ".bedgraph"))
        df %>%
          dplyr::select(qvalue_chrom, qvalue_start, qvalue_end, qvalue_score) %>%
          readr::write_tsv(paste0("reports/detect_rdc/qvalues-", df$tlx_group[1], "-", df$qvalue_strand_name[1], ".bedgraph"), col_names=F, append=T)


        writeLines('track color="255,153,51" altColor="255,0,0"', con=paste0("reports/detect_rdc/baseline-", df$tlx_group[1], "-", df$qvalue_strand_name[1], ".bedgraph"))
        df %>%
          dplyr::select(qvalue_chrom, qvalue_start, qvalue_end, bgmodel_signal) %>%
          readr::write_tsv(paste0("reports/detect_rdc/baseline-", df$tlx_group[1], "-", df$qvalue_strand_name[1], ".bedgraph"), col_names=F, append=T)
      })(.))

    offtargets_df %>%
      dplyr::mutate(score=1, strand="*") %>%
      dplyr::select(offtarget_chrom, offtarget_start, offtarget_end, offtarget_bait_chrom, score, strand) %>%
      readr::write_tsv(paste0("reports/detect_rdc/offtargets_dkfz.bed"), col_names=F)


    tlxcov_rdc_df %>%
      tlx_remove_rand_chromosomes() %>%
      dplyr::filter(tlx_copynumber==1 & !tlx_duplicated & !tlx_is_bait_junction & tlx_is_offtarget & (!grepl("prom", group) | !tlx_is_bait_chrom)) %>%
      dplyr::mutate(tlx_group=ifelse(tlx_is_bait_chrom, "Intra", "Inter")) %>%
      tlx_coverage(group="group", extsize=params_rdc$extsize, exttype=params_rdc$exttype, libfactors_df=libfactors_df, ignore.strand=T) %>%
      tlxcov_write_bedgraph(path="reports/detect_rdc/offtargets", group="group")

    tlxcov_rdc_df %>% tlxcov_write_bedgraph(path="reports/detect_rdc/bedgraph", group="group")


    tlx_rdc_df %>% tlx_write_bed(path="reports/detect_rdc/bed", group="group", ignore.strand=T)
    tlx_rdc_df %>%
      tlx_coverage(group="group", extsize=params_rdc$extsize, exttype=params_rdc$exttype, libfactors_df=libfactors_df, ignore.strand=F) %>%
      tlxcov_write_bedgraph(path="reports/detect_rdc/bedgraph", group="group")
  }


# islands_combined_reduced_df %>%
#       dplyr::group_by(island_combined_group) %>%
#       dplyr::do((function(df){
#         df %>%
#           dplyr::mutate(score=1, name="Test", strand="*") %>%
#           dplyr::select(island_combined_chrom, island_combined_start, island_combined_end, name, score, strand) %>%
#           readr::write_tsv(paste0("reports/detect_rdc/islands-combined-", df$island_combined_group[1], ".bed"), col_names=F)
#       })(.))

# islands_combined_df.bck = islands_combined_df
# qvalues_combined_df.bck = qvalues_combined_df
# qvalues_combined_df = rbind(qvalues_combined_df, qvalues_combined_df.bck)
# islands_combined_df = rbind(islands_combined_df, islands_combined_df.bck)

  #
  # Reduce stranded/non-stranded extended peaks to a single set of detected RDC (stratified by tlx_group)
  #
  rdc_maxgap = 500e3
  rdc_minlen = 100e3
  islands_combined_df = dplyr::bind_rows(macs_rdc$islands, macs_strand_rdc$islands)
  qvalues_combined_df = dplyr::bind_rows(macs_rdc$qvalues, macs_strand_rdc$qvalues)
  islands_combined_reduced_df = islands_combined_df %>%
    dplyr::group_by(tlx_group) %>%
    dplyr::do(GenomicRanges::reduce(df2ranges(., island_chrom, island_extended_start, island_extended_end), min.gapwidth=rdc_maxgap) %>% as.data.frame()) %>%
    dplyr::ungroup()  %>%
    dplyr::select(island_combined_group=tlx_group, island_combined_chrom=seqnames, island_combined_start=start, island_combined_end=end)

  #
  # Extract start/end positions from overlapping islands and count islands (stratified by strand)
  #
#   rdc_df = islands_combined_reduced_df %>%
#     df2ranges(island_combined_chrom, island_combined_start, island_combined_end) %>%
#     leftJoinByOverlaps(islands_combined_df %>% df2ranges(island_chrom, island_extended_start, island_extended_end)) %>%
#     dplyr::filter(island_combined_group==tlx_group) %>%
#     dplyr::group_by(tlx_group, island_combined_chrom, island_combined_start, island_combined_end) %>%
#     dplyr::filter(any(island_strand=="*")) %>%
#     dplyr::group_by(tlx_group, rdc_chrom=island_combined_chrom, rdc_extended_start=island_combined_start, rdc_extended_end=island_combined_end) %>%
#     dplyr::summarize(
#       rdc_start=min(island_start),
#       rdc_end=max(island_end),
#       rdc_significant_nplus=sum(island_strand=="+"),
#       rdc_significant_nminus=sum(island_strand=="-"),
#       rdc_significant_ncombined=sum(island_strand=="*"),
#       rdc_significant_strands=paste0(rdc_significant_nplus, "+", rdc_significant_nminus, "-", rdc_significant_ncombined, "*"))
# plot(rdc_df$rdc_start-rdc_df$rdc_extended_start, rdc_df$rdc_extended_end-rdc_df$rdc_end)

rdc_df %>%
  dplyr::rename(rdc_tlx_group="tlx_group") %>%
  dplyr::ungroup() %>%
  dplyr::filter(rdc_start-rdc_extended_start >= 5e6) %>%
  dplyr::slice(1) %>%
  df2ranges(rdc_chrom, rdc_extended_start, rdc_extended_end) %>%
  innerJoinByOverlaps(islands_combined_df %>% df2ranges(island_chrom, island_start, island_end)) %>%
  dplyr::filter(tlx_group==rdc_tlx_group) %>%
  dplyr::select(tlx_group, island_chrom, island_strand, island_start, island_end, island_extended_start, island_extended_end )
# rdc_df %>%
#   dplyr::rename(rdc_tlx_group="tlx_group") %>%
#   dplyr::ungroup() %>%
#   dplyr::filter(rdc_start-rdc_extended_start >= 1e6) %>%
#   dplyr::slice(1) %>%
#   df2ranges(rdc_chrom, rdc_extended_start, rdc_extended_end)%>%
#   innerJoinByOverlaps(qvalues_combined_df %>%df2ranges(qvalue_chrom, qvalue_start, qvalue_end)) %>%
#   dplyr::filter(tlx_group==rdc_tlx_group)%>%
#   dplyr::select(rdc_chrom, rdc_extended_start, rdc_extended_end, dplyr::starts_with("qvalue")) %>%
#   reshape2::melt(measure.vars=c("qvalue_start", "qvalue_end"))  %>%
#   ggplot() +
#     geom_line(aes(x=value, y=qvalue_pvalue, color=paste0(rdc_chrom, ":", rdc_extended_start, "-", rdc_extended_end)))

  #
  # Find overlap between stranded and no-stranded peak detection and reduce
  #
  rdc_df = islands_combined_reduced_df %>%
    df2ranges(island_combined_chrom, island_combined_start, island_combined_end) %>%
    leftJoinByOverlaps(islands_combined_df %>% df2ranges(island_chrom, island_extended_start, island_extended_end)) %>%
    dplyr::filter(island_combined_group==tlx_group) %>%
    dplyr::group_by(tlx_group, island_combined_chrom, island_combined_start, island_combined_end) %>%
    dplyr::filter(any(island_strand=="*")) %>%
    dplyr::group_by(tlx_group, rdc_chrom=island_combined_chrom, rdc_extended_start=island_combined_start, rdc_extended_end=island_combined_end) %>%
    dplyr::summarize(
      rdc_start=min(island_start),
      rdc_end=max(island_end),
      n=dplyr::n(),
      rdc_significant_nplus=sum(island_strand=="+"),
      rdc_significant_nminus=sum(island_strand=="-"),
      rdc_significant_ncombined=sum(island_strand=="*"),
      rdc_significant_strands=paste0(rdc_significant_nplus, "+", rdc_significant_nminus, "-", rdc_significant_ncombined, "*")) %>%
    df2ranges(rdc_chrom, rdc_start, rdc_end) %>%
    leftJoinByOverlaps(qvalues_combined_df %>% dplyr::rename(qvalue_tlx_group="tlx_group") %>% df2ranges(qvalue_chrom, qvalue_start, qvalue_end)) %>%
    dplyr::filter(tlx_group==qvalue_tlx_group)%>%
    dplyr::group_by(tlx_group, rdc_chrom, rdc_start, rdc_end) %>%
    dplyr::do((function(z){
      zz<<-z
      z_reduced = z %>%
        dplyr::filter(qvalue_score>=-log10(params_rdc$minsignif)) %>%
        df2ranges(qvalue_chrom, qvalue_start, qvalue_end) %>%
        GenomicRanges::reduce(min.gapwidth=params_rdc$extsize)
      z %>%
        dplyr::select(dplyr::starts_with("rdc_")) %>%
        dplyr::slice(1) %>%
        dplyr::mutate(
          rdc_significant_sumarea=sum(GenomicRanges::width(z_reduced)),
          rdc_significant_maxarea=max(GenomicRanges::width(z_reduced)),
          rdc_significant_sumprop=rdc_significant_sumarea/(z$rdc_end - z$rdc_start)[1],
          rdc_significant_maxprop=rdc_significant_maxarea/(z$rdc_end - z$rdc_start)[1])
    })(.)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(rdc_subset=gsub(".*\\((.*)\\)", "\\1", tlx_group), tlx_group=gsub(" ?\\(.*", "", tlx_group))
  rdc_df = dplyr::bind_rows(rdc_df,  rdc_df %>% dplyr::filter(grepl("-(Inter|Intra)", tlx_group)) %>% dplyr::mutate(tlx_group=gsub("Inter|Intra", "Combined", tlx_group)))

  rdc_df %>%
    # dplyr::filter(rdc_significant_sumarea>=100000) %>%
    dplyr::group_by(tlx_group, rdc_subset) %>%
    dplyr::do((function(df){
      dff<<-df
      df %>%
        dplyr::mutate(score=1, strand=".") %>%
        dplyr::distinct(rdc_chrom, rdc_start, rdc_end, .keep_all=T) %>%
        dplyr::select(rdc_chrom, rdc_start, rdc_end, score, strand) %>%
        readr::write_tsv(file=paste0("reports/detect_rdc/rdc_", df$tlx_group[1], " (", df$rdc_subset[1], ").bed"), col_names=F)
    })(.))

  # x = x %>%
  #   leftJoinByOverlaps(genes_ranges) %>%
  #   dplyr::arrange(dplyr::desc(tidyr::replace_na(gene_length, 0))) %>%
  #   dplyr::rename(rdc_chrom="island_chrom", rdc_start="island_start", rdc_end="island_end", rdc_extended_start="island_extended_start", rdc_extended_end="island_extended_end", rdc_snr="island_snr", rdc_qvalue_log10="island_qvalue", rdc_pileup="island_summit_abs") %>%
  #   dplyr::mutate(rdc_cluster=paste0("RDC_", stringr::str_pad((0:(dplyr::n()))[-1], 3, pad="0"))) %>%
  #   dplyr::group_by(tlx_group, rdc_cluster, rdc_chrom, rdc_start, rdc_end, rdc_extended_start, rdc_extended_end, rdc_snr, rdc_pileup, rdc_qvalue_log10, island_score_mean, island_score_maxwidth, island_score_allwidth) %>%
  #   dplyr::summarize(
  #     genes_count=sum(!is.na(gene_id)),
  #     rdc_strand=ifelse(length(unique(na.omit(gene_strand)))==1, unique(na.omit(gene_strand)), "."),
  #     rdc_gene=ifelse(all(is.na(gene_id)), "-", paste0(na.omit(unique(gene_id))[pmin(5, length(na.omit(unique(gene_id))))], collapse=",")),
  #     rdc_group=dplyr::case_when(genes_count==1~1, genes_count>1~2, T~3),
  #     rdc_cluster_display=paste0(rdc_chrom, ":", rdc_cluster, rdc_strand)[1]) %>%
  #   dplyr::mutate(rdc_extended_start=rdc_start, rdc_extended_end=rdc_end) %>%
  #   dplyr::select(tlx_group, rdc_chrom, rdc_start, rdc_end, rdc_extended_start, rdc_extended_end, rdc_cluster,	rdc_group, rdc_strand, rdc_gene, rdc_cluster_display, rdc_snr, rdc_qvalue_log10, rdc_pileup, island_score_mean, island_score_maxwidth, island_score_allwidth) %>%
  #   dplyr::ungroup() %>%
  #   dplyr::mutate(rdc_subset=gsub(".*\\((.*)\\)", "\\1", tlx_group), tlx_group=gsub(" ?\\(.*", "", tlx_group))
  # rdc_df = dplyr::bind_rows(rdc_df,  rdc_df %>% dplyr::filter(grepl("-", tlx_group)) %>% dplyr::mutate(tlx_group=gsub("Inter|Intra", "Combined", tlx_group))  )
  # table(rdc_df$tlx_group, rdc_df$rdc_subset)

  # rdc_df %>%
  #   dplyr::filter(grepl("-(Inter|Intra)", tlx_group)) %>%
  #   dplyr::mutate(val=log10(island_score_maxwidth)) %>%
  #   dplyr::group_by(tlx_group, rdc_subset) %>%
  #   dplyr::arrange(val) %>%
  #   dplyr::mutate(i=(1:dplyr::n())) %>%
  #   dplyr::ungroup() %>%
  #   ggplot() +
  #     geom_line(aes(x=i, y=val, color=rdc_subset)) +
  #     facet_wrap(~tlx_group, scales="free")

  #
  # Load public RDC dataset
  #
  # readr::read_tsv("~/Workspace/Datasets/HTGTS/rdc_macs2_mm10.tsv") %>%
  #   dplyr::filter(RDC_050)
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

# Strong
# Tenm3 (Intra), Auts2 (Inter), Grid2 (Intra), Ctnna2 (Intra),
# Ptpn14 (Inter), Wls (Inter, telomere), Cdkn2a (Inter), Lrig1 (Intra), Tead1 (Intra), Col4a1/Col4a2, Nedd4 (Inter), Amotl2 (Inter) Tubd1/Vmp1 (Inter), Rbfox2 (inter), Zfp608 (Inter), Tcf4 (inter), Malat1 (Inter), Mid1 (Inter, telomere)

  pubrdc_raw_df %>%
    dplyr::group_by(tlx_group) %>%
    dplyr::do((function(df){
      dff<<-df
      df %>%
        dplyr::mutate(score=1, strand=".") %>%
        dplyr::distinct(pubrdc_chrom, pubrdc_start, pubrdc_end, .keep_all=T) %>%
        dplyr::select(pubrdc_chrom, pubrdc_start, pubrdc_end, pubrdc_gene, score, strand) %>%
        readr::write_tsv(file=paste0("data/pubrdc_", df$tlx_group[1], ".bed"), col_names=F)
    })(.)) %>%
    dplyr::summarize()


  pubrdc_df = dplyr::bind_rows(
    pubrdc_raw_df,
    pubrdc_raw_df %>% dplyr::mutate(tlx_group=dplyr::case_when(pubrdc_source=="Wei2018_DMSO" ~ "DMSO",          pubrdc_source %in% c("Tena2020", "Wei2018") ~ "APH")),
    pubrdc_raw_df %>% dplyr::mutate(tlx_group=dplyr::case_when(pubrdc_source=="Wei2018_DMSO" ~ "DMSO-Combined", pubrdc_source %in% c("Tena2020", "Wei2018") ~ "APH-Combined")))

  pubrdc_reduced_df = pubrdc_df %>%
    dplyr::group_by(pubrdc_source, pubrdc_celline, tlx_group) %>%
    dplyr::do(GenomicRanges::reduce(df2ranges(., pubrdc_chrom, pubrdc_start, pubrdc_end)) %>% as.data.frame()) %>%
    dplyr::ungroup() %>%
    dplyr::select(tlx_group, pubrdc_source, pubrdc_celline, pubrdc_chrom=seqnames, pubrdc_start=start, pubrdc_end=end)

  #
  # Combine public and our RDC datasets
  #
  allrdc_df = dplyr::bind_rows(
    pubrdc_reduced_df %>%
      dplyr::mutate(rdc_source=paste0(pubrdc_source, "-", pubrdc_celline), rdc_significant_sumarea=1e7, rdc_significant_maxarea=1e7, rdc_significant_sumprop=1, rdc_significant_maxprop=1e7) %>%
      dplyr::select(tlx_group, rdc_chrom=pubrdc_chrom, rdc_start=pubrdc_start, rdc_end=pubrdc_end, rdc_source, rdc_significant_sumarea, rdc_significant_maxarea, rdc_significant_sumprop, rdc_significant_maxprop) %>%
      tidyr::crossing(data.frame(rdc_subset=unique(rdc_df$rdc_subset))),
    rdc_df %>%
        tidyr::crossing(rdc_source="DKFZ") %>%
        dplyr::select(tlx_group, rdc_chrom, rdc_start, rdc_end, rdc_source, rdc_subset, rdc_significant_sumarea, rdc_significant_maxarea, rdc_significant_sumprop, rdc_significant_maxprop)) %>%
    dplyr::mutate(rdc_bait_chrom=as.character(ifelse(grepl("Intra", tlx_group), as.character(rdc_chrom), "Other chromosome")))

  allrdc_ranges = allrdc_df %>% df2ranges(rdc_chrom, rdc_start, rdc_end)
  allrdc_reduced_df = allrdc_ranges %>%
    GenomicRanges::reduce(min.gapwidth=500e3) %>%
    as.data.frame() %>%
    dplyr::mutate(rdc_reduced_id=paste0(seqnames, ":", start, "-", end)) %>%
    dplyr::select(rdc_reduced_chrom=seqnames, rdc_reduced_start=start, rdc_reduced_end=end, rdc_reduced_id)

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

  overlaps_df = allrdc_reduced_df %>%
    df2ranges(rdc_reduced_chrom, rdc_reduced_start, rdc_reduced_end) %>%
    innerJoinByOverlaps(allrdc_ranges) %>%
    df2ranges(rdc_reduced_chrom, rdc_reduced_start, rdc_reduced_end) %>%
    leftJoinByOverlaps(offtargets_ranges) %>%
    dplyr::filter(is.na(offtarget_bait_chrom) | offtarget_bait_chrom!=rdc_bait_chrom)

  # overlaps_venn = overlaps_df %>%
  #   dplyr::filter(rdc_source != "Tena2020-NPC" & rdc_bait_chrom %in% c("chr5", "chr8", "chr17", "Other chromosome")) %>%
  #   dplyr::mutate(rdc_source=factor(rdc_source)) %>%
  #   split(.$tlx_group) %>%
  #   lapply(FUN=function(odf) {
  #     ooo <<- odf
  #     odf_wide = odf %>%
  #       reshape2::dcast(rdc_reduced_id~rdc_source, fun.aggregate=function(x) as.numeric(length(x)>0)) %>%
  #       tibble::column_to_rownames("rdc_reduced_id")
  #     ComplexUpset::upset(odf_wide, colnames(odf_wide), name='genre', wrap=F, height_ratio=0.5) +
  #       ggtitle(odf$tlx_group[1])
  #   })
  # gridExtra::grid.arrange(grobs=overlaps_venn, ncol = 6)

#   x = overlaps_df %>%
#     dplyr::mutate(rdc_size=rdc_end-rdc_start) %>%
#     dplyr::filter(rdc_subset=="Wei+DKFZ" & grepl("APH-Inter|APH-Intra", tlx_group) & rdc_bait_chrom %in% c("chr5", "chr6", "chr8", "chr17", "Other chromosome")) %>%
#     # dplyr::filter(tlx_group+rdc_reduced_id~rdc_subset)
#     # dplyr::filter(rdc_significant_sumarea<200e3)
#     # dplyr::filter("Wei2018-NPC"==rdc_source | rdc_significant_sumarea<=700e3) %>%
#     dplyr::group_by(rdc_reduced_id, tlx_group, rdc_subset) %>%
#     dplyr::filter(sum("Wei2018-NPC"==rdc_source)>0 & sum("DKFZ"==rdc_source)==0) %>%
#     dplyr::select(rdc_reduced_id, rdc_subset, rdc_significant_sumarea, rdc_source, rdc_size) %>%
#     dplyr::ungroup() %>%
#     dplyr::group_by(rdc_reduced_id) %>%
#     dplyr::filter(sum(rdc_source=="Wei2018-NPC" & tlx_group=="APH-Inter")>=1 & sum(rdc_source=="Wei2018-NPC" & tlx_group=="APH-Intra")>=1) %>%
#     dplyr::ungroup() %>%
#     dplyr::group_by(rdc_subset) %>%
#     dplyr::arrange(rdc_size) %>%
#     dplyr::mutate(i=1:dplyr::n()) %>%
#     dplyr::ungroup() %>%
#     data.frame() %>%
#     dplyr::group_by(tlx_group, rdc_reduced_id) %>%
#     dplyr::summarize(rdc_subset=paste0(rdc_subset, collapse=","), rdc_size=mean(rdc_size))
# table(x$rdc_subset, x$rdc_source)
# x %>% dplyr::filter(i==120) %>% data.frame()
# overlaps_df %>% dplyr::filter(rdc_reduced_id=="chrX:119166391-121296390" & rdc_source=="DKFZ" & grepl("Inter|Intra", tlx_group))
# overlaps_df %>% dplyr::filter(rdc_reduced_id=="chr17:8567135-11461715" & rdc_source=="DKFZ" & grepl("Inter|Intra", tlx_group))
# rdc_df %>% dplyr::filter(rdc_chrom=="chr1" & rdc_start>=130e6-10e6 & rdc_end<=132e6+11e6 & grepl("Inter|Intra", tlx_group))
# macs_strand_rdc$islands %>% dplyr::filter(island_chrom=="chr1" & island_start>=130e6-10e6 & island_end<=132e6+11e6 & grepl("Inter|Intra", tlx_group))

#
#
# ggplot(x) +
#   geom_line(aes(x=i, y=rdc_size, color=rdc_subset)) +
#   facet_wrap(~rdc_subset)

# seq(0.1, max(rdc_df$rdc_significant_sumarea), length.out=100)
  allrdc_reduced_df %>%
    dplyr::select(rdc_reduced_chrom, rdc_reduced_start, rdc_reduced_end) %>%
    readr::write_tsv(file="reports/detect_rdc/allrdc_reduced_df.bed", col_names=F)

  overlaps_venn_mean = overlaps_df %>%
    dplyr::filter(rdc_bait_chrom %in% c("chr5", "chr6", "chr8", "chr17", "Other chromosome")) %>%
    # dplyr::filter(grepl("APH-", tlx_group)) %>%
    tidyr::crossing(cur_cutoff=seq(0.1, max(rdc_df$rdc_significant_sumarea), length.out=100)) %>%
    dplyr::group_by(tlx_group, rdc_subset, cur_cutoff) %>%
    dplyr::do((function(odf) {
      ooo <<- odf
      is_aph = grepl("APH", odf$tlx_group[1])
      tlx_group_short = gsub("(APH|DMSO)-?", "", odf$tlx_group[1])
      tlx_group_short = ifelse(tlx_group_short=="", "All", tlx_group_short)

      if(all(odf$rdc_significant_sumarea<odf$cur_cutoff)) {
        d = data.frame(tlx_group_short=tlx_group_short, condition=ifelse(is_aph, "APH", "DMSO"), fpr=0, tpr=0)
      } else {
        x = odf %>%
          dplyr::mutate(tlx_group_short=tlx_group_short) %>%
          dplyr::filter(is_aph & rdc_source=="Wei2018-NPC" | !is_aph & rdc_source=="Wei2018_DMSO-NPC" | rdc_source=="DKFZ") %>%
          dplyr::filter(rdc_significant_sumarea>=cur_cutoff | rdc_source!="DKFZ") %>%
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
    ggplot(overlaps_venn_mean) +
      geom_line(aes(x=cur_cutoff, y=fpr, color=tlx_group_short, linetype="fpr")) +
      geom_line(aes(x=cur_cutoff, y=tpr, color=tlx_group_short, linetype="tpr")) +
      geom_vline(xintercept=0.6) +
      coord_cartesian(xlim=c(0, 1e6), ylim=c(0, 1)) +
      facet_grid(condition~rdc_subset, scales="free") +
      labs(title="pvalue, gap50, minlen50")

  pdf("reports/rdc_compare_with_published2.pdf", width=8.27, height=8.27)
  overlaps_venn = overlaps_df %>%
    dplyr::filter(rdc_significant_sumarea>=250000) %>%
    dplyr::filter(rdc_bait_chrom %in% c("chr5", "chr6", "chr8", "chr17", "Other chromosome")) %>%
    dplyr::filter(grepl("APH-", tlx_group) & grepl("^(Wei2018-NPC|DKFZ|Wei)$", rdc_source) | grepl("DMSO-", tlx_group) & grepl("^(Wei2018_DMSO-NPC|DKFZ|Wei)$", rdc_source)) %>%
    split(.$tlx_group) %>%
    lapply(FUN=function(odf) {
      ooo <<- odf
      odf_filter = odf %>% dplyr::mutate(rdc_source=factor(rdc_source))
      odf_list = split(odf_filter$rdc_reduced_id, odf_filter$rdc_source, drop=F)
      odf_list = odf_list[sort(names(odf_list))]
      ggvenn::ggvenn(odf_list, fill_color=RColorBrewer::brewer.pal(8, "Pastel2")[1:length(odf_list)], stroke_size=0.5, set_name_size=2, show_percentage=F) +
        ggtitle(gsub("Intra", "Intra (chr5, chr6, chr8, chr17)", odf$tlx_group[1]))
    })
  gridExtra::grid.arrange(grobs=overlaps_venn, ncol=2, top=grid::textGrob("Overlap between RDC sources", gp=gpar(fontsize=20,font=3)))


  overlaps_venn = overlaps_df %>%
    # dplyr::filter(rdc_bait_chrom %in% c("chr5", "chr8", "chr17", "Other chromosome")) %>%
    dplyr::filter(grepl("Inter|Intra", tlx_group) & grepl("Wei2018-NPC|DKFZ", rdc_source)) %>%
    split(.$rdc_source) %>%
    lapply(FUN=function(odf) {
      ooo <<- odf
      odf_filter = odf %>% dplyr::mutate(tlx_group=factor(tlx_group))
      odf_list = split(odf_filter$rdc_reduced_id, odf_filter$tlx_group, drop=F)
      odf_list = odf_list[sort(names(odf_list))]
      ggvenn::ggvenn(odf_list, fill_color=RColorBrewer::brewer.pal(8, "Pastel2")[1:length(odf_list)], stroke_size=0.5, set_name_size=2, show_percentage=F) +
        ggtitle(odf$rdc_source[1])
    })
  gridExtra::grid.arrange(grobs=overlaps_venn, ncol=2, top=grid::textGrob("Overlap between intra- and inter-chromosomal RDC", gp=gpar(fontsize=20,font=3)))
  dev.off()

}