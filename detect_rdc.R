setwd("~/Workspace/DSB_Paper")

library(readr)
library(dplyr)
library(ggplot2)
library(cowplot)
library(stringr)
library(randomcoloR)
devtools::load_all("~/Workspace/breaktools/")



detect_rdc = function()
{
  debug=T
  dir.create("reports/detect_rdc/offtargets", recursive=T)
  # params = macs2_params(extsize=1e5, exttype="symmetrical", llocal=1e7, minqvalue=0.01, effective_size=1.87e9, maxgap=2e5, minlen=1e5)
  #
  #
  # baits_df = readr::read_tsv("~/Workspace/Datasets/HTGTS/dkfz_baits.tsv")

  #
  # Load gene annotations
  #
  genes_df = gtf_read('~/Workspace/genomes/mm10/annotation/mm10.ncbiRefSeq.gtf.gz')
  genes_ranges = genes_df %>%
    dplyr::filter(gene_length>=1e5) %>%
    df2ranges(gene_chrom, gene_start, gene_end)

  #
  # Load RDC
  #
  samples_df = tlx_read_samples("~/Workspace/Datasets/HTGTS/samples/All_samples.tsv", "~/Workspace/Datasets/HTGTS") %>%
    dplyr::filter(tlx_exists & celltype=="NPC" & organism=="mouse" & sample!="VI035" & (
      grepl("(Csmd1|Ctnna2|Nrxn1) promoter/enhancer", experiment) |
      grepl("concentration", experiment) & concentration==0.4 |
      grepl("Wei|Tena", experiment))
    )

  tlx_all_df = tlx_read_many(samples_df, threads=10) %>%
    tlx_extract_bait(bait_size=19, bait_region=12e6) %>%
    tlx_calc_copynumber(bowtie2_index="~/Workspace/genomes/mm10/mm10", max_hits=100, threads=24)
    # dplyr::group_by(tlx_sample) %>%
    # dplyr::filter(dplyr::n()>5000) %>%
    # dplyr::ungroup()
  
  libfactors_df = tlx_all_df %>% tlx_libsizes()
  
  #
  # Detect offtargets
  #
  offtargets_params = macs2_params(extsize=20, exttype="opposite", llocal=1e7, minqvalue=0.01, baseline=2, effective_size=1.87e9, maxgap=10e3, minlen=4)
  tlx_offtarget_df = tlx_all_df %>%
    tlx_remove_rand_chromosomes() %>%
    dplyr::mutate(Qname=bait_name, tlx_group=bait_name, tlx_control=F) %>%
    dplyr::filter(tlx_copynumber==1 & !tlx_duplicated & !tlx_is_bait_junction)


  tlxcov_offtargets_df = tlx_offtarget_df %>%
    tlx_coverage(group="group", extsize=offtargets_params$extsize, exttype=offtargets_params$exttype, libfactors_df=libfactors_df, ignore.strand=T, min_sample_pileup=0)
  macs_offtargets = tlxcov_macs2(tlxcov_df=tlxcov_offtargets_df, group="group", params=offtargets_params)


  if(debug) {
    chrom_names = unique(macs_offtargets$islands$tlx_group)
    chrom_colors = apply(col2rgb(randomcoloR::distinctColorPalette(length(chrom_names))), 2, paste, collapse=",")
    names(chrom_colors) = chrom_names

    tlx_offtarget_df %>%
      dplyr::filter(tlx_group=="Chr17_41Mb" & Rname=="chr4" & Junction>=89217093 & Junction<=89217407) %>%
      tlx_write_bed("reports/detect_rdc/offtargets/test", "all", mode="alignment", ignore.strand=T)
    tlx_write_bed(tlx_offtarget_df, "reports/detect_rdc/offtargets/off", "group", mode="alignment", ignore.strand=T)
    tlxcov_write_bedgraph(tlxcov_offtargets_df, "reports/detect_rdc/offtargets/off", "group")
    macs_offtargets$islands %>%
      dplyr::mutate(score=1, strand="*", island_name=paste0(island_name, " (", tlx_group, ")")) %>%
      dplyr::mutate(thickStart=island_start, thickEnd=island_end, score=1, rgb=chrom_colors[tlx_group]) %>%
      dplyr::select(island_chrom, island_start, island_end, island_name, score, strand, thickStart, thickEnd, rgb) %>%
      readr::write_tsv("reports/detect_rdc/offtargets/offtargets.bed", col_names=F)
    macs_offtargets$qvalues %>%
      dplyr::select(qvalue_chrom, qvalue_start, qvalue_end, qvalue_score) %>%
      readr::write_tsv("reports/detect_rdc/offtargets/offtargets2-qvalues.bedgraph", col_names=F)
  }

  offtargets_df = macs_offtargets$islands %>%
    dplyr::group_by(tlx_group) %>%
    dplyr::do((function(df){
      dff<<-df
      tlx_group_ranges = tlx_offtarget_df %>%
        dplyr::filter(tlx_group==df$tlx_group[1]) %>%
        dplyr::select(tlx_group, Rname, Junction) %>%
        df2ranges(Rname, Junction, Junction)
      df %>%
        df2ranges(island_chrom, island_start, island_end) %>%
        GenomicRanges::reduce(min.gapwidth=10e3) %>%
        as.data.frame() %>%
        dplyr::select(offtarget_chrom=seqnames, offtarget_start=start, offtarget_end=end) %>%
        df2ranges(offtarget_chrom, offtarget_start, offtarget_end) %>%
        innerJoinByOverlaps(tlx_group_ranges) %>%
        dplyr::group_by(offtarget_chrom, offtarget_start, offtarget_end) %>%
        dplyr::summarize(offtarget_center=round(mean(Junction))) %>%
        dplyr::mutate(offtarget_start=offtarget_center-50, offtarget_end=offtarget_center+50) %>%
        dplyr::select(offtarget_chrom, offtarget_start, offtarget_end)
    })(.)) %>%
    dplyr::ungroup() %>%
    dplyr::rename(offtarget_bait_name="tlx_group")
  table(offtargets_df$offtarget_bait_name)
  readr::write_tsv(offtargets_df, path="data/offtargets_dkfz.tsv")

  #
  # Merge all off-targets
  #
  offtargets_all_df = macs_offtargets$islands %>%
    df2ranges(island_chrom, island_start, island_end) %>%
    GenomicRanges::reduce(min.gapwidth=1e4) %>%
    as.data.frame() %>%
    dplyr::select(offtarget_chrom=seqnames, offtargets_extended_start=start, offtargets_extended_end=end) %>%
    df2ranges(offtarget_chrom, offtargets_extended_start, offtargets_extended_end) %>%
    innerJoinByOverlaps(macs_offtargets$islands %>% df2ranges(island_chrom, island_start, island_end)) %>%
    dplyr::arrange(offtargets_extended_end-offtargets_extended_start) %>%
    dplyr::distinct(offtarget_chrom, offtargets_extended_start, offtargets_extended_end, .keep_all=T) %>%
    dplyr::mutate(offtarget_start=island_summit_pos-1e5, offtarget_end=island_summit_pos+1e5) %>%
    dplyr::select(dplyr::starts_with("offtarget_")) %>%
    dplyr::mutate(offtarget_name=paste0(offtarget_chrom, ":", offtarget_start))

  #
  # Plot heatmap with all off-targets
  #
  tlx_offtarget_ranges = tlx_offtarget_df %>%
    dplyr::group_by(tlx_sample) %>%
    dplyr::mutate(tlx_sample_size=dplyr::n()) %>%
    dplyr::ungroup() %>%
    dplyr::filter(!tlx_is_bait_junction) %>%
    df2ranges(Rname, Junction, Junction)
  offtargets_map = offtargets_all_df %>%
    df2ranges(offtarget_chrom, offtarget_start, offtarget_end) %>%
    innerJoinByOverlaps(tlx_offtarget_ranges) %>%
    dplyr::group_by(tlx_group, tlx_sample, offtarget_name, tlx_sample_size) %>%
    dplyr::summarize(tlx_count=pmin(sum(tlx_strand=="+"), sum(tlx_strand=="-")), tlx_prop=tlx_count/tlx_sample_size[1], tlx_present=tlx_prop>0.0001)
  offtargets_ann = offtargets_map %>%
    dplyr::group_by(offtarget_name, tlx_group) %>%
    dplyr::summarize(tlx_group_count=sum(tlx_count>0), tlx_group_present=sum(tlx_present), tlx_group_present=factor(ifelse(tlx_group_present>=2, "Yes", "No"))) %>%
    reshape2::dcast(offtarget_name~tlx_group, value.var="tlx_group_present") %>%
    tibble::column_to_rownames("offtarget_name") %>%
    replace(is.na(.), "No")
  offtargets_ann_colors = sapply(colnames(offtargets_ann), simplify=F, FUN=function(z) c("No"="#FFFFFF", "Yes"="#666666"))
  samples_ann = samples_df %>%
    tibble::column_to_rownames("sample") %>%
    dplyr::select(bait_name, experiment)
  offtargets_pheatmap = offtargets_map %>%
    reshape2::dcast(offtarget_name ~ tlx_sample, value.var="tlx_present") %>%
    tibble::column_to_rownames("offtarget_name") %>%
    replace(is.na(.), 0) %>%
    as.matrix()
  pdf("reports/offtargets_map.pdf", width=11.69, height=8.27, paper="a4r")
  pheatmap::pheatmap(
    offtargets_pheatmap,
    annotation_col=samples_ann,
    annotation_row=offtargets_ann,
    annotation_colors =offtargets_ann_colors,
    cutree_cols=6,
    fontsize=5)
  dev.off()
  
  

  #
  # Detect RDC
  #
  tlx_rdc_df = tlx_all_df %>%
    tlx_remove_rand_chromosomes() %>%
    dplyr::filter(tlx_copynumber==1 & !tlx_duplicated) %>%
    dplyr::mutate(tlx_group=ifelse(tlx_is_bait_chrom, "Intra", "Inter"))


  params_rdc = macs2_params(extsize=5e4, exttype="symmetrical", llocal=1e7, minpvalue=0.01, effective_size=1.87e9, maxgap=2e5, minlen=2e5, baseline=2)
  # params_rdc = macs2_params(extsize=5e4, exttype="symmetrical", llocal=1e7, minpvalue=0.01, effective_size=1.87e9, maxgap=0, minlen=1e3, baseline=2)
  tlxcov_rdc_df = tlx_rdc_df %>%
    dplyr::filter(!tlx_is_bait_junction & (!grepl("prom", group) | !tlx_is_bait_chrom)) %>%
    tlx_coverage(group="group", extsize=params_rdc$extsize, exttype=params_rdc$exttype, libfactors_df=libfactors_df, ignore.strand=T)

  macs_rdc = tlxcov_macs2(tlxcov_rdc_df, group="group", params=params_rdc)
  if(debug)
  {
    macs_rdc$qvalues %>%
      dplyr::group_by(tlx_group) %>%
      dplyr::do((function(df){
        df %>%
          dplyr::select(qvalue_chrom, qvalue_start, qvalue_end, qvalue_score) %>%
          readr::write_tsv(paste0("reports/detect_rdc/qvalues-", df$tlx_group[1], ".bedgraph"), col_names = F)
      })(.))
    macs_rdc$islands %>%
      dplyr::mutate(strand="*") %>%
      dplyr::group_by(tlx_group) %>%
      dplyr::do((function(df){
        df %>%
          dplyr::select(island_chrom, island_extended_start, island_extended_end, island_name, island_baseline, strand) %>%
          readr::write_tsv(paste0("reports/detect_rdc/extislands-", df$tlx_group[1], ".bed"), col_names=F)
        df %>%
          dplyr::select(island_chrom, island_start, island_end, island_name, island_baseline, strand) %>%
          readr::write_tsv(paste0("reports/detect_rdc/islands-", df$tlx_group[1], ".bed"), col_names=F)
      })(.))

    tlxcov_rdc_df %>% tlxcov_write_bedgraph(path="reports/detect_rdc/bedgraph", group="group")
    tlx_rdc_df %>% tlx_write_bed(path="reports/detect_rdc/bed", group="group")

    # tlx_rdc_df %>%
    #   tlx_coverage(group="group", extsize=params_rdc$extsize, exttype=params_rdc$exttype, libfactors_df=libfactors_df, ignore.strand=F) %>%
    #   tlxcov_write_bedgraph(path="reports/pdetect_rdc/bedgraph", group="group")
  }


  # tlxcov_df = tlx_df %>%
  #   dplyr::filter(B_Rname==bait_chrom & tlx_is_bait_chrom & !tlx_is_bait_junction) %>%
  #   tlx_coverage(group="group", extsize=params$extsize, exttype=params$exttype, libfactors_df=libfactors_df, ignore.strand=T)

  #
  # Find clusters using MACS3
  #
  macs_results = tlxcov_macs2(tlxcov_df, group="group", params)
  islands_df = macs_results[["islands"]] %>%
    dplyr::inner_join(samples_df %>% dplyr::distinct(group, bait_chrom), by=c("tlx_group"="group"))


  #
  # Remove clusters overlapping with offtarget regions
  #
  islands_ranges = islands_df %>%
    dplyr::mutate(seqnames=island_chrom, start=island_start, end=island_end) %>%
    GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns=T)
  offtarget_df = readr::read_tsv("~/Workspace/Datasets/HTGTS/offtargets_pnas_mm10.tsv")
  offtarget_ranges = offtarget_df %>%
    dplyr::rename(offtarget_bait_chrom="bait_chrom") %>%
    dplyr::mutate(offtarget_start=offtarget_start-5e4, offtarget_end=offtarget_end+5e4) %>%
    GenomicRanges::makeGRangesFromDataFrame(seqnames.field="offtarget_chrom", start.field="offtarget_start", end.field="offtarget_end", ignore.strand=T, keep.extra.columns=T)
  islands_df = as.data.frame(leftJoinByOverlaps(islands_ranges, offtarget_ranges)) %>%
    dplyr::select(-dplyr::matches("_ranges\\.")) %>%
    dplyr::mutate(island_is_offtarget=!is.na(offtarget_bait_chrom) & offtarget_bait_chrom==bait_chrom) %>%
    dplyr::filter(!island_is_offtarget) %>%
    dplyr::arrange(as.numeric(gsub("chr", "", island_chrom)), island_start) %>%
    dplyr::mutate(island_name=stringr::str_glue("MACS_{stringr::str_pad(i, 3, pad='0')}", i=1:dplyr::n()))

  #
  # Write final results
  #
  islands_ranges = islands_df %>%
    dplyr::mutate(seqnames=island_chrom, start=island_start, end=island_end) %>%
    GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns=T)
  island_rdc_df = as.data.frame(leftJoinByOverlaps(islands_ranges, genes_ranges)) %>%
    dplyr::select(-dplyr::matches("_ranges\\.")) %>%
    dplyr::arrange(dplyr::desc(tidyr::replace_na(gene_length, 0))) %>%
    dplyr::distinct(island_name, .keep_all=T)  %>%
    dplyr::arrange(as.numeric(gsub("[^0-9]*(\\d+).*", "\\1", island_name))) %>%
    dplyr::mutate(gene_strand=tidyr::replace_na(gene_strand, "."), gene_id=tidyr::replace_na(gene_id, "-")) %>%
    dplyr::mutate(rdc_group=1, rdc_cluster_display=paste0(island_chrom, ":", island_name, gene_strand)) %>%
    dplyr::select(rdc_chrom=island_chrom, rdc_start=island_start, rdc_end=island_end, rdc_cluster=island_name,	rdc_group, rdc_strand=gene_strand, rdc_gene=gene_id, rdc_cluster_display, rdc_stn=island_snr, rdc_qvalue_log10=island_summit_qvalue, rdc_pileup=island_summit_abs)

  readr::write_tsv(island_rdc_df, "~/Workspace/Datasets/HTGTS/rdc_macs2_mm10.tsv", col_names=T)
  island_rdc_df %>%
    dplyr::select(rdc_chrom, rdc_start, rdc_end, rdc_cluster_display, rdc_qvalue_log10, rdc_strand) %>%
    readr::write_tsv(file="~/Workspace/Datasets/HTGTS/rdc_macs2_mm10.bed", col_names=F)

  macs_results[["qvalues"]] %>%
    dplyr::select(qvalue_chrom, qvalue_start, qvalue_end, qvalue_score) %>%
    readr::write_tsv(file="~/Workspace/Datasets/HTGTS/rdc_macs2_mm10_qvalues.bdg", col_names=F)

  tlxcov_df %>%
    dplyr::select(tlxcov_chrom, tlxcov_start, tlxcov_end, tlxcov_pileup) %>%
    readr::write_tsv(file="~/Workspace/Datasets/HTGTS/rdc_macs2_mm10.bdg", col_names=F)


  islands_df %>% dplyr::filter(grepl("MACS_003", island_name))
}