setwd("~/Workspace/DSB_Paper")

library(readr)
library(dplyr)
library(ggplot2)
library(cowplot)
library(stringr)
library(randomcoloR)
library(ComplexHeatmap)
devtools::load_all("~/Workspace/breaktools/")


detect_rdc = function()
{
  debug=T
  dir.create("reports/detect_offtargets", recursive=T)

  #
  # Load samples
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
  
  libfactors_df = tlx_all_df %>% tlx_libsizes()
  
  #
  # Detect offtargets
  #
  offtargets_params = macs2_params(extsize=50, exttype="opposite", llocal=1e7, minqvalue=0.01, baseline=2, effective_size=1.87e9, maxgap=10e3, minlen=2)
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

    tlx_write_bed(tlx_offtarget_df, "reports/detect_offtargets/off", "group", mode="alignment", ignore.strand=T)
    tlxcov_write_bedgraph(tlxcov_offtargets_df, "reports/detect_offtargets/off", "group")
    macs_offtargets$islands %>%
      dplyr::mutate(score=1, strand="*", island_name=paste0(island_name, " (", tlx_group, ")")) %>%
      dplyr::mutate(thickStart=island_start, thickEnd=island_end, score=1, rgb=chrom_colors[tlx_group]) %>%
      dplyr::select(island_chrom, island_start, island_end, island_name, score, strand, thickStart, thickEnd, rgb) %>%
      readr::write_tsv("reports/detect_offtargets/offtargets.bed", col_names=F)
    # macs_offtargets$qvalues %>%
    #   dplyr::select(qvalue_chrom, qvalue_start, qvalue_end, qvalue_score) %>%
    #   readr::write_tsv("reports/detect_offtargets/off-qvalues.bedgraph", col_names=F)
  }

  #
  # Write results
  #
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
  readr::write_tsv(offtargets_df, file="data/offtargets_dkfz.tsv")

  #
  # Plot heatmap with all off-targets
  #
  offtargets_ranges = offtargets_df %>% df2ranges(offtarget_chrom, offtarget_start, offtarget_end)
  offtargets_all_df = offtargets_ranges %>%
    GenomicRanges::reduce(min.gapwidth=1e4) %>%
    as.data.frame() %>%
    dplyr::mutate(offtarget_extended_chrom=seqnames, offtarget_extended_start=start-1e5, offtarget_extended_end=end+1e5, offtarget_name=paste0(offtarget_extended_chrom, ":", offtarget_extended_start)) %>%
    dplyr::select(dplyr::starts_with("offtarget_"))
  tlx_offtarget_ranges = tlx_offtarget_df %>%
    dplyr::group_by(tlx_sample) %>%
    dplyr::mutate(tlx_sample_size=dplyr::n()) %>%
    dplyr::ungroup() %>%
    dplyr::filter(!tlx_is_bait_junction) %>%
    df2ranges(Rname, Junction, Junction)
  offtargets_map = offtargets_all_df %>%
    df2ranges(offtarget_extended_chrom, offtarget_extended_start, offtarget_extended_end) %>%
    innerJoinByOverlaps(tlx_offtarget_ranges) %>%
    dplyr::group_by(tlx_group, tlx_sample, offtarget_name, offtarget_extended_chrom, offtarget_extended_start, offtarget_extended_end, tlx_sample_size) %>%
    dplyr::summarize(tlx_count=pmin(sum(tlx_strand=="+"), sum(tlx_strand=="-")), tlx_prop=tlx_count/tlx_sample_size[1], tlx_present=tlx_prop>0.0001) %>%
    dplyr::group_by(tlx_group, offtarget_name) %>%
    dplyr::filter(sum(tlx_present) > 3) %>%
    dplyr::ungroup()
  offtargets_pheatmap = offtargets_map %>%
    reshape2::dcast(offtarget_name ~ tlx_sample, value.var="tlx_present") %>%
    tibble::column_to_rownames("offtarget_name") %>%
    replace(is.na(.), 0) %>%
    as.matrix()
  offtargets_ann = offtargets_map %>%
    df2ranges(offtarget_extended_chrom, offtarget_extended_start, offtarget_extended_end) %>%
    innerJoinByOverlaps(offtargets_ranges) %>%
    dplyr::mutate(has_bait=T) %>%
    reshape2::dcast(offtarget_name~offtarget_bait_name, value.var="has_bait", fun.aggregate=any) %>%
    replace(is.na(.), 0) %>%
    tibble::column_to_rownames("offtarget_name")
  samples_ann = samples_df %>%
    dplyr::mutate(order=match(sample, colnames(offtargets_pheatmap)), bait_name=gsub("_.*", "", bait_name)) %>%
    dplyr::filter(!is.na(order)) %>%
    dplyr::arrange(order) %>%
    tibble::column_to_rownames("sample") %>%
    dplyr::select(bait_name, experiment)
  pdf("reports/offtargets_map.pdf", width=11.69, height=8.27, paper="a4r")
  offtargets_anncol = sapply(colnames(offtargets_ann), function(z) c("FALSE"="#FFFFFF", "TRUE"="#666666"), simplify=F)
  ComplexHeatmap::Heatmap(offtargets_pheatmap, row_names_gp = gpar(fontsize = 6), column_names_gp = gpar(fontsize = 6),
    top_annotation = ComplexHeatmap::HeatmapAnnotation(bait=ComplexHeatmap::anno_text(samples_ann$bait_name, gp=gpar(fontsize = 6)), experiment=samples_ann$experiment),
    right_annotation = do.call(ComplexHeatmap::rowAnnotation, c(as.list(offtargets_ann), list(col=offtargets_anncol, show_legend=F)))
  )
  dev.off()
}