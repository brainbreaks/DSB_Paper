setwd("~/Workspace/DSB_Paper")

library(readr)
library(dplyr)
library(ggplot2)
library(cowplot)
library(stringr)
library(randomcoloR)
library(Biostrings)
library(ComplexHeatmap)
devtools::load_all("~/Workspace/breaktools/")


detect_offtargets = function()
{
  debug=T
  dir.create("reports/detect_offtargets", recursive=T)

  #
  # Load baits
  #
  baits_df = readr::read_tsv("~/Workspace/Datasets/HTGTS/dkfz_baits.tsv")

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
  # save(tlx_all_df, samples_df, baits_df, libfactors_df, file="detect_offtargets.rda")
  # load("detect_offtargets.rda")

  
  #
  # Detect offtargets
  #
  offtargets_params = macs2_params(extsize=50, exttype="opposite", llocal=1e7, minqvalue=0.01, effective_size=1.87e9, seedlen=2, seedgap=1e2, minlen=2, baseline=2)
  tlx_offtarget_df = tlx_all_df %>%
    tlx_remove_rand_chromosomes() %>%
    dplyr::mutate(Qname=bait_name, tlx_group=bait_name, tlx_control=F) %>%
    dplyr::filter(tlx_copynumber==1 & !tlx_duplicated & !tlx_is_bait_junction)

  tlxcov_offtargets_strand_df = tlx_offtarget_df %>%
    tlx_coverage(group="group", extsize=offtargets_params$extsize, exttype=offtargets_params$exttype, libfactors_df=libfactors_df, ignore.strand=F, min_sample_pileup=0)

  tlxcov_offtargets_corrected_df = tlxcov_offtargets_strand_df %>%
    dplyr::group_by(tlx_group) %>%
    dplyr::do((function(tlx_g, tlx_control) {
      tlx_gg <<- tlx_g
      tlx_g_ranges = tlx_g %>% df2ranges(tlxcov_chrom, tlxcov_start, tlxcov_end)
      res_df = tlx_g %>%
        reshape2::melt(measure.vars=c("tlxcov_start", "tlxcov_end"), value.name="tlxcov_pos") %>%
        dplyr::distinct(tlxcov_sum_chrom=tlxcov_chrom, tlx_sum_control=tlx_control, tlxcov_pos) %>%
        dplyr::arrange(tlxcov_sum_chrom, tlxcov_pos) %>%
        dplyr::mutate(tlxcov_sum_start=dplyr::lag(tlxcov_pos), tlxcov_sum_end=tlxcov_pos-1, tlxcov_sum_start=ifelse(tlxcov_sum_start>tlxcov_sum_end, 1, tlxcov_sum_start)) %>%
        dplyr::filter(!is.na(tlxcov_sum_start)) %>%
        dplyr::select(tlxcov_sum_chrom, tlxcov_sum_start, tlxcov_sum_end) %>%
        df2ranges(tlxcov_sum_chrom, tlxcov_sum_start, tlxcov_sum_end) %>%
        innerJoinByOverlaps(tlx_g_ranges) %>%
        reshape2::dcast(tlx_group+tlx_control+tlxcov_sum_chrom+tlxcov_sum_start+tlxcov_sum_end ~ tlx_strand, value.var="tlxcov_pileup") %>%
        dplyr::mutate(`+`=tidyr::replace_na(`+`,0), `-`=tidyr::replace_na(`-`,0)) %>%
        dplyr::mutate(tlxcov_pileup=pmin(`+`, `-`)*2) %>%
        dplyr::select(tlx_group, tlx_control, tlxcov_chrom=tlxcov_sum_chrom, tlxcov_start=tlxcov_sum_start, tlxcov_end=tlxcov_sum_end, tlxcov_pileup)
      res_df
    })(.)) %>%
    dplyr::mutate(tlx_strand="*") %>%
    dplyr::ungroup()

  bgmodel_df = data.frame(bgmodel_distr="pois", bgmodel_lambda=2) %>% tidyr::crossing(tlxcov_offtargets_df %>% dplyr::distinct(tlx_group, bgmodel_chrom=tlxcov_chrom, bgmodel_strand=tlx_strand))
  macs_offtargets = tlxcov_macs2(tlxcov_df=tlxcov_offtargets_corrected_df, bgmodel_df=bgmodel_df, group="group", extended_islands=F, params=offtargets_params)

  tlx_offtarget_ranges = tlx_offtarget_df %>% df2ranges(Rname, Rstart, Rend)
  islands_strand_test_df = macs_offtargets$islands %>%
    dplyr::rename(island_tlx_group="tlx_group") %>%
    dplyr::mutate(island_center_start=island_summit_pos-100, island_center_end=island_summit_pos+100) %>%
    df2ranges(island_chrom, island_center_start, island_center_end) %>%
    innerJoinByOverlaps(tlx_offtarget_ranges) %>%
    dplyr::filter(island_tlx_group==tlx_group) %>%
    dplyr::select(tlx_group, dplyr::starts_with("island_"), Junction, tlx_strand) %>%
    dplyr::group_by(tlx_group, island_name) %>%
    dplyr::do((function(z){
      zz<<-z
      # z = islands_strand_test_df %>% dplyr::filter(island_name=="RDC_108" & tlx_group=="Chr6_70Mb")
      # z = islands_strand_test_df %>% dplyr::filter(island_name=="RDC_064" & tlx_group=="Chr15_86Mb")
      # z = islands_strand_test_df %>% dplyr::filter(island_name=="RDC_116" & tlx_group=="Chr4_51Mb")
      # z = islands_strand_test_df %>% dplyr::filter(island_name=="RDC_001" & tlx_group=="Chr1_41Mb")

      z.right = z$Junction>z$island_summit_pos
      z.plus = z$tlx_strand=="+"
      z.table = matrix(c(sum(!z.plus & !z.right), sum(z.plus & !z.right), sum(!z.plus & z.right), sum(z.plus & z.right)), ncol=2)
      z.fisher = fisher.test(z.table)
      data.frame(island_strand_pvalue=z.fisher$p.value, island_strand_odds=z.fisher$estimate)
    })(.)) %>%
    dplyr::ungroup()
  macs_offtargets$islands = macs_offtargets$islands %>%
    dplyr::select(-dplyr::matches("island_strand_pvalue|island_strand_odds")) %>%
    dplyr::inner_join(islands_strand_test_df, by=c("tlx_group", "island_name"))
  table(macs_offtargets$islands$tlx_group, macs_offtargets$islands$island_strand_pvalue<=0.01 & macs_offtargets$islands$island_strand_odds>=2)


  if(debug) {
    chrom_names = unique(macs_offtargets$islands$tlx_group)
    chrom_colors = apply(col2rgb(randomcoloR::distinctColorPalette(length(chrom_names))), 2, paste, collapse=",")
    names(chrom_colors) = chrom_names

    tlx_write_bed(tlx_offtarget_df, "reports/detect_offtargets/off", "group", mode="alignment", ignore.strand=T)
    tlxcov_write_bedgraph(tlxcov_df=tlxcov_offtargets_corrected_df, path="reports/detect_offtargets/off", group="group")
    macs_offtargets$islands %>%
      dplyr::filter(island_strand_pvalue<=0.01 & island_strand_odds>=2) %>%
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
    dplyr::filter(island_strand_pvalue<=0.01 & island_strand_odds>=2) %>%
    dplyr::group_by(tlx_group) %>%
    dplyr::do((function(df){
      dff<<-df
      tlx_group_ranges = tlx_offtarget_df %>%
        dplyr::filter(tlx_group==df$tlx_group[1]) %>%
        dplyr::select(tlx_group, Rname, Junction) %>%
        df2ranges(Rname, Junction, Junction)
      df %>%
        # df2ranges(island_chrom, island_start, island_end) %>%
        # GenomicRanges::reduce(min.gapwidth=10e3) %>%
        # as.data.frame() %>%
        # dplyr::select(offtarget_chrom=seqnames, offtarget_start=start, offtarget_end=end) %>%
        dplyr::select(offtarget_chrom=island_chrom, offtarget_start=island_start, offtarget_end=island_end) %>%
        df2ranges(offtarget_chrom, offtarget_start, offtarget_end) %>%
        innerJoinByOverlaps(tlx_group_ranges) %>%
        dplyr::group_by(offtarget_chrom, offtarget_start, offtarget_end) %>%
        dplyr::summarize(offtarget_center=round(mean(Junction))) %>%
        dplyr::mutate(offtarget_region_start=offtarget_center-50, offtarget_region_end=offtarget_center+50) %>%
        dplyr::select(offtarget_chrom, offtarget_region_start, offtarget_region_end)
    })(.)) %>%
    dplyr::ungroup() %>%
    dplyr::rename(offtarget_bait_name="tlx_group")

  #
  # Export offtargets
  #
  offtargets_exported_df = offtargets_df %>%
    dplyr::inner_join(baits_df, by=c("offtarget_bait_name"="bait_name")) %>%
    dplyr::mutate(sgRNA_sign=ifelse(bait_strand_sgRNA=="+", 1, -1)) %>%
    dplyr::mutate(offtarget_bait_sequence_sgRNA=bait_sequence_sgRNA, coord_begin=ifelse(bait_strand=="+", bait_end-16*sgRNA_sign, bait_start-16*sgRNA_sign), coord_end=ifelse(bait_strand=="+", bait_end+3*sgRNA_sign, bait_start+3*sgRNA_sign), offtarget_bait_start=pmin(coord_begin, coord_end), offtarget_bait_end=pmax(coord_begin, coord_end)) %>%
    dplyr::rename(offtarget_bait_chrom="bait_chrom", offtarget_bait_strand="bait_strand_sgRNA") %>%
    dplyr::do(as.data.frame(get_seq("~/Workspace/genomes/mm10/mm10.fa", df2ranges(., offtarget_chrom, offtarget_region_start, offtarget_region_end))) %>% dplyr::rename(offtarget_region_sequence="sequence")) %>%
    dplyr::mutate(offtarget_region_sequence=toupper(offtarget_region_sequence)) %>%
    dplyr::select(dplyr::starts_with("offtarget_")) %>%
    dplyr::bind_cols(get_pairwise_alignment(paste0(.$offtarget_bait_sequence_sgRNA, "NGG"), .$offtarget_region_sequence, gapOpening=2, gapExtension=0.5)) %>%
    dplyr::mutate(offtarget_start=offtarget_region_start+start-1, offtarget_end=offtarget_region_start+end-3) %>%
    dplyr::select(offtarget_bait_name, offtarget_bait_chrom, offtarget_bait_start, offtarget_bait_end, offtarget_bait_strand, offtarget_chrom, offtarget_start, offtarget_end)

  table(offtargets_exported_df$offtarget_bait_name)
  readr::write_tsv(offtargets_exported_df, file="data/offtargets_dkfz.tsv")
  offtargets_exported_df %>%
      dplyr::mutate(score=1, strand="*", thickStart=offtarget_start, thickEnd=offtarget_end, score=1, rgb=chrom_colors[offtarget_bait_name]) %>%
      dplyr::select(offtarget_chrom, offtarget_start, offtarget_end, offtarget_bait_name, score, strand, thickStart, thickEnd, rgb) %>%
      readr::write_tsv("data/offtargets_dkfz.bed", col_names=F)
  readr::write_tsv(offtargets_exported_df, file="data/offtargets_dkfz.tsv")

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