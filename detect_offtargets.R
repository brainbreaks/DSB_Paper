setwd("~/Workspace/DSB_Paper")

library(readr)
library(dplyr)
library(ggplot2)
library(cowplot)
library(stringr)
library(randomcoloR)
library(Biostrings)
library(ComplexHeatmap)
library(igraph)
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
  # Set-up parameters
  #
  offtargets_params = macs2_params(extsize=150, exttype="opposite", llocal=1e7, minpvalue=0.01, effective_size=1.87e9, seedlen=2, seedgap=10, minlen=2, baseline=2)

  #
  # Filter out data suitable for analysis
  #
  tlx_offtarget_df = tlx_all_df %>%
    tlx_remove_rand_chromosomes() %>%
    dplyr::mutate(tlx_group=bait_name, tlx_control=F) %>%
    dplyr::filter(tlx_copynumber==1 & !tlx_duplicated & !tlx_is_bait_junction)
  tlx_offtarget_ranges = tlx_offtarget_df %>% df2ranges(Rname, Junction, Junction)

  #
  # Use maximal matching to connect allreads to reads with a single reads of opposite direction. Reads without a match are discarded
  #
  tlx_offtarget_plus = tlx_offtarget_df %>%
    dplyr::filter(tlx_strand=="+") %>%
    dplyr::mutate(Tstart=Junction-offtargets_params$extsize, Tend=Rend) %>%
    dplyr::select(tlx_group_plus=tlx_group, Rname_plus=Rname, Tstart_plus=Tstart, Tend_plus=Rend, tlx_sample_plus=tlx_sample, Qname_plus=Qname) %>%
    df2ranges(Rname_plus, Tstart_plus, Tend_plus)
  tlx_offtarget_minus = tlx_offtarget_df %>%
    dplyr::filter(tlx_strand=="-") %>%
    dplyr::mutate(Tstart=Rstart, Tend=Junction+offtargets_params$extsize) %>%
    dplyr::select(tlx_group_minus=tlx_group, Rname_minus=Rname, Tstart_minus=Tstart, Tend_minus=Tend, tlx_sample_minus=tlx_sample, Qname_minus=Qname) %>%
    df2ranges(Rname_minus, Tstart_minus, Tend_minus)
  tlx_offtarget_connected_graph =  innerJoinByOverlaps(tlx_offtarget_minus, tlx_offtarget_plus) %>%
    dplyr::filter(tlx_group_minus==tlx_group_plus) %>%
    dplyr::distinct(Qname_minus, Qname_plus) %>%
    igraph::graph_from_data_frame(directed=T)
  igraph::V(tlx_offtarget_connected_graph)$type = table(tlx_offtarget_df$tlx_strand[match(igraph::V(tlx_offtarget_connected_graph)$name, tlx_offtarget_df$Qname)])
  tlx_offtarget_matching = igraph::max_bipartite_match(tlx_offtarget_connected_graph)$matching
  qnames_paired = unique(unname(c(names(tlx_offtarget_matching), tlx_offtarget_matching))[!is.na(tlx_offtarget_matching)])
  tlx_offtarget_paired_df = tlx_offtarget_df %>%
    dplyr::filter(Qname %in% qnames_paired)


  #
  # Calculate coverage of aligned reads in opposite directions, but remove areas where the pileup consists only from a single direction alignments
  #
  tlxcov_offtargets_paired_df = tlx_offtarget_paired_df %>%
    dplyr::mutate(tlx_group_i=1, tlx_control=F) %>%
    tlx_coverage(group="group", extsize=offtargets_params$extsize, exttype="opposite", libfactors_df=libfactors_df, ignore.strand=F, min_sample_pileup=0)
  tlxcov_offtargets_paired_corrected_df = tlxcov_offtargets_paired_df %>%
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
    })(.)) %>%ucs
    dplyr::mutate(tlx_strand="*") %>%
    dplyr::ungroup()


  #
  # Detect peaks (from now on islands) in prepared pileups
  #
  bgmodel_df = data.frame(bgmodel_distr="pois", bgmodel_lambda=2) %>% tidyr::crossing(tlx_offtarget_df %>% dplyr::distinct(tlx_group, bgmodel_chrom=Rname, bgmodel_strand="*"))
  macs_offtargets = tlxcov_macs2(tlxcov_df=tlxcov_offtargets_paired_corrected_df, bgmodel_df=bgmodel_df, group="group", extended_islands=F, params=offtargets_params)

  #
  # For each island perform fisher test to detect off-ballance in centromeric and telomeric reads.
  # There should be more telomeric alignments downstream from the island center and more centromeric reads upstream.
  # Also perform Kolmogorov–Smirnov test to test for difference in telomeric ad centromeric distribution junctions
  #
  islands_strand_test_df = macs_offtargets$islands %>%
    dplyr::rename(island_tlx_group="tlx_group") %>%
    df2ranges(island_chrom, island_start-100, island_end+100) %>%
    innerJoinByOverlaps(tlx_offtarget_df %>% df2ranges(Rname, Rstart, Rend)) %>%
    dplyr::filter(island_tlx_group==tlx_group)%>%
    dplyr::select(tlx_group, dplyr::starts_with("island_"), tlx_sample, bait_name, Qname, Junction, Rstart, Rend, tlx_strand) %>%
    dplyr::group_by(tlx_group, island_name) %>%
    dplyr::do((function(z){
      z.right = z$Junction>z$island_summit_pos
      z.plus = z$tlx_strand=="+"
      z.table = matrix(c(sum(!z.plus & !z.right), sum(z.plus & !z.right), sum(!z.plus & z.right), sum(z.plus & z.right)), ncol=2, dimnames=list(c("-", "+"), c("<-", "->")))
      z.fisher = fisher.test(z.table)
      if(sum(z$tlx_strand=="+")>2 && sum(z$tlx_strand=="-")>2) {
        z.ks = ks.test(z$Junction[z$tlx_strand=="-"], z$Junction[z$tlx_strand=="+"])$p.value
      } else {
        z.ks = 1
      }

      data.frame(island_strand_pvalue=z.fisher$p.value, island_strand_odds=pmin(1000, z.fisher$estimate), island_strand_ks_pvalue=z.ks, island_strand_n=nrow(z))
    })(.)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(island_strand_pvalue.adj=p.adjust(island_strand_pvalue))
  macs_offtargets$islands = macs_offtargets$islands %>%
    dplyr::select(-dplyr::matches("island_strand_pvalue|island_strand_odds|island_strand_ks_pvalue|island_is_offtarget")) %>%
    dplyr::inner_join(islands_strand_test_df, by=c("tlx_group", "island_name")) %>%
    dplyr::mutate(island_is_offtarget=island_strand_pvalue<=0.01 & island_strand_odds>=2 & island_strand_ks_pvalue<=0.01)
  table(macs_offtargets$islands$tlx_group, macs_offtargets$islands$island_strand_pvalue<=0.01 & macs_offtargets$islands$island_strand_odds>=2 & macs_offtargets$islands$island_strand_ks_pvalue<=0.01)
  table(macs_offtargets$islands$island_strand_pvalue<=0.01 & macs_offtargets$islands$island_strand_odds>=2 & macs_offtargets$islands$island_strand_ks_pvalue<=0.01)

  #
  # Write debugging information
  #
  if(debug) {
    chrom_names = unique(macs_offtargets$islands$tlx_group)
    chrom_colors = apply(col2rgb(randomcoloR::distinctColorPalette(length(chrom_names))), 2, paste, collapse=",")
    names(chrom_colors) = chrom_names

    tlx_write_bed(tlx_offtarget_df, "reports/detect_offtargets/off-raw", group="group", mode="alignment", ignore.strand=T, ignore.treatment=T)
    tlx_write_bed(tlx_offtarget_paired_df, "reports/detect_offtargets/off-paired", group="group", mode="alignment", ignore.strand=T, ignore.treatment=T)
    tlxcov_write_bedgraph(tlxcov_df=tlxcov_offtargets_paired_corrected_df, path="reports/detect_offtargets/off-paired", group="group")

    macs_offtargets$islands %>%
      dplyr::filter(island_is_offtarget) %>%
      dplyr::mutate(score=1, strand="*", island_name=paste0(island_name, " (", tlx_group, ")")) %>%
      dplyr::mutate(thickStart=island_start, thickEnd=island_end, score=1, rgb=chrom_colors[tlx_group]) %>%
      dplyr::select(island_chrom, island_start, island_end, island_name, score, strand, thickStart, thickEnd, rgb) %>%
      readr::write_tsv(paste0("reports/detect_offtargets/offt-islands.bed"), col_names=F)
    macs_offtargets$qvalues %>%
      dplyr::select(qvalue_chrom, qvalue_start, qvalue_end, qvalue_score) %>%
      readr::write_tsv("reports/detect_offtargets/off-qvalues.bedgraph", col_names=F)
  }


  #
  # Find best aligned sequence inside detected offtarget island and prepare final table
  #
  offtargets_exported_df = macs_offtargets$islands %>%
    dplyr::filter(island_is_offtarget) %>%
    dplyr::mutate(offtarget_chrom=island_chrom, offtarget_bait_name=tlx_group, island_region_start=island_start, island_region_end=island_end, offtarget_island_start=island_start, offtarget_island_end=island_end, offtarget_center=island_summit_pos, offtarget_strand_pvalue=island_strand_pvalue, offtarget_strand_odds=island_strand_odds) %>%
    dplyr::inner_join(baits_df, by=c("offtarget_bait_name"="bait_name")) %>%
    dplyr::mutate(sgRNA_sign=ifelse(bait_strand_sgRNA=="+", 1, -1)) %>%
    dplyr::mutate(offtarget_id=paste0(offtarget_chrom, ":", offtarget_center), offtarget_bait_sequence_sgRNA=bait_sequence_sgRNA, coord_begin=ifelse(bait_strand=="+", bait_end-16*sgRNA_sign, bait_start-16*sgRNA_sign), coord_end=ifelse(bait_strand=="+", bait_end+3*sgRNA_sign, bait_start+3*sgRNA_sign), offtarget_bait_start=pmin(coord_begin, coord_end), offtarget_bait_end=pmax(coord_begin, coord_end)) %>%
    dplyr::rename(offtarget_bait_chrom="bait_chrom", offtarget_bait_strand="bait_strand_sgRNA")
  offtargets_exported_df$island_region_sequence = toupper(get_seq("~/Workspace/genomes/mm10/mm10.fa", df2ranges(offtargets_exported_df, offtarget_chrom, island_region_start, island_region_end))$sequence)
  offtargets_exported_sequences_df = get_pairwise_alignment(paste0(offtargets_exported_df$offtarget_bait_sequence_sgRNA, "NGG"), offtargets_exported_df$island_region_sequence, gapOpening=2, gapExtension=0.5) %>%
    dplyr::select(offtarget_alignment_score=score, offtarget_alignment_pid=pid, offtarget_alignment_start=start, offtarget_alignment_end=end) %>%
    dplyr::bind_cols(offtargets_exported_df) %>%
    dplyr::mutate(offtarget_sequence_start=island_region_start+offtarget_alignment_start-1, offtarget_sequence_end=island_region_start+offtarget_alignment_end-3)

  # offtargets_exported_sequences_df = dplyr::bind_rows(offtargets_exported_df %>% dplyr::mutate(subset="signal"), offtargets_exported_df %>% tidyr::crossing(i=1:100, subset="random")) %>%
  #   dplyr::summarize(cbind(subset=subset, island_region_sequence=island_region_sequence, get_pairwise_alignment(paste0(offtarget_bait_sequence_sgRNA, "NGG"), stringi::stri_rand_shuffle(island_region_sequence), gapOpening=2, gapExtension=0.5))) %>%
  #   dplyr::select(subset, island_region_sequence, offtarget_alignment_score=score, offtarget_alignment_pid=pid, offtarget_alignment_start=start, offtarget_alignment_end=end) %>%
  #   dplyr::inner_join(offtargets_exported_df, by="island_region_sequence") %>%
  #   dplyr::mutate(offtarget_sequence_start=island_region_start+offtarget_alignment_start-1, offtarget_sequence_end=island_region_start+offtarget_alignment_end-3)
  # ggplot(offtargets_exported_sequences_df) +
  #   geom_density(aes(x=offtarget_alignment_score, fill=subset), alpha=0.5, data=offtargets_exported_sequences_df)
  # ggplot(offtargets_exported_sequences_df) +
  #   geom_point(aes(x=offtarget_alignment_score, y=offtarget_id, fill=subset), data=offtargets_exported_sequences_df %>% dplyr::filter(subset=="signal")) +
  #   ggridges::geom_density_ridges(aes(x=offtarget_alignment_score, y=offtarget_id, fill=subset), alpha=0.2, data=offtargets_exported_sequences_df %>% dplyr::filter(subset=="random"))
  #   facet_grid(offtarget_bait_name~., scales="free")

  offtargets_exported_sequences_df %>%
    dplyr::mutate(score=1, strand="*", thickStart=offtarget_island_start, thickEnd=offtarget_island_end, score=1, rgb=chrom_colors[offtarget_bait_name]) %>%
    dplyr::select(offtarget_chrom, offtarget_island_start, offtarget_island_end, offtarget_bait_name, score, strand, thickStart, thickEnd, rgb) %>%
    readr::write_tsv("data/offtargets_dkfz.bed", col_names=F)

  offtargets_exported_sequences_df %>%
    dplyr::select(
      offtarget_bait_name, offtarget_bait_chrom, offtarget_bait_start, offtarget_bait_end, offtarget_bait_strand, offtarget_chrom,
      offtarget_island_start, offtarget_island_end, offtarget_strand_pvalue, offtarget_strand_odds,
      offtarget_sequence_start, offtarget_sequence_end, offtarget_alignment_score, offtarget_alignment_pid) %>%
      readr::write_tsv(file="data/offtargets_dkfz.tsv")

  #
  # Plot heatmap with all off-targets
  #
  offtargets_exported_sequences_ranges = offtargets_exported_sequences_df %>% df2ranges(offtarget_chrom, offtarget_island_start, offtarget_island_end)
  offtargets_all_df = offtargets_exported_sequences_ranges %>%
    GenomicRanges::reduce(min.gapwidth=1e4) %>%
    as.data.frame() %>%
    dplyr::mutate(offtarget_extended_chrom=seqnames, offtarget_extended_start=start-1e4, offtarget_extended_end=end+1e4, offtarget_name=paste0(offtarget_extended_chrom, ":", offtarget_extended_start)) %>%
    dplyr::select(dplyr::starts_with("offtarget_"))
  tlx_offtarget_ranges = tlx_offtarget_df %>%
    dplyr::group_by(tlx_sample, Rname) %>%
    dplyr::mutate(tlx_sample_size=dplyr::n()) %>%
    dplyr::ungroup() %>%
    dplyr::filter(!tlx_is_bait_junction) %>%
    df2ranges(Rname, Junction, Junction)
  offtargets_map = offtargets_all_df %>%
    df2ranges(offtarget_extended_chrom, offtarget_extended_start, offtarget_extended_end) %>%
    innerJoinByOverlaps(tlx_offtarget_ranges) %>%
    dplyr::group_by(tlx_group, tlx_sample, bait_name, offtarget_name, offtarget_extended_chrom, offtarget_extended_start, offtarget_extended_end, tlx_sample_size) %>%
    dplyr::summarize(tlx_count=pmin(sum(tlx_strand=="+"), sum(tlx_strand=="-")), tlx_prop=tlx_count/tlx_sample_size[1], tlx_present=tlx_prop>0.0001) %>%
    dplyr::group_by(tlx_group, offtarget_name) %>%
    dplyr::filter(sum(tlx_present) > 3) %>%
    dplyr::ungroup()
  offtargets_pheatmap = offtargets_map %>%
    reshape2::dcast(offtarget_name ~ tlx_sample, value.var="tlx_prop") %>%
    tibble::column_to_rownames("offtarget_name") %>%
    replace(is.na(.), 0) %>%
    as.matrix()

  offtarget_name_order = offtargets_map %>%
    tidyr::separate(offtarget_name, into=c("offtarget_name_chrom", "offtarget_name_pos"), sep=":", remove=F) %>%
    dplyr::group_by(offtarget_name, offtarget_name_chrom, offtarget_name_pos, bait_name) %>%
    dplyr::summarize(bait_mean=max(tlx_prop)) %>%
    dplyr::arrange(dplyr::desc(bait_mean)) %>%
    dplyr::distinct(offtarget_name, .keep_all=T) %>%
    dplyr::arrange(bait_name, offtarget_name_chrom, offtarget_name_pos) %>%
    .$offtarget_name
  sample_order = unique(offtargets_map %>% dplyr::arrange(bait_name) %>% .$tlx_sample)
  offtargets_pheatmap = offtargets_pheatmap[offtarget_name_order, sample_order]

  samples_ann = tlx_offtarget_df %>%
    dplyr::filter(!tlx_duplicated & !tlx_is_bait_junction & !tlx_is_bait_chrom & tlx_copynumber==1) %>%
    dplyr::group_by(sample=tlx_sample, experiment, bait_name) %>%
    dplyr::summarize(size=dplyr::n()) %>%
    dplyr::filter(sample %in% sample_order) %>%
    dplyr::mutate(order=match(sample, sample_order), bait_name=gsub("_.*", "", bait_name)) %>%
    dplyr::filter(!is.na(order)) %>%
    dplyr::arrange(order) %>%
    tibble::column_to_rownames("sample") %>%
    dplyr::select(bait_name, experiment, bait_name, size)

  # bait=ComplexHeatmap::anno_text(samples_ann$bait_name, gp=gpar(fontsize = 6)),
  pdf("reports/offtargets_map_new.pdf", width=2*11.69, height=2*8.27)
  ComplexHeatmap::Heatmap(offtargets_pheatmap, row_names_gp=gpar(fontsize=6), column_names_gp=gpar(fontsize = 6),
    cluster_columns=F, cluster_rows=F, column_split=samples_ann$bait_name,
    top_annotation = ComplexHeatmap::HeatmapAnnotation(experiment=samples_ann$experiment, sample_size=samples_ann$size)
  )
  dev.off()
}