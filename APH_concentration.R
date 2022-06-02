setwd("~/Workspace/DSB_Paper")

library(readr)
library(dplyr)
library(ggplot2)
library(cowplot)
library(stringr)
devtools::load_all('~/Workspace/breaktools/')

APH_concentration = function()
{
  group_palette = c("APH 0.2 uM 96h"="#C6DBEF", "APH 0.3 uM 96h"="#6DAACE", "APH 0.4 uM 96h"="#317BA5", "APH 0.6 uM 96h"="#335E9D", "DMSO"="#CCCCCC")
  group_alpha = c("APH 0.2 uM 96h"=0.7, "APH 0.3 uM 96h"=0.8, "APH 0.4 uM 96h"=0.9, "APH 0.6 uM 96h"=1, "DMSO"=0.6)

  baits_df = readr::read_tsv("~/Workspace/Datasets/HTGTS/wei_pnas2018_baits.tsv")

  #
  # Read TLX
  #
  samples_df = tlx_read_samples("~/Workspace/Datasets/HTGTS/samples/All_samples.tsv", "~/Workspace/Datasets/HTGTS") %>%
    dplyr::filter(grepl("concentration", experiment))

  tlx_all_df = tlx_read_many(samples_df, threads=30)
  tlx_df = tlx_all_df %>%
    tlx_extract_bait(bait_size=19, bait_region=12e6) %>%
    tlx_remove_rand_chromosomes() %>%
    tlx_mark_dust() %>%
    tlx_calc_copynumber(bowtie2_index="~/Workspace/genomes/mm10/mm10", max_hits=100, threads=24) %>%
    dplyr::filter(tlx_copynumber==1 & !tlx_duplicated) %>%
    dplyr::group_by(tlx_sample) %>%
    dplyr::filter(dplyr::n()>5000) %>%
    dplyr::ungroup()




  libfactors_bait_df = tlx_libfactors(tlx_df %>% dplyr::mutate(tlx_group=bait_chrom), normalize_within="none", normalize_between="none", normalization_target="min")

  #
  # Indentify offtarget candidates
  #
  offtargets_params = macs2_params(extsize=50, exttype="opposite", llocal=1e7, minqvalue=0.01, baseline=2, effective_size=1.87e9, maxgap=1e3, minlen=10)
  tlx_offtarget_df = tlx_df %>%
    dplyr::filter(!tlx_control & tlx_is_bait_chrom & !tlx_is_bait_junction) %>%
    dplyr::mutate(tlx_group=bait_chrom)
  tlxcov_offtargets_df = tlx_offtarget_df %>%
    tlx_coverage(group="group", extsize=offtargets_params$extsize, exttype=offtargets_params$exttype, libfactors_df=libfactors_bait_df, ignore.strand=T)
  macs_offtargets = tlxcov_macs2(tlxcov_offtargets_df, group="group", offtargets_params)


  # Export debuging info
  tlx_write_bed(tlx_offtarget_df, "reports/APH_concentration/offtargets", "all", mode="alignment")
  tlxcov_write_bedgraph(tlxcov_offtargets_df, "reports/APH_concentration/offtargets", "all")
  macs_offtargets$islands %>%
    dplyr::mutate(score=1, strand="*") %>%
    dplyr::select(island_chrom, island_start, island_end, island_name, score, strand) %>%
    readr::write_tsv("reports/APH_concentration/offtargets-islands.bed", col_names=F)
  macs_offtargets$qvalues %>%
    dplyr::select(qvalue_chrom, qvalue_start, qvalue_end, qvalue_score) %>%
    readr::write_tsv("reports/APH_concentration/offtargets-qvalues.bedgraph", col_names = F)


  #
  # Mark offtargets
  #
  offtargets_df = readr::read_tsv("~/Workspace/Datasets/HTGTS/offtargets_pnas_mm10.tsv")
  offtargets_extended_df = dplyr::bind_rows(
    offtargets_df,
    macs_offtargets$islands %>%
      dplyr::inner_join(offtargets_df %>% dplyr::distinct(offtarget_bait_chrom, offtarget_bait_start, offtarget_bait_end, offtarget_bait_strand), by=c("island_chrom"="offtarget_bait_chrom")) %>%
      dplyr::mutate(offtarget_is_primary=F) %>%
      dplyr::select(
        offtarget_bait_chrom=island_chrom, offtarget_bait_start, offtarget_bait_end, offtarget_bait_strand,
        offtarget_chrom=island_chrom, offtarget_start=island_start, offtarget_end=island_end, offtarget_is_primary)
  )

  tlx_clean_df = tlx_df %>%
    tlx_mark_offtargets(offtargets_extended_df, offtarget_region=1e6, bait_region=1e4) %>%
    dplyr::filter(!tlx_is_offtarget & tlx_is_bait_chrom)

  #
  # Use all concentrations to call RDC
  #
  params_clean = macs2_params(extsize=5e4, exttype="symmetrical", llocal=1e7, minqvalue=0.01, effective_size=1.87e9, maxgap=1e6, minlen=1e5)
  tlxcov_clean_df = tlx_clean_df %>%
    dplyr::filter(!tlx_control) %>%
    dplyr::mutate(tlx_group=bait_chrom) %>%
    tlx_coverage(group="group", extsize=params_clean$extsize, exttype=params_clean$exttype, libfactors_df=libfactors_bait_df, ignore.strand=T)
  macs_clean = tlxcov_macs2(tlxcov_clean_df, group="group", params_clean)

  # Export debuging info
  macs_clean$islands %>%
    dplyr::mutate(score=1, strand="*") %>%
    dplyr::select(island_chrom, island_extended_start, island_extended_end, island_name, score, strand) %>%
    readr::write_tsv("reports/APH_concentration/all-islands.bed", col_names=F)
  macs_clean$qvalues %>%
    dplyr::select(qvalue_chrom, qvalue_start, qvalue_end, qvalue_score) %>%
    readr::write_tsv("reports/APH_concentration/all-qvalues.bedgraph", col_names = F)
  tlxcov_write_bedgraph(tlxcov_clean_df, "reports/APH_concentration/all", "all")
  tlx_clean_df %>%
    dplyr::mutate(score=1) %>%
    dplyr::select(Rname, Rstart, Rend, Qname, score, tlx_strand) %>%
    readr::write_tsv("reports/APH_concentration/all-reference.bed", col_names=F)


  #
  # Each sample separately
  #
  params_baitconcentration = macs2_params(extsize=100e3, exttype="symmetrical")
  tlx_baitconcentration_df = tlx_clean_df %>% dplyr::mutate(tlx_group=paste0(tlx_group, " (", bait_chrom, ")"))
  tlx_baitconcentration_ranges = tlx_clean_df %>% df2ranges(Rname, Junction, Junction)
  libfactors_baitconcentration_df = tlx_libfactors(tlx_baitconcentration_df, normalize_within="none", normalize_between="none", normalization_target="min")
  tlxcov_baitconcentration_strand_df = tlx_clean_df %>%
    tlx_coverage(group="group", extsize=params_baitconcentration$extsize, exttype=params_baitconcentration$exttype, libfactors_df=libfactors_baitconcentration_df, ignore.strand=F)
  tlxcov_baitconcentration_strand_ranges = tlxcov_baitconcentration_strand_df %>% df2ranges(tlxcov_chrom, tlxcov_start, tlxcov_end)


  #
  # Load RDC
  #
  # rdc_df = readr::read_tsv("~/Workspace/Datasets/HTGTS/rdc_macs2_mm10.tsv")
  rdc_df = macs_clean$islands %>%
    dplyr::select(rdc_chrom=island_chrom, rdc_start=island_extended_start, rdc_end=island_extended_end, rdc_cluster=island_name, island_summit_abs)

  #
  # Write strand-stratified bedgraph
  #
  tlxcov_baitconcentration_strand_df %>%
    dplyr::mutate(tlx_group=gsub(" \\(.*", "", tlx_group)) %>%
    tlxcov_write_bedgraph(path="reports/APH_concentration/concentration", group="group")

  rdc_junctions_df = rdc_df %>%
    df2ranges(rdc_chrom, rdc_start, rdc_end) %>%
    innerJoinByOverlaps(tlx_baitconcentration_ranges) %>%
    dplyr::inner_join(baits_df, by=c("rdc_chrom"="bait_chrom")) %>%
    dplyr::group_by(rdc_cluster, rdc_chrom, rdc_start, rdc_end, bait_start, rdc_tlx_group=tlx_group, island_summit_abs) %>%
    dplyr::summarize(n_sense=sum(tlx_strand=="+"), n_anti=sum(tlx_strand=="-")) %>%
    dplyr::ungroup()
    # dplyr::mutate(n_sense=n_sense+n_sense*(rdc_start<bait_start)*1.1, n_anti=n_anti+n_anti*(rdc_start>bait_start)*1.1)
  ccs_df = rdc_junctions_df %>%
    df2ranges(rdc_chrom, rdc_start, rdc_end) %>%
    innerJoinByOverlaps(tlxcov_baitconcentration_strand_ranges) %>%
    dplyr::filter(rdc_tlx_group==tlx_group) %>%
    dplyr::group_by(rdc_cluster) %>%
    dplyr::mutate(pileup_sense=max(tlxcov_pileup[tlx_strand=="+"], na.rm=T), pileup_anti=max(tlxcov_pileup[tlx_strand=="-"], na.rm=T)) %>%
    dplyr::group_by(rdc_cluster, rdc_chrom, tlx_group, n_sense, n_anti, island_summit_abs, pileup_sense, pileup_anti) %>%
    summarize(tlx_strand_crosscorrelation(dplyr::cur_data(), step=5000)) %>%
    dplyr::mutate(tlx_group_int=as.numeric(factor(tlx_group))) %>%
    dplyr::filter(!is.na(crosscorrelation_lag) & crosscorrelation_value>0.5 & pmin(pileup_sense, pileup_anti) >= 10)

  tlx_shift_df = rdc_junctions_df %>%
    df2ranges(rdc_chrom, rdc_start, rdc_end) %>%
    innerJoinByOverlaps(tlx_baitconcentration_ranges) %>%
    dplyr::filter(rdc_tlx_group==tlx_group) %>%
    dplyr::group_by(rdc_cluster, rdc_chrom, tlx_group, n_sense, n_anti) %>%
    dplyr::summarize(tlx_strand_shift=mean(Junction[tlx_strand=="+"])-mean(Junction[tlx_strand=="-"]), tlx_strand_relshift=tlx_strand_shift/(max(Junction)-min(Junction))) %>%
    dplyr::filter(pmin(n_sense, n_anti)>10)

  pdf("reports/APH_concentration5.pdf", width=11.69, height=8.27, paper="a4r")
  ggplot(rdc_junctions_df, aes(y=log2(n_sense/n_anti), x=rdc_tlx_group)) +
    geom_hline(yintercept=0, linetype="dashed") +
    geom_boxplot(aes(fill=rdc_tlx_group), outlier.shape = NA, outlier.alpha = 0) +
    geom_point(aes(x=rdc_tlx_group, size=n_sense+n_anti), position=position_jitter(width=0.2), alpha=0.3) +
    # geom_text(aes(x=rdc_tlx_group, label=gsub("MACS3_", "", rdc_cluster), size=n_sense+n_anti), color="#FF0000", position=position_jitter(width=0.2)) +
    labs(y="Proportion of right-moving fork breaks, log2", size="No. of breaks", fill="APH concentration") +
    scale_fill_manual(values=group_palette) +
    theme_bw(base_size=12) +
    theme(axis.title.x=element_blank(), axis.text.x=element_text(angle=45, vjust=1, hjust=1, size=16))
  p = ggplot(ccs_df, aes(y=crosscorrelation_rellag, x=tlx_group)) +
    geom_hline(yintercept=0, linetype="dashed") +
    geom_boxplot(aes(fill=tlx_group), outlier.shape = NA, outlier.alpha = 0) +
    # geom_text(aes(x=tlx_group, label=gsub("MACS3_", "", rdc_cluster), color=rdc_chrom, size=n_sense+n_anti), position=position_jitter(width=0.2)) +
    geom_point(aes(x=tlx_group, size=n_sense+n_anti), position=position_jitter(width=0.2), alpha=0.3) +
    labs(y="Relative distance between centromeric(+) and telomeric(-)\nstrand breaks peaks (cross-correlation)", size="No. of breaks", fill="APH concentration") +
    scale_fill_manual(values=group_palette) +
    # scale_alpha_manual(values=group_palette) +
    theme_bw(base_size=12) +
    theme(axis.title.x=element_blank(), axis.text.x=element_text(angle=45, vjust=1, hjust=1, size=16))
  print(p)
  p$data = ccs_df %>% dplyr::filter(tlx_group!="DMSO")
  print(p)
  ggplot(tlx_shift_df, aes(y=tlx_strand_shift, x=tlx_group)) +
    geom_boxplot(aes(fill=tlx_group), outlier.shape = NA, outlier.alpha = 0) +
    geom_point(aes(x=tlx_group, size=n_sense+n_anti), position=position_jitter(width=0.2), alpha=0.3) +
    # geom_text(aes(x=tlx_group, label=gsub("MACS3_", "", rdc_cluster), size=n_sense+n_anti), color="#FF0000", position=position_jitter(width=0.2)) +
    labs(y="Relative distance between centromeric(+) and telomeric(-)\nstrand breaks mean position", size="No. of breaks", fill="APH concentration") +
    scale_fill_manual(values=group_palette) +
    theme_bw(base_size=12) +
    theme(axis.title.x=element_blank(), axis.text.x=element_text(angle=45, vjust=1, hjust=1, size=16))
  dev.off()
}