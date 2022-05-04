bedgraph = function()
{
  samples_df = tlx_read_samples("~/Workspace/Datasets/HTGTS/samples/All_samples.tsv", "~/Workspace/Datasets/HTGTS") %>%
    dplyr::filter(celltype=="NPC" & tlx_exists) %>%
    dplyr::filter(grepl("promoter/enhancer", experiment) & alleles==2 | "APH concentration"==experiment & concentration==0.4 | grepl("Wei", experiment) | grepl("Tena", experiment)) %>%
    dplyr::mutate(group=paste0("APH (", bait_chrom, ")"))

  tlx_all_df = tlx_read_many(samples_df, threads=30) %>%
    tlx_remove_rand_chromosomes() %>%
    tlx_extract_bait(bait_size=19, bait_region=2e6) %>%
    tlx_mark_dust() %>%
    tlx_calc_copynumber(bowtie2_index="~/Workspace/genomes/mm10/mm10", max_hits=100, threads=24)

  used_samples = tlx_all_df %>%
    dplyr::filter(tlx_copynumber==1 & !tlx_duplicated) %>%
    dplyr::group_by(tlx_group, tlx_control, tlx_sample) %>%
    dplyr::summarize(tlx_count=dplyr::n()) %>%
    dplyr::group_by(tlx_group, tlx_control) %>%
    dplyr::filter(dplyr::n()==1 | tlx_count>=max(tlx_count)/2) %>%
    dplyr::ungroup() %>%
    .$tlx_sample
  tlx_used_df = tlx_all_df %>%
    dplyr::filter(tlx_sample %in% used_samples) %>%
    dplyr::filter(tlx_copynumber==1 & !tlx_duplicated)

  #
  # Get offtargets
  #
  libfactors_offtargets_df = tlx_libfactors(tlx_used_df, normalize_within="none", normalize_between="group", normalization_target="min")
  params_offtargets = macs2_params(extsize=50, exttype="opposite", llocal=1e7, minqvalue=0.05, effective_size=1.87e9, maxgap=1e3, minlen=10)
  tlx_offtarget_df = tlx_used_df %>%
    dplyr::filter(tlx_is_bait_chrom & !tlx_is_bait_junction) %>%
    dplyr::mutate(tlx_group=bait_chrom, tlx_control=F)
  tlxcov_offtargets_df = tlx_offtarget_df %>%
    tlx_coverage(group="group", extsize=params_offtargets$extsize, exttype=params_offtargets$exttype, libfactors_df=libfactors_offtargets_df, ignore.strand=T)
  macs_offtargets = tlxcov_macs2(tlxcov_offtargets_df, group="group", params_offtargets)
  macs_offtargets$qvalues %>% dplyr::select(qvalue_chrom, qvalue_start, qvalue_end, qvalue_score) %>% readr::write_tsv("reports/all/offtargets.bedgraph", col_names=F)
  macs_offtargets$islands %>% df2ranges(island_chrom, island_start, island_end) %>% rtracklayer::export.bed("reports/all/offtargets.bed")

  #
  # Mark offtargets
  #

  baits_df = readr::read_tsv("~/Workspace/Datasets/HTGTS/wei_pnas2018_baits.tsv")
  offtargets_df = baits_df %>%
    dplyr::left_join(macs_offtargets$islands, by=c("bait_chrom"="island_chrom")) %>%
    dplyr::mutate(offtarget_is_primary=F, offtarget_start=ifelse(is.na(island_start), bait_start, island_start), offtarget_end=ifelse(is.na(bait_end), bait_end, island_end)) %>%
    dplyr::select(
        offtarget_bait_chrom=bait_chrom, offtarget_bait_start=bait_start, offtarget_bait_end=bait_end, offtarget_bait_strand=bait_strand,
        offtarget_chrom=bait_chrom, offtarget_start, offtarget_end, offtarget_is_primary)


  devtools::load_all('~/Workspace/breaktools/')
  libfactors_df = tlx_libfactors(tlx_used_df, normalize_within="treatment", normalize_between="treatment", normalization_target="min") %>%
    dplyr::arrange(tlx_group, tlx_control) %>%
    data.frame()
  tlx_clean_df = tlx_used_df %>%
    tlx_mark_offtargets(offtargets_df, offtarget_region=1e5, bait_region=1e5) %>%
    dplyr::filter(!tlx_is_offtarget & !tlx_is_bait_junction & tlx_is_bait_chrom)


  table(tlx_clean_df %>% dplyr::filter(tlx_is_offtarget & Junction>=138172000 & Junction<=138178000) %>% .$tlx_control, tlx_clean_df %>% dplyr::filter(tlx_is_offtarget & Junction>=138172000 & Junction<=138178000) %>% .$tlx_strand)
  # offtargets_df = offtargets_df %>% dplyr::filter(offtarget_bait_chrom=="chr5")
  # tlx_df = tlx_clean_df %>% dplyr::filter(Rname=="chr5")
  tlx_write_bed(tlx_clean_df, path="reports/all/all", group="all", mode="alignment", ignore.strand=T)
  tlx_clean_df %>% dplyr::filter(Rname=="chr5" & tlx_is_offtarget & !tlx_is_bait_junction & tlx_is_bait_chrom)

  params_clean = macs2_params(extsize=5e4, exttype="symmetrical", llocal=1e7, minqvalue=0.01, effective_size=1.87e9, maxgap=2e5, minlen=1e5)
  tlxcov_clean_df = tlx_clean_df %>%
    tlx_coverage(group="group", extsize=params_clean$extsize, exttype=params_clean$exttype, libfactors_df=libfactors_df, ignore.strand=F)
  tlxcov_write_bedgraph(tlxcov_clean_df, path="reports/all/all", group="all")
  tlx_write_bed(tlx_clean_df, path="reports/all/all", group="all", mode="junction", ignore.strand=F)
}