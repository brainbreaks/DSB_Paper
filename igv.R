

main = function()
{
  extention = 1e2

  # file = "~/Workspace/Datasets/ORM/allfiber_All.HeLa.async.gtf"
  # file = "~/Workspace/Datasets/ORM/allfiber_C.0.gtf"
  for(file in paste0("~/Workspace/Datasets/ORM/", c("D.20.gtf", "D.30.gtf", "D.45.gtf", "D.60.gtf", "D.90.gtf"))) {
    data_df = readr::read_tsv(file, col_names=c("chrom", "pipeline", "type", "start", "end")) %>%
      dplyr::mutate(name="", strand="*", score="1")
    data_df$duplicated = duplicated(data_df[,c("chrom", "start", "end")])
    data_df = data_df %>%
      dplyr::mutate(transcript_id=gsub(".*transcript_id ", "", X9)) %>%
      dplyr::group_by(transcript_id) %>%
      dplyr::filter(all(!duplicated)) %>%
      dplyr::ungroup()

    tracks_df = data_df %>%
      dplyr::filter(type=="exon") %>%
      dplyr::group_by(chrom) %>%
      dplyr::mutate(length=max(end)) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(start=pmax(1, start-extention), end=pmin(end+extention, length))
    tracks_ranges = GenomicRanges::makeGRangesFromDataFrame(tracks_df, keep.extra.columns=T)
    tracks_coverage_df = as.data.frame(as(GenomicRanges::coverage(tracks_ranges), "GRanges")) %>%
      dplyr::rename(track_coverage="score") %>%
      dplyr::mutate(track_id=1:n(), track_chrom=seqnames, track_start=start, track_end=end)
    tracks_coverage_ranges = GenomicRanges::makeGRangesFromDataFrame(tracks_coverage_df, keep.extra.columns=T)
    tracks_coverage_path = gsub("\\.gtf", paste0(".", extention, ".track_coverage.bedGraph"), file)
    tracks_coverage_df %>%
      dplyr::mutate(name="") %>%
      dplyr::select(track_chrom, track_start, track_end, track_coverage) %>%
      dplyr::mutate(track_coverage=ifelse(is.infinite(track_coverage), 0, track_coverage)) %>%
      readr::write_tsv(file=tracks_coverage_path, col_names=F)

    fiber_df = data_df %>%
      dplyr::filter(type=="transcript") %>%
      dplyr::group_by(chrom) %>%
      dplyr::mutate(length=max(end)) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(start=pmax(1, start-extention), end=pmin(end+extention, length))
    fiber_ranges = GenomicRanges::makeGRangesFromDataFrame(fiber_df, keep.extra.columns=T)
    fiber_coverage_df = as.data.frame(as(GenomicRanges::coverage(fiber_ranges), "GRanges")) %>%
      dplyr::rename(fiber_coverage="score") %>%
      dplyr::mutate(fiber_id=1:n(), fiber_chrom=seqnames, fiber_start=start, fiber_end=end)
    fiber_coverage_ranges = GenomicRanges::makeGRangesFromDataFrame(fiber_coverage_df, keep.extra.columns=T)

    effectivenesss_df = as.data.frame(IRanges::mergeByOverlaps(tracks_coverage_ranges, fiber_coverage_ranges))  %>%
      dplyr::select(-dplyr::matches("_ranges\\.")) %>%
      dplyr::mutate(effectiveness=track_coverage/fiber_coverage)
    effectivenesss_df %>%
      dplyr::mutate(name="", effectiveness=tidyr::replace_na(effectiveness, 0)) %>%
      dplyr::select(fiber_chrom, fiber_start, fiber_end, effectiveness) %>%
      dplyr::mutate(effectiveness=ifelse(is.infinite(effectiveness), 0, effectiveness)) %>%
      readr::write_tsv(file=gsub("\\.gtf", paste0(".", extention, ".effectiveness.bedGraph"), file), col_names=F)
  }





    # dplyr::filter(end-start>0) %>%
   tracks_df %>%
    dplyr::select(chrom, start, end, name, score, strand) %>%
    readr::write_tsv(file=gsub("\\.gtf", ".bed", file), col_names=F)

  tracks_ranges = tracks_df


  fiber_ranges = fiber_df %>%
    dplyr::group_by(chrom) %>%
    dplyr::mutate(length=max(end)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(start=pmax(1, start-extention), end=pmin(end+extention, length)) %>%
    GenomicRanges::makeGRangesFromDataFrame()

  GenomicRanges::merge

as.data.frame()
  hela_coverage %>%
    dplyr::mutate(name="") %>%
    dplyr::select(seqnames, start, end, score) %>%
    readr::write_tsv(file=gsub("\\.gtf", paste0(".", extention, ".2.bedGraph"), file), col_names=F)


  tlx_df = tlx_read_many(samples_df, threads=30)
  tlx_df = tlx_remove_rand_chromosomes(tlx_df)
  tlx_df = tlx_extract_bait(tlx_df, bait_size=19, bait_region=2e6)
  tlx_df = tlx_df %>% dplyr::inner_join(samples_df, by=c("tlx_sample"="sample"))
  tlx_df = tlx_mark_dust(tlx_df)
  # tlx_df = tlx_df %>% dplyr::filter(!tlx_is_bait_junction & tlx_is_bait_chromosome)

  #
  # Export bedgraph
  #
  tlx_df %>%
    dplyr::filter(tlx_bait_chrom %in% chromosomes) %>%
    dplyr::group_by(tlx_sample) %>%
    dplyr::filter(dplyr::n() >= 2000) %>%
    dplyr::ungroup() %>%
    dplyr::filter(!tlx_control) %>%
    dplyr::mutate(tlx_group="APH") %>%
    dplyr::mutate(dbscan_cluster=dbscan::dbscan(matrix(Junction), minPts=20, eps=100)$cluster) %>%
    dplyr::filter(dbscan_cluster==0) %>%
    tlx_write_bedgraph(path="reports/bedgraph-1e5", group="group", exttype="symmetrical", extsize=1e5)


  readr::read_tsv("~/Workspace/Datasets/HTGTS/rdc_pnas_mm10.tsv") %>%
    dplyr::mutate(score=1) %>%
    dplyr::select(rdc_chrom, rdc_start,  rdc_end, rdc_cluster, score, rdc_strand) %>%
    readr::write_tsv("~/Workspace/Datasets/HTGTS/rdc_pnas_mm10.bed", col_names=F)

  readr::read_tsv("~/Workspace/Datasets/zhao_bmc_repliseq_2020/preprocessed/repliseq_NPC.tsv") %>%
    dplyr::arrange(repliseq_fraction, repliseq_chrom, repliseq_start) %>%
    dplyr::mutate(repliseq_fraction=paste0("repliseq_", repliseq_fraction)) %>%
    dplyr::mutate(repliseq_fraction=factor(repliseq_fraction, unique(repliseq_fraction))) %>%
    dplyr::mutate(repliseq_dtype="repliseq") %>%
    reshape2::dcast(repliseq_chrom+repliseq_start+repliseq_end+repliseq_dtype ~ repliseq_fraction, value.var="repliseq_value") %>%
    dplyr::rename(`#type=GENE_EXPRESSION\nrepliseq_chrom`="repliseq_chrom") %>%
    readr::write_tsv("~/Workspace/Datasets/zhao_bmc_repliseq_2020/preprocessed/repliseq_NPC.igv")

}