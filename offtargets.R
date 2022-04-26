library(readr)
library(dplyr)
library(ggplot2)
library(cowplot)
library(stringr)
# library(plyranges)
devtools::load_all('~/Workspace/breaktools/')

offtargets = function()
{
  params = macs2_params(extsize=1e5, exttype="symmetrical", llocal=1e7, minqvalue=0.01, effective_size=1.87e9, maxgap=2e5, minlen=1e5)

  # Load TLX
  samples_df = tlx_read_samples("~/Workspace/Datasets/HTGTS/samples/All_samples.tsv", "~/Workspace/Datasets/HTGTS/TLX") %>%
    # dplyr::filter(experiment %in% c("APH concentration", "Wei et al. PNAS 2018")) %>%
    # TODO: add "Chr6 and offtarget bait|Chr8 and offtarget baits" when they are sequenced
    # dplyr::filter(grepl("Csmd1 promoter/enhancer|Ctnna2 promoter/enhancer|Nrxn1 promoter/enhancer|Wei et al", experiment) & alleles==2 | grepl("concentration", experiment) & concentration==0.4) %>%
    # dplyr::filter(grepl("Chr6 and offtarget bait|Chr8 and offtarget baits", experiment) & alleles==2) %>%
    dplyr::filter(!grepl("Sonic hedgehog|Hydroxyurea|Nocodazole", experiment) & !control) %>%
    dplyr::mutate(group=dplyr::case_when(
      control ~ "DMSO",
      grepl("Wei", experiment) & !control ~ "Treatment",
      grepl("Wei", experiment) & control ~ "DMSO",
      grepl("All", experiment) & control ~ "DMSO",
      T ~ gsub(" ?\\(.*", "", group)
    ))
  tlx_df = tlx_read_many(samples_df %>% dplyr::filter(tlx_exists), threads=30)
  tlx_df = tlx_extract_bait(tlx_df, bait_size=19, bait_region=2e6)
  tlx_df = tlx_df %>%
    dplyr::mutate(tlx_group=tlx_bait_chrom) %>%
    dplyr::group_by(tlx_sample) %>%
    dplyr::filter(dplyr::n() >= 2000) %>%
    dplyr::ungroup()

  libfactors_df = tlx_libfactors(tlx_df, group="group", normalize_within="group", normalize_between="none", normalization_target="smallest")
  tlxcov_df = tlx_df %>%
    dplyr::filter(!tlx_is_bait_junction) %>%
    tlx_coverage(group="group", exttype=params$exttype, extsize=params$extsize, libfactors_df=libfactors_df, ignore.strand=F)
  tlxcov_ranges = tlxcov_df %>% df2ranges(tlxcov_chrom, tlxcov_start, tlxcov_end, tlx_strand)

  offtarget_df = readr::read_tsv("~/Workspace/Datasets/HTGTS/offtargets_pnas_mm10.tsv")
  offtarget_ranges = offtarget_df %>%
    dplyr::rename(offtarget_bait_chrom="bait_chrom") %>%
    df2ranges(offtarget_chrom, offtarget_start-5e5, offtarget_end+5e5)


  pdf("reports/offtargets.pdf", width=6*11.69, height=3*11.69)
  tlxcov2offtargets_df = innerJoinByOverlaps(tlxcov_ranges, offtarget_ranges) %>%
    dplyr::filter(tlx_group==offtarget_bait_chrom) %>%
    dplyr::mutate(rdc_chrom=tlxcov_chrom, rdc_cluster=paste0("bait-", offtarget_bait_chrom, " off-", offtarget_chrom, "[", offtarget_start, "-", offtarget_end, "]"), rdc_cluster_display=rdc_cluster)
  ggplot() +
    geom_tlxcov(tlxcov2offtargets_df) +
    facet_wrap(~rdc_cluster, scales="free")
  dev.off()

  genes_df = gtf_read('~/Workspace/genomes/mm10/annotation/mm10.ncbiRefSeq.gtf.gz')
  ctnna2_ranges = genes_df %>%
    dplyr::filter(gene_id=="Ctnna2") %>%
    df2ranges(gene_chrom, gene_start-5e5, gene_end+5e5)
  pdf("reports/Ctnna2.pdf", width=5.69, height=3.69)
  tlxcov2ctnna2_df = innerJoinByOverlaps(tlxcov_ranges, ctnna2_ranges) %>%
    dplyr::mutate(rdc_chrom=tlxcov_chrom, rdc_cluster="Ctnna2", rdc_cluster_display=rdc_cluster)
  ggplot() +
    geom_tlxcov(tlxcov2ctnna2_df) +
    facet_wrap(~rdc_cluster, scales="free") +
    coord_cartesian(xlim=c(76.5e6, 78.5e6))
  dev.off()
}