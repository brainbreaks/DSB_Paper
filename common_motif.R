library(dplyr)
devtools::load_all('breaktools/')

common_motif = function()
{
  rdc_df = readr::read_tsv("~/Workspace/Datasets/HTGTS/rdc_pnas_mm10.tsv") %>%
    dplyr::mutate(rdc_cluster=paste0(rdc_chrom, ":", rdc_cluster, rdc_strand)) %>%
    dplyr::do(as.data.frame(get_seq("~/Workspace/genomes/mm10/mm10.fa", df2ranges(., rdc_chrom, rdc_start, rdc_end)))) %>%
    dplyr::select(dplyr::matches("rdc_"), rdc_sequence=sequence)


  writeLines(paste0(">", rdc_df$rdc_cluster, "\n", rdc_df$rdc_sequence), "rdc.fasta")

  colnames(rdc_df)

  rdc_ranges = df2ranges(rdc_df, rdc_chrom, rdc_start, rdc_end)
  x = get_seq("~/Workspace/genomes/mm10/mm10.fa", rdc_ranges)
  writeLines(x)
}