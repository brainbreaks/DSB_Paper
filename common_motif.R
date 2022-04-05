library(dplyr)
devtools::load_all('~/Workspace/breaktools/')

common_motif = function()
{
  rdc_df = readr::read_tsv("~/Workspace/Datasets/HTGTS/rdc_macs2_mm10.tsv") %>%
    dplyr::mutate(rdc_width=rdc_end-rdc_start) %>%
    dplyr::mutate(rdc_cluster=paste0(rdc_chrom, ":", rdc_cluster, rdc_strand)) %>%
    dplyr::do(as.data.frame(get_seq("~/Workspace/genomes/mm10/mm10.fa", df2ranges(., rdc_chrom, rdc_start, rdc_end)))) %>%
    dplyr::select(dplyr::matches("rdc_"), rdc_sequence=sequence)
  control_df = rdc_df %>%
    dplyr::select(-rdc_sequence) %>%
    tidyr::crossing(data.frame(offset=c(-1, 1))) %>%
    dplyr::mutate(rdc_start=rdc_start+rdc_width*offset) %>%
    dplyr::mutate(rdc_cluster=paste0(rdc_chrom, ":", rdc_cluster, rdc_strand, "__", offset)) %>%
    dplyr::do(as.data.frame(get_seq("~/Workspace/genomes/mm10/mm10.fa", df2ranges(., rdc_chrom, rdc_start, rdc_end)))) %>%
    dplyr::select(dplyr::matches("rdc_"), rdc_sequence=sequence)


  output = "output1"
  dir.create(output, recursive="T")
  writeLines(paste0(">", rdc_df$rdc_cluster, "\n", rdc_df$rdc_sequence), file.path(output, "input.fasta"))
  writeLines(paste0(">", control_df$rdc_cluster, "\n", control_df$rdc_sequence), file.path(output, "control.fasta"))
  writeLines(c("#!/bin/bash", paste0("meme -objfun se -searchsize ", max(rdc_df$rdc_width)+1, " -dna -nmotifs 100 -revcomp input.fasta -neg control.fasta -oc results")), con=file.path(output, "run.sh"))
}