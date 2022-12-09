library(dplyr)
devtools::load_all('breaktools/')

common_motif = function()
{
  rdc_df = readr::read_tsv("data/rdc.tsv") %>%
    dplyr::filter(tlx_group=="APH-Inter" & rdc_subset=="Wei+DKFZ" & rdc_is_significant) %>%
    dplyr::select(-tlx_group, -rdc_subset) %>%
    dplyr::do(as.data.frame(get_seq("genomes/mm10/mm10.fa", df2ranges(., rdc_chrom, rdc_start, rdc_end)))) %>%
    dplyr::select(dplyr::matches("rdc_"), rdc_sequence=sequence)
  control_df = rdc_df %>%
    dplyr::select(-rdc_sequence) %>%
    tidyr::crossing(data.frame(offset=c(-1, 1))) %>%
    dplyr::mutate(rdc_start=rdc_start+rdc_length*offset) %>%
    dplyr::do(as.data.frame(get_seq("genomes/mm10/mm10.fa", df2ranges(., rdc_chrom, rdc_start, rdc_end)))) %>%
    dplyr::select(dplyr::matches("rdc_"), rdc_sequence=sequence)


  output = "output1"
  dir.create(output, recursive=T, showWarnings=F)
  writeLines(paste0(">", rdc_df$rdc_name, "\n", rdc_df$rdc_sequence), file.path(output, "input.fasta"))
  writeLines(paste0(">", control_df$rdc_name, "\n", control_df$rdc_sequence), file.path(output, "control.fasta"))
  writeLines(c("#!/bin/bash", paste0("meme -objfun se -searchsize ", max(rdc_df$rdc_length+1), " -dna -nmotifs 100 -revcomp input.fasta -neg control.fasta -oc results")),file.path(output, "run.sh"))
}