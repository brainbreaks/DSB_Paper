crossbait_contamination = function() {
  samples_df = tlx_read_samples("~/Workspace/Datasets/HTGTS/samples/All_samples.tsv", "~/Workspace/Datasets/HTGTS/TLX") %>%
    dplyr::filter(grepl("concentration", experiment)) %>%
    dplyr::mutate(run=dplyr::case_when(
      grepl("PW103|PW104|PW105|PW106|PW107|PW108|PW109|PW110|PW111|PW112|PW113|PW114", sample) ~ "B400_017",
      grepl("PW115|PW116|PW117|PW118|PW119|PW120|PW121|PW122|PW123|PW124|PW125|PW126|PW127|PW128|PW129|PW130|PW131|PW132|PW133|PW134", sample) ~ "B400_018",
      T ~ "other"
    )) %>%
    dplyr::filter(grepl("B400", run)) %>%
    dplyr::mutate(fasta_r1=paste0("~/Workspace/", run, "/preprocess/", sample, "_R1.fq.gz"), fasta_r2=paste0("~/Workspace/", run, "/preprocess/", sample, "_R2.fq.gz"))
  samples_long_df = samples_df %>%
    reshape2::melt(measure.vars=c("fasta_r1", "fasta_r2"), value.name="fasta" )

  tlx_df = tlx_read_many(samples_df, threads=30)
  tlx_df = tlx_extract_bait(tlx_df, bait_size=19, bait_region=2e6)

  baits_df = tlx_identify_baits(tlx_df, genome_fasta="~/Workspace/genomes/mm10/mm10.fa") %>%
    dplyr::distinct(bait_chrom, .keep_all=T) %>%
    dplyr::select(-bait_group, -bait_sample)

  bait_presence_df = samples_long_df %>%
    dplyr::rowwise() %>%
    dplyr::do((function(z){
      zz<<-z
      fasta = ShortRead::readFastq(z$fasta)
      fasta_reads = ShortRead::sread(fasta)
      baits_df %>%
        tidyr::crossing(as.data.frame(z)) %>%
        dplyr::rowwise() %>%
        dplyr::mutate(complement=sum(Biostrings::vcountPattern(Biostrings::complement(Biostrings::DNAString(bait_sequence)), fasta_reads)>0), total=length(fasta_reads))  %>%
        dplyr::mutate(reverse=sum(Biostrings::vcountPattern(Biostrings::reverse(Biostrings::DNAString(bait_sequence)), fasta_reads)>0), total=length(fasta_reads))  %>%
        dplyr::mutate(reverse_complement=sum(Biostrings::vcountPattern(Biostrings::reverseComplement(Biostrings::DNAString(bait_sequence)), fasta_reads)>0), total=length(fasta_reads))  %>%
        dplyr::mutate(direct=sum(Biostrings::vcountPattern(bait_sequence, fasta_reads)>0), total=length(fasta_reads)) %>%
        dplyr::mutate(matches=complement+reverse+reverse_complement+direct) %>%
        dplyr::ungroup()
    })(.))

  bait_presence_dff = bait_presence_df %>%
    dplyr::group_by(sample, treatment, chrom, bait_chrom) %>%
    dplyr::summarise(matches=sum(matches)) %>%
    dplyr::mutate(expected=ifelse(bait_chrom==chrom, "expected bait", "unexpected bait"))

  ggplot(bait_presence_dff) +
    geom_bar(aes(x=paste(sample, treatment), y=matches, fill=bait_chrom), position="dodge", size=0.1, stat="identity") +
    facet_grid(chrom~expected, scales="free") +
    coord_flip()
}