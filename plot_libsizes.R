setwd("~/Workspace/Everything")

library(readr)
library(dplyr)
library(ggplot2)
library(cowplot)
library(stringr)
# library(plyranges)
devtools::load_all('~/Workspace/breaktools/')


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

plot_libsizes = function()
{
  # Load TLX
  samples_df = tlx_read_samples("~/Workspace/Datasets/HTGTS/samples/All_samples.tsv", "~/Workspace/Datasets/HTGTS/TLX") %>%
    # dplyr::filter(experiment %in% c("APH concentration", "Wei et al. PNAS 2018")) %>%
    # TODO: add "Chr6 and offtarget bait|Chr8 and offtarget baits" when they are sequenced
    # dplyr::filter(grepl("Csmd1 promoter/enhancer|Ctnna2 promoter/enhancer|Nrxn1 promoter/enhancer|Wei et al", experiment) & alleles==2 | grepl("concentration", experiment) & concentration==0.4) %>%
    # dplyr::filter(grepl("Chr6 and offtarget bait|Chr8 and offtarget baits", experiment) & alleles==2) %>%
    dplyr::filter(!grepl("Sonic hedgehog|Hydroxyurea|Nocodazole", experiment)) %>%
    dplyr::select(-chrom) %>%
    dplyr::mutate(group=dplyr::case_when(
      control ~ "DMSO",
      grepl("Wei", experiment) & !control ~ "Treatment",
      grepl("Wei", experiment) & control ~ "DMSO",
      grepl("All", experiment) & control ~ "DMSO",
      T ~ gsub(" ?\\(.*", "", group)
    ))

  samples_df.missing = samples_df %>%
    dplyr::filter(!tlx_exists) %>%
    dplyr::group_by(experiment) %>%
    dplyr::summarise(samples=paste0(sample, collapse=", "))

  writeLines(paste0("TLX files missing: ", paste0(samples_df.missing$experiment, ":\n    ", samples_df.missing$samples, collapse="\n")))

  tlx_df = tlx_read_many(samples_df %>% dplyr::filter(tlx_exists), threads=30)
  tlx_df = tlx_extract_bait(tlx_df, bait_size=19, bait_region=2e6)

  # Display library sizes
  libsizes_df = tlx_df %>%
    dplyr::group_by(tlx_sample, .drop=F) %>%
    dplyr::summarize(library_size=dplyr::n(), usefull_size=sum(!tlx_is_bait_junction & tlx_is_bait_chrom)) %>%
    dplyr::inner_join(samples_df, by=c("tlx_sample"="sample"))




  plist = lapply(split(libsizes_df, f=libsizes_df$experiment), FUN=function(df) {
    dff <<- df
    palette = c(
      "Parental cell"="#FCBBA1", "Allelic deletion"="#EF3B2C", "Allelic+promoter deletion"="#67000D",
      "APH 0.2 uM 96h"="#C6DBEF", "APH 0.4 uM 96h"="#4292C6", "APH 0.6 uM 96h"="#08306B",
      "Treatment"="#4292C6", "DMSO"="#666666")[unique(df$group)]

    ggplot(df) +
      geom_bar(aes(x=tlx_sample, y=library_size, fill=group, group=tlx_sample), position="dodge", color="#EEEEEE", size=0.1, stat="identity") +
      scale_y_continuous(labels=scales::label_number(accuracy=1, scale=1e-3, suffix="K")) +
      labs(y="Library size", fill="Group") +
      facet_grid(~experiment, scales="free", space="free_x") +
      theme_grey(base_size=10) +
      theme(axis.text.x = ggplot2::element_text(angle=45, hjust=1, vjust=1, size=4), axis.title.x=ggplot2::element_blank(), legend.position=ifelse(!grepl("Csmd1|Ctnna2", df$experiment[1]), "bottom", "none")) +
      scale_fill_manual(values=palette)
  })

  pdf("reports/libsizes.pdf", width=8.27, height=11.69)
  cowplot::plot_grid(plotlist=plist, align="v", axis="t", ncol=1)
  dev.off()

}