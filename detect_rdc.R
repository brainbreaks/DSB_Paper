setwd("~/Workspace/DSB_Paper")

library(readr)
library(dplyr)
library(ggplot2)
library(cowplot)
library(stringr)
library(UpSetR)
library(grid)
devtools::load_all("~/Workspace/breaktools/")



detect_rdc = function()
{
  debug=F
  dir.create("reports/detect_rdc", recursive=T, showWarnings=F)

  #
  # Load offtargets
  #
  offtargets_df = readr::read_tsv("data/offtargets_dkfz.tsv")


  #
  # Load gene annotations
  #
  genes_df = gtf_read('~/Workspace/genomes/mm10/annotation/mm10.ncbiRefSeq.gtf.gz')
  genes_ranges = genes_df %>%
    dplyr::filter(gene_length>=1e5) %>%
    df2ranges(gene_chrom, gene_start, gene_end)

  #
  # Load RDC
  #
  samples_df = tlx_read_samples("~/Workspace/Datasets/HTGTS/samples/All_samples.tsv", "~/Workspace/Datasets/HTGTS") %>%
    dplyr::filter(tlx_exists & celltype=="NPC" & organism=="mouse" & sample!="VI035" & (
      grepl("(Csmd1|Ctnna2|Nrxn1) promoter/enhancer", experiment) |
      grepl("concentration", experiment) & concentration==0.4 |
      grepl("Wei", experiment))
    )

  tlx_all_df = tlx_read_many(samples_df, threads=10) %>%
    tlx_extract_bait(bait_size=19, bait_region=12e6) %>%
    tlx_calc_copynumber(bowtie2_index="~/Workspace/genomes/mm10/mm10", max_hits=100, threads=24) %>%
    tlx_mark_offtargets(offtargets_df, offtarget_region=1e6, bait_region=1e4)
  libfactors_df = tlx_all_df %>%
    tlx_libsizes()

  #
  # Select junctions suitable for analysis
  #
  tlx_rdc_df = tlx_all_df %>%
    tlx_remove_rand_chromosomes() %>%
    dplyr::filter(tlx_copynumber==1 & !tlx_duplicated & !tlx_is_bait_junction & !tlx_is_offtarget & (!grepl("prom", group) | !tlx_is_bait_chrom)) %>%
    dplyr::mutate(tlx_group=dplyr::case_when(
      !tlx_control & tlx_is_bait_chrom ~ "APH-Intra",
      !tlx_control & !tlx_is_bait_chrom ~ "APH-Inter",
      tlx_control & tlx_is_bait_chrom ~ "DMSO-Intra",
      tlx_control & !tlx_is_bait_chrom ~ "DMSO-Inter",
    ), tlx_control=F)


  #
  # Evaluate extsize
  #
  if(F) {
    extsize_df = data.frame(extsize=c(seq(100, 500, 100), seq(1000, 9000, 1000), seq(10000, 90000, 10000), seq(1e5, 3e5, 5e5))) %>%
      dplyr::rowwise() %>%
      dplyr::do((function(z){
        zz<<-z
        params_extsize = macs2_params(extsize=z$extsize, exttype="symmetrical", llocal=1e7, minqvalue=0.01, effective_size=1.87e9, maxgap=2e5, minlen=2e5, baseline=2)
        tlxcov_rdc_df = tlx_rdc_df %>%
          dplyr::filter(!tlx_is_bait_junction & (!grepl("prom", group) | !tlx_is_bait_chrom)) %>%
          tlx_coverage(group="group", extsize=params_extsize$extsize, exttype=params_extsize$exttype, libfactors_df=libfactors_df, ignore.strand=T)
        macs_rdc = tlxcov_macs2(tlxcov_rdc_df, group="group", params=params_extsize)
        macs_rdc$islands %>%
          tidyr::crossing(as.data.frame(z))
      })(.)) %>%
      dplyr::ungroup()
    extsize_sumdf = extsize_df %>%
      dplyr::group_by(extsize, tlx_group) %>%
      dplyr::summarize(count=dplyr::n(), width=mean(island_end-island_start))
    pdf("reports/extsize_selection.pdf", width=11.69, height=8.27, paper="a4r")
    ggplot(extsize_sumdf) +
      geom_line(aes(x=extsize, y=count)) +
      geom_vline(xintercept=5e4, color="#FF0000") +
      geom_smooth(aes(x=extsize, y=count)) +
      labs(y="RDC count") +
      facet_wrap(~tlx_group, scales="free")
    ggplot(extsize_sumdf) +
      geom_line(aes(x=extsize, y=width)) +
      geom_vline(xintercept=5e4, color="#FF0000") +
      geom_smooth(aes(x=extsize, y=width)) +
      labs(y="RDC width") +
      facet_wrap(~tlx_group, scales="free")
    dev.off()
  }


  #
  # Detect RDC
  #
  params_rdc = macs2_params(extsize=5e4, exttype="symmetrical", llocal=1e7, minqvalue=0.01, effective_size=1.87e9, maxgap=250e3, minlen=200e3, baseline=2)
  tlxcov_rdc_df = tlx_rdc_df %>%
    tlx_coverage(group="group", extsize=params_rdc$extsize, exttype=params_rdc$exttype, libfactors_df=libfactors_df, ignore.strand=T)
  macs_rdc = tlxcov_macs2(tlxcov_rdc_df, group="group", params=params_rdc)

  #
  # Write debugging information from RDC calling
  #
  if(debug)
  {
    macs_rdc$qvalues %>%
      dplyr::group_by(tlx_group) %>%
      dplyr::do((function(df){
        df %>%
          dplyr::select(qvalue_chrom, qvalue_start, qvalue_end, qvalue_score) %>%
          readr::write_tsv(paste0("reports/detect_rdc/qvalues-", df$tlx_group[1], ".bedgraph"), col_names = F)
      })(.))
    macs_rdc$islands %>%
      dplyr::mutate(strand="*") %>%
      dplyr::group_by(tlx_group) %>%
      dplyr::do((function(df){
        df %>%
          dplyr::select(island_chrom, island_extended_start, island_extended_end, island_name, island_baseline, strand) %>%
          readr::write_tsv(paste0("reports/detect_rdc/extislands-", df$tlx_group[1], ".bed"), col_names=F)
        df %>%
          dplyr::mutate(score=1) %>%
          dplyr::select(island_chrom, island_start, island_end, island_name, score, strand) %>%
          readr::write_tsv(paste0("reports/detect_rdc/islands-", df$tlx_group[1], ".bed"), col_names=F)
      })(.))

    offtargets_df %>%
      dplyr::mutate(score=1, strand="*") %>%
      dplyr::select(offtarget_chrom, offtarget_start, offtarget_end, offtarget_bait_chrom, score, strand) %>%
      readr::write_tsv(paste0("reports/detect_rdc/offtargets_dkfz.bed"), col_names=F)


    tlx_all_df %>%
      tlx_remove_rand_chromosomes() %>%
      dplyr::filter(tlx_copynumber==1 & !tlx_duplicated & !tlx_is_bait_junction & tlx_is_offtarget & (!grepl("prom", group) | !tlx_is_bait_chrom)) %>%
      dplyr::mutate(tlx_group=ifelse(tlx_is_bait_chrom, "Intra", "Inter")) %>%
      tlx_coverage(group="group", extsize=params_rdc$extsize, exttype=params_rdc$exttype, libfactors_df=libfactors_df, ignore.strand=T) %>%
      tlxcov_write_bedgraph(path="reports/detect_rdc/offtargets", group="group")
    tlxcov_rdc_df %>% tlxcov_write_bedgraph(path="reports/detect_rdc/bedgraph", group="group")
    tlx_rdc_df %>% tlx_write_bed(path="reports/detect_rdc/bed", group="group")

    tlx_rdc_df %>%
      tlx_coverage(group="group", extsize=params_rdc$extsize, exttype=params_rdc$exttype, libfactors_df=libfactors_df, ignore.strand=F) %>%
      tlxcov_write_bedgraph(path="reports/pdetect_rdc/bedgraph", group="group")
  }


  #
  # Create final RDC table
  #
  rdc_df = macs_rdc$islands %>%
    df2ranges(island_chrom, island_start, island_end) %>%
    leftJoinByOverlaps(genes_ranges) %>%
    dplyr::arrange(dplyr::desc(tidyr::replace_na(gene_length, 0))) %>%
    dplyr::rename(rdc_chrom="island_chrom", rdc_start="island_start", rdc_end="island_end", rdc_extended_start="island_extended_start", rdc_extended_end="island_extended_end", rdc_stn="island_snr", rdc_qvalue_log10="island_summit_qvalue", rdc_pileup="island_summit_abs") %>%
    dplyr::mutate(rdc_cluster=paste0("RDC_", stringr::str_pad((0:(dplyr::n()))[-1], 3, pad="0"))) %>%
    dplyr::group_by(tlx_group, rdc_cluster, rdc_chrom, rdc_start, rdc_end, rdc_extended_start, rdc_extended_end, rdc_stn, rdc_pileup, rdc_qvalue_log10) %>%
    dplyr::summarize(
      genes_count=sum(!is.na(gene_id)),
      rdc_strand=ifelse(length(unique(na.omit(gene_strand)))==1, unique(na.omit(gene_strand)), "."),
      rdc_gene=ifelse(all(is.na(gene_id)), "-", paste0(na.omit(unique(gene_id))[pmin(5, length(na.omit(unique(gene_id))))], collapse=",")),
      rdc_group=dplyr::case_when(genes_count==1~1, genes_count>1~2, T~3),
      rdc_cluster_display=paste0(rdc_chrom, ":", rdc_cluster, rdc_strand)[1]) %>%
    dplyr::select(tlx_group, rdc_chrom, rdc_start, rdc_end, rdc_extended_start, rdc_extended_end, rdc_cluster,	rdc_group, rdc_strand, rdc_gene, rdc_cluster_display, rdc_stn, rdc_qvalue_log10, rdc_pileup) %>%
    dplyr::ungroup()


  #
  # Overlap with published datasets
  #
  pubrdc_df = readr::read_tsv("data/pubrdc.tsv") %>%
    dplyr::filter(pubrdc_source %in% c("Wei2018", "Wei2018_DMSO", "Tena2020") & pubrdc_celline %in% c("NPC", "NSPC")) %>%
    tidyr::separate_rows(rdcpub_bait_chrom, rdcpub_bait_chrom, sep=", ?") %>%
    dplyr::mutate(tlx_group=dplyr::case_when(
      pubrdc_source=="Wei2018_DMSO" & rdcpub_bait_chrom==pubrdc_chrom ~ "DMSO-Intra",
      pubrdc_source=="Wei2018_DMSO" & rdcpub_bait_chrom!=pubrdc_chrom ~ "DMSO-Inter",
      rdcpub_bait_chrom==pubrdc_chrom ~ "APH-Intra",
      rdcpub_bait_chrom!=pubrdc_chrom ~ "APH-Inter")) %>%
    dplyr::filter(!is.na(pubrdc_start))
  table(pubrdc_df$tlx_group)

  pubrdc_reduced_df = pubrdc_df %>%
    dplyr::group_by(pubrdc_source, pubrdc_celline, tlx_group) %>%
    dplyr::do(GenomicRanges::reduce(df2ranges(., pubrdc_chrom, pubrdc_start, pubrdc_end)) %>% as.data.frame()) %>%
    dplyr::ungroup() %>%
    dplyr::select(tlx_group, pubrdc_source, pubrdc_celline, pubrdc_chrom=seqnames, pubrdc_start=start, pubrdc_end=end)


  allrdc_df = dplyr::bind_rows(
    pubrdc_reduced_df %>%
      dplyr::mutate(rdc_source=paste0(pubrdc_source, "-", pubrdc_celline), rdc_source_location=paste0(rdc_source, " (", tlx_group, ")")) %>%
      dplyr::select(tlx_group, rdc_chrom=pubrdc_chrom, rdc_start=pubrdc_start, rdc_end=pubrdc_end, rdc_source, rdc_source_location),
    rdc_df %>%
        dplyr::mutate(rdc_source="DKFZ-NPC", rdc_source_location=paste0(rdc_source, " (", tlx_group, ")")) %>%
        dplyr::select(tlx_group, rdc_chrom, rdc_start, rdc_end, rdc_source, rdc_source_location))
  allrdc_ranges = allrdc_df %>% df2ranges(rdc_chrom, rdc_start, rdc_end)
  allrdc_reduced_df = allrdc_ranges %>%
    GenomicRanges::reduce(min.gapwidth=500e3) %>%
    as.data.frame() %>%
    dplyr::mutate(rdc_reduced_id=paste0(seqnames, ":", start, "-", end)) %>%
    dplyr::select(rdc_reduced_chrom=seqnames, rdc_reduced_start=start, rdc_reduced_end=end, rdc_reduced_id)

  # allrdc_reduced_df %>%
  #   dplyr::select(rdc_reduced_chrom, rdc_reduced_start, rdc_reduced_end) %>%
  #   readr::write_tsv("reports/detect_rdc/allrdc_reduced.bed", col_names=F)
  #
  #   allrdc_df %>%
  #     dplyr::group_by(tlx_group) %>%
  #     dplyr::do((function(df){
  #       df %>%
  #         dplyr::select(rdc_chrom, rdc_start, rdc_end) %>%
  #         readr::write_tsv(paste0("reports/detect_rdc/allrdc_-", df$tlx_group[1], ".bed"), col_names=F)
  #     })(.))


  overlaps_df = allrdc_reduced_df %>%
    df2ranges(rdc_reduced_chrom, rdc_reduced_start, rdc_reduced_end) %>%
    innerJoinByOverlaps(allrdc_ranges) %>%
    reshape2::dcast(rdc_chrom+rdc_reduced_start+rdc_reduced_end+rdc_reduced_id+tlx_group ~ rdc_source_location, value.var="rdc_source_location", fun.aggregate=function(z) pmin(1, length(z)))


  pdf("reports/rdc_compare_with_published.pdf", width=8.27, height=8.27)
  overlaps_upset = split(overlaps_df, overlaps_df$tlx_group) %>%
    lapply(FUN=function(odf) {
      ooo <<- odf
      odf_datasets = colnames(odf)[colSums(odf==1)>0]
      UpSetR::upset(as.data.frame(odf), sets=odf_datasets, order.by=c("freq", "degree"), decreasing=c(T,F))
      # grid::grid.text(paste0("RDC (", odf$tlx_group[1], "-chromosomal)"), x=0.65, y=0.95, gp=grid::gpar(fontsize=20))
    })
  gridExtra::grid.arrange(grobs=overlaps_upset, ncol = 2)


  plot.new()
  p = VennDiagram::venn.diagram(
    x = list(
      DKFZ=overlaps_inter_df %>% dplyr::filter(`DKFZ-NPC`>0) %>% .$rdc_reduced_id,
      Wei2018_Tena2020=overlaps_inter_df[rowSums(overlaps_inter_df %>% dplyr::select(dplyr::matches("(Tena2020-NPC|Wei2018-NPC|Wei2018_DMSO-NPC)"))) > 0, "rdc_reduced_id"]),
    lwd = 2, lty = 'blank', fill = RColorBrewer::brewer.pal(3, "Pastel2")[1:2], cat.cex=2, cat.fontface = "bold", filename=NULL)
  grid::grid.draw(p)

  overlapping_baits = intersect(with(tlx_rdc_df, names(table(bait_chrom))[table(bait_chrom)>5000]), pubrdc_df$rdcpub_bait_chrom)
  overlaps_intra_df = overlaps_df %>% dplyr::filter(tlx_group=="Intra" & rdc_chrom %in% overlapping_baits)
  UpSetR::upset(overlaps_intra_df, sets=colnames(overlaps_intra_df)[colSums(overlaps_intra_df==1)>0], order.by=c("freq", "degree"), decreasing=c(T,F))
  grid::grid.text(paste0("RDC (Intra-chromosomal: ", paste(overlapping_baits, collapse=", "), ")"), x=0.65, y=0.95, gp=grid::gpar(fontsize=20))

  plot.new()
  p = VennDiagram::venn.diagram(
    x = list(
      DKFZ=overlaps_intra_df %>% dplyr::filter(`DKFZ-NPC`>0) %>% .$rdc_reduced_id,
      Wei2018_Tena2020=overlaps_intra_df[rowSums(overlaps_intra_df %>% dplyr::select(dplyr::matches("(Tena2020-NPC|Wei2018-NPC|Wei2018_DMSO-NPC)"))) > 0, "rdc_reduced_id"]),
    lwd = 2, lty = 'blank', fill = RColorBrewer::brewer.pal(3, "Pastel2")[1:2], cat.cex=2, cat.fontface = "bold", filename=NULL)
  grid::grid.draw(p)
  dev.off()


}