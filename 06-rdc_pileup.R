library(readr)
library(dplyr)
library(ggplot2)
library(cowplot)
library(stringr)
devtools::load_all('~/Workspace/breaktools/')


main = function()
{
  #
  # Load offtargets
  #
  offtargets_df = readr::read_tsv("data/offtargets_dkfz.tsv")

  #
  # Load genes
  #
  genes_df = gtf_read('genomes/mm10/annotation/mm10.ncbiRefSeq.gtf.gz') %>%
    dplyr::filter(gene_cluster_i<=3)
  genes_ranges = genes_df %>% df2ranges(gene_chrom, gene_start, gene_end)

  #
  # Replication and replication termination site (annotated repli-seq data)
  #
  replication_df = readr::read_tsv("data/replication_reduced_subsets.tsv")
  replication_ranges = replication_df %>% df2ranges(replication_chrom, pmin(replication_start, replication_end), pmax(replication_start, replication_end))
  rts_df = replication_df %>%
    dplyr::filter(replication_end>replication_start) %>%
    dplyr::select(rts_chrom=replication_chrom, rts_pos=replication_end)


  #
  # Load RDC and add the most middle termination site position
  #
  rdc_df = readr::read_tsv("data/rdc.tsv") %>%
    dplyr::filter(tlx_group=="APH-Inter" & rdc_subset=="Wei+DKFZ" & rdc_is_significant) %>%
    dplyr::mutate(rdc_region_start=rdc_extended_start-rdc_extended_length/4, rdc_region_end=rdc_extended_end+rdc_extended_length/4) %>%
    df2ranges(rdc_chrom, rdc_start, rdc_end) %>%
    leftJoinByOverlaps(rts_df %>% df2ranges(rts_chrom, rts_pos, rts_pos))  %>%
    dplyr::arrange(dplyr::desc(pmax(abs(rdc_extended_start-rts_pos), abs(rdc_extended_end-rts_pos)))) %>%
    dplyr::distinct(tlx_group, rdc_subset, rdc_name, .keep_all=T) %>%
    dplyr::select(dplyr::matches("rdc"), rts_pos, -rdc_subset)

  #
  # Load TLX
  #
  samples_df = tlx_read_samples("data/htgts_samples.tsv", "~/Workspace/Datasets/HTGTS") %>%
    dplyr::filter(!control & (grepl("(Nrxn1|Ctnna2) promoter/enhancer", experiment) | grepl("APH concentration", experiment) & concentration==0.4 | grepl("Wei", experiment)))

  tlx_all_df = tlx_read_many(samples_df, threads=30) %>%
    tlx_extract_bait(bait_size=19, bait_region=12e6) %>%
    tlx_mark_dust() %>%
    tlx_mark_offtargets(offtargets_df, offtarget_region=1e5, bait_region=1e4) %>%
    tlx_calc_copynumber(bowtie2_index="genomes/mm10/mm10", max_hits=100, threads=24)
  tlx_df = tlx_all_df %>%
    tlx_remove_rand_chromosomes() %>%
    dplyr::filter(tlx_copynumber==1 & !tlx_duplicated) %>%
    dplyr::group_by(tlx_sample) %>%
    dplyr::filter(dplyr::n()>5000) %>%
    dplyr::ungroup()
  libfactors_df = tlx_df %>%
    tlx_libsizes() %>%
    tlx_libfactors_within(min(library_size)/library_size)
  tlx_viz_df = tlx_df %>%
    dplyr::mutate(tlx_group=ifelse(tlx_is_bait_chrom, "Inter", "Intra")) %>%
    dplyr::filter(!tlx_is_offtarget & !tlx_is_bait_junction & !(alleles==1 & grepl("promoter/enhancer", experiment) & tlx_is_bait_chrom))

  #
  # TLX coverage
  #
  params_viz = macs2_params(extsize=5e4, exttype="symmetrical", llocal=1e7, minqvalue=0.01, effective_size=1.87e9, maxgap=2e5, minlen=1e5)
  tlxcov_viz_df = tlx_viz_df %>%
    tlx_coverage(group="group", extsize=params_viz$extsize, exttype=params_viz$exttype, libfactors_df=libfactors_df, ignore.strand=F)

  rdc2tlxcov_df = rdc_df %>%
    df2ranges(rdc_chrom, rdc_region_start, rdc_region_end) %>%
    innerJoinByOverlaps(tlxcov_viz_df %>% df2ranges(tlxcov_chrom, tlxcov_start, tlxcov_end))
  rdc_filter = rdc_df %>%
    df2ranges(rdc_chrom, rdc_start, rdc_end) %>%
    innerJoinByOverlaps(tlx_viz_df %>% df2ranges(Rname, Junction, Junction)) %>%
    dplyr::filter(!is.na(rdc_gene_strand) & !is.na(rts_pos)) %>%
    dplyr::group_by(rdc_name, tlx_group) %>%
    dplyr::summarise(junctions_count=dplyr::n()) %>%
    dplyr::filter(junctions_count > 50) %>%
    dplyr::distinct(tlx_group, rdc_name)

  #
  # Prepare ggplot2 data
  #
  rdc2tlxcov_ggplot = rdc2tlxcov_df %>%
    dplyr::inner_join(rdc_filter, by=c("tlx_group", "rdc_name")) %>%
    dplyr::mutate(tlxcov_rel_start=tlxcov_start-rts_pos, tlxcov_rel_end=tlxcov_end-rts_pos) %>%
    dplyr::mutate(rdc_name_i=as.numeric(as.factor(rdc_name))) %>%
    dplyr::group_by(rdc_name) %>%
    dplyr::mutate(tlxcov_rel_pileup=tlxcov_pileup/max(tlxcov_pileup, na.rm=T)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(tlxcov_start=tlxcov_rel_start, tlxcov_end=tlxcov_rel_end, tlxcov_pileup=tlxcov_rel_pileup)

  rdc2replication_ggplot = rdc_df %>%
    dplyr::inner_join(rdc_filter %>% dplyr::distinct(rdc_name), by="rdc_name") %>%
    df2ranges(rdc_chrom, rdc_region_start, rdc_region_end) %>%
    innerJoinByOverlaps(replication_ranges) %>%
    dplyr::mutate(replication_start=pmax(rdc_region_start, replication_start), replication_end=pmin(rdc_region_end, replication_end)) %>%
    dplyr::filter(pmax(replication_start, replication_end)>=rdc_region_start & pmin(replication_start, replication_end)<=rdc_region_end) %>%
    dplyr::mutate(replication_rel_start=replication_start-rts_pos, replication_rel_end=replication_end-rts_pos)

  rdc2genes_ggplot = rdc_df %>%
    dplyr::inner_join(rdc_filter %>% dplyr::distinct(rdc_name), by="rdc_name") %>%
    df2ranges(rdc_chrom, rdc_region_start, rdc_region_end) %>%
    innerJoinByOverlaps(genes_ranges) %>%
    dplyr::mutate(gene_start=pmax(rdc_region_start, gene_start), gene_end=pmin(rdc_region_end, gene_end)) %>%
    dplyr::mutate(gene_rel_start=gene_start-rts_pos, gene_rel_end=gene_end-rts_pos) %>%
    dplyr::mutate(segment_start=ifelse(gene_strand=="+", gene_rel_start, gene_rel_end), segment_end=ifelse(gene_strand=="+", gene_rel_end, gene_rel_start))

  plist = rdc2tlxcov_ggplot %>%
    # dplyr::filter(rdc_name=="RDC_166" & tlx_group=="Inter") %>%
    dplyr::group_split(tlx_group, rdc_gene_strand) %>%
    lapply(function(rdc2tlxcov_ggplot.s) {
      key = rdc2tlxcov_ggplot.s %>% dplyr::distinct(tlx_group, rdc_gene_strand)
      rdc_filter.s = rdc2tlxcov_ggplot.s %>% dplyr::distinct(rdc_name, rdc_gene_strand)
      rdc2genes_ggplot.s = rdc2genes_ggplot %>% dplyr::inner_join(rdc_filter.s, by=c("rdc_name", "rdc_gene_strand"))
      rdc2replication_ggplot.s = rdc2replication_ggplot %>% dplyr::inner_join(rdc_filter.s, by=c("rdc_name", "rdc_gene_strand"))

      strand_colors = c("+"="#FF7A7A", "-"="#7AFF7A")
      g = ggplot() +
        geom_vline(xintercept=0, size=0.3, alpha=0.4) +
        geom_tlxcov(rdc2tlxcov_ggplot.s) +
        geom_segment(aes(x=segment_start, xend=segment_end, y=-0.1*gene_cluster_i, yend=-0.1*gene_cluster_i), size=0.3, arrow=grid::arrow(length=unit(4,"pt")), data=rdc2genes_ggplot.s) +
        geom_segment(aes(x=segment_start, xend=segment_start, y=-0.1*gene_cluster_i-0.05, yend=-0.1*gene_cluster_i+0.05), size=0.3, data=rdc2genes_ggplot.s) +
        geom_text(aes(x=(gene_rel_end+gene_rel_start)/2, y=-0.1*(gene_cluster_i+1)+0.05, label=gene_id), data=rdc2genes_ggplot.s %>% dplyr::distinct(rdc_name, .keep_all=T), size=3, vjust=1, hjust=0, color="#FF6500") +
        geom_segment(aes(x=replication_rel_start, xend=replication_rel_end, y=-0.5+y, yend=-0.5+y, color=replication_strand), size=0.1, arrow=grid::arrow(length=unit(4,"pt")), data=rdc2replication_ggplot.s %>% dplyr::mutate(y=ifelse(replication_strand=="+", 0.02, -0.02))) +
        facet_grid(rdc_chrom+rdc_name~., scales="free_y") +
        labs(x="Relative distance from RTS (100kb)", y="", color="Replication direction", fill="Translocation orientation", title=paste0(names(key), ": ", key, collapse="    ")) +
        scale_x_continuous(labels=scales::label_number(accuracy=1, scale=1e-5, suffix="")) +
        scale_fill_manual(values=strand_colors) +
        scale_color_manual(values=strand_colors) +
        coord_cartesian(xlim=c(-1e6, 1e6)) +
        theme_classic(base_size=12) +
        theme(strip.background=element_blank(), strip.text.y=element_text(size=10, angle=0), legend.position="bottom", axis.title.y=ggplot2::element_blank(), legend.key.size=unit(1.0, 'cm'))

      # if(s<length(plist_strands)) {
      #   g = g + theme(axis.title.x=ggplot2::element_blank())
      # }
      g
  })


  pdf("reports/06-rdc_pileup-tlxcov_relative.pdf", width=8.27, height=11.69)
  plist
  dev.off()
}