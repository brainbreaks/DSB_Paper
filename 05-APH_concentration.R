library(readr)
library(dplyr)
library(ggplot2)
library(stringr)
library(ggpattern)
devtools::load_all('breaktools/')
source("00-utils.R")

APH_concentration = function()
{
  dir.create("reports/05-APH_concentration", recursive=T, showWarnings=F)

  debug = F
  group_palette = c("APH 0.2 uM 96h"="#C49A6C", "APH 0.3 uM 96h"="#C49A6C", "APH 0.4 uM 96h"="#C49A6C", "APH 0.6 uM 96h"="#C49A6C", "DMSO"="#CCCCCC")
  region_palette = c("RDC"="#00B9D3", "Gene > 100kb"="#D34B00", "Gene < 100kb"="#D39B00")
  params_concentration = macs2_params(extsize=50e3, exttype="symmetrical")

  #
  # Load genes
  #
  genes_df = gtf_read('genomes/mm10/annotation/mm10.ncbiRefSeq.gtf.gz')

  #
  # Load offtargets
  #
  offtargets_df = readr::read_tsv("data/offtargets_dkfz.tsv")

  #
  # Repliseq data
  #
  forks_df = readr::read_tsv("data/replication_forks_NPC.tsv") %>%
    dplyr::mutate(fork_strand=ifelse(fork_direction=="telomeric", "+", "-"))

  #
  # Load GRO-seq data
  #
  groseq_df = data.frame()
  for(f in Sys.glob("data/groseq/bedgraph50k/*.bedgraph")) {
    g = readr::read_tsv(f, col_names=c("groseq_chrom", "groseq_start", "groseq_end", "groseq_score")) %>%
      dplyr::mutate(groseq_filename=basename(f)) %>%
      dplyr::mutate(groseq_strand=gsub(".*(\\+|\\-).bedgraph", "\\1", groseq_filename)) %>%
      dplyr::mutate(groseq_sample=gsub("_bin.*$", "", groseq_filename)) %>%
      dplyr::mutate(groseq_primary=dplyr::case_when(grepl("NXP010", groseq_filename, ignore.case=T)~"NXP010", grepl("NXP047", groseq_filename, ignore.case=T)~"NXP047", T~NA_character_))
    groseq_df = dplyr::bind_rows(groseq_df, g)
  }
  groseq_sumdf = groseq_df %>%
    dplyr::filter(groseq_end-groseq_start==50e3) %>%
    dplyr::filter(groseq_primary=="NXP047") %>%
    dplyr::group_by(groseq_chrom, groseq_start, groseq_end) %>%
    dplyr::summarise(groseq_score=sum(groseq_score))


  #
  # Read TLX
  #
  samples_df = tlx_read_samples("data/htgts_samples.tsv", "data/TLX") %>%
    dplyr::filter(subset_aph_concentration=="Y") %>%
    dplyr::mutate(control=F)

  #
  # Extract additional information about translocations
  #
  tlx_all_df = tlx_read_many(samples_df, threads=24) %>%
    tlx_extract_bait(bait_size=19, bait_region=12e6) %>%
    tlx_calc_copynumber(bowtie2_index="genomes/mm10/mm10", max_hits=100, threads=24) %>%
    tlx_mark_offtargets(offtargets_df, offtarget_region=1e5, bait_region=1e4) %>%
    dplyr::mutate(tlx_sample_raw=tlx_sample, tlx_sample=paste(tlx_group, bait_name), tlx_control=F)

  #
  # Clean-up data
  #
  tlx_concentration_df = tlx_all_df %>%
    tlx_remove_rand_chromosomes() %>%
    dplyr::filter(tlx_copynumber==1 & !tlx_duplicated) %>%
    dplyr::group_by(tlx_sample) %>%
    dplyr::filter(dplyr::n()>5000) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(tlx_group=paste0(tlx_group, " (", ifelse(tlx_is_bait_chrom, "Intra", "Inter"), ")")) %>%
    dplyr::arrange(match(concentration, c(0.2, 0.3,0.4,0.6,0))) %>%
    dplyr::mutate(tlx_translocation=gsub(".*(Inter|Intra).*", "\\1", tlx_group), tlx_concentration=gsub("(.*96h|DMSO) .*", "\\1", tlx_group)) %>%
    dplyr::mutate(tlx_concentration=factor(tlx_concentration, unique(tlx_concentration)))

  #
  # Off-target based normalization calculation
  #
  libfactors_centration_df = tlx_offtarget_libfactor(tlx_all_df, offtargets_df)
  libfactors_centration_df$libfactors %>%
    dplyr::inner_join(tlx_all_df %>% dplyr::distinct(tlx_sample, tlx_sample_raw), by="tlx_sample") %>%
    dplyr::group_by(group=tlx_sample, library_size, library_factor) %>%
    dplyr::summarise(sample_count=length(unique(tlx_sample_raw))) %>%
    readr::write_tsv("reports/05-APH_concentration/APH_concentration_normalization.tsv")
  libfactors_centration_df$offtargets %>%
    readr::write_tsv("reports/05-APH_concentration/APH_concentration_normalization_offtargets.tsv")

  #
  # Load RDC and test to filter out only the RDC that are significant in concentration samples
  #
  rdc_df = readr::read_tsv("data/rdc.tsv") %>%
    dplyr::filter(tlx_group=="APH-Inter" & rdc_subset=="Wei+DKFZ" & rdc_is_significant) %>%
    dplyr::select(-tlx_group, -rdc_subset)

  #
  # Simple shift in mean
  #
  tlx_shift_df = rdc_df %>%
    df2ranges(rdc_chrom, rdc_start, rdc_end) %>%
    innerJoinByOverlaps(tlx_concentration_df %>% df2ranges(Rname, Junction, Junction)) %>%
    dplyr::filter(!tlx_is_offtarget) %>%
    dplyr::ungroup()%>%
    dplyr::group_by(tlx_group, tlx_concentration, tlx_translocation, rdc_name, rdc_gene_strand) %>%
    dplyr::summarize(
      junctions_count=sum(tlx_strand %in% c("+", "-")),
      junctions_sense_count=sum(tlx_strand=="+"),
      junctions_anti_count=sum(tlx_strand=="-"),
      tlx_strand_shift=mean(Junction[tlx_strand=="+"])-mean(Junction[tlx_strand=="-"]),
      tlx_strand_relshift=tlx_strand_shift/(max(rdc_end)-min(rdc_start)),
      tlx_fc=junctions_sense_count/(junctions_count),
      tlx_fc=dplyr::case_when(is.finite(tlx_fc)~tlx_fc, T~NA_real_),
      tlx_relfc=tlx_fc/(max(rdc_end)-min(rdc_start))
    )

  tlx_proportion_df = rdc_df %>%
    # dplyr::filter(!is.na(rdc_gene)) %>%
    # dplyr::inner_join(genes_df, by=c("rdc_gene"="gene_id")) %>%
    # dplyr::mutate(rdc_start=gene_start, rdc_end=gene_end) %>%
    df2ranges(rdc_chrom, rdc_start, rdc_end) %>%
    innerJoinByOverlaps(tlx_concentration_df %>% df2ranges(Rname, Junction, Junction)) %>%
    dplyr::filter(!tlx_is_offtarget) %>%
    df2ranges(rdc_chrom, Junction, Junction) %>%
    innerJoinByOverlaps(forks_df %>% df2ranges(fork_chrom, fork_start, fork_end)) %>%
    dplyr::mutate(
      collision=dplyr::case_when(rdc_gene_strand==fork_strand~"co-directional", T~"head-on"),
      replication_overlap_start=pmax(rdc_start, fork_start),
      replication_overlap_end=pmin(rdc_end, fork_end),
      replication_overlap_length=replication_overlap_end-replication_overlap_start
    ) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(tlx_group, tlx_concentration, tlx_translocation, rdc_name, rdc_gene, rdc_gene_strand, collision) %>%
    dplyr::summarize(
      replication_overlap_length=sum(replication_overlap_length),
      junctions_count=sum(tlx_strand %in% c("+", "-")),
      junctions_sense_count=sum(tlx_strand=="+"),
      junctions_anti_count=sum(tlx_strand=="-"),
      tlx_fc=dplyr::case_when(rdc_gene_strand[1]=="+"~junctions_sense_count/(junctions_sense_count+junctions_anti_count), rdc_gene_strand[1]=="-"~junctions_anti_count/(junctions_sense_count+junctions_anti_count)),
      tlx_relfc=tlx_fc/(max(rdc_end)-min(rdc_start))
    ) %>%
    dplyr::filter(pmin(junctions_sense_count, junctions_anti_count)>=1) %>%
    dplyr::filter(!is.na(rdc_gene_strand)) %>%
    dplyr::ungroup()

  #
  # Export to file
  #
  x = tlx_proportion_df %>%
    dplyr::group_by(tlx_translocation, tlx_concentration, rdc_name, rdc_gene) %>%
    dplyr::summarise(
      n=dplyr::n(),
      fork_length.head_on=sum(replication_overlap_length[collision=="head-on"]),
      fork_length.co_directional=sum(replication_overlap_length[collision=="co-directional"]),
      junctions_all_count.head_on=sum(junctions_count[collision=="head-on"]),
      junctions_all_count.co_directional=sum(junctions_count[collision=="co-directional"]),
      junctions_telomeric_count.head_on=sum(junctions_sense_count[collision=="head-on"]),
      junctions_telomeric_count.co_directional=sum(junctions_sense_count[collision=="co-directional"]),
      junctions_centromeric_count.head_on=sum(junctions_anti_count[collision=="head-on"]),
      junctions_centromeric_count.co_directional=sum(junctions_anti_count[collision=="co-directional"])) %>%
      readr::write_tsv("reports/05-APH_concentration/APH_junctions_count.tsv")


  genes_nonrdc_df = genes_df %>%
    # dplyr::filter(gene_length>=min(rdc_df$rdc_length)) %>%
    df2ranges(gene_chrom, gene_start, gene_end) %>%
    leftJoinByOverlaps(rdc_df %>% df2ranges(rdc_chrom, rdc_start, rdc_end)) %>%
    dplyr::filter(is.na(rdc_name)) %>%
    dplyr::select(dplyr::matches("gene_")) %>%
    dplyr::mutate(region_subset=ifelse(gene_length>=100e3, "Gene > 100kb", "Gene < 100kb"))

  tlx_count_df = dplyr::bind_rows(
      genes_nonrdc_df %>% dplyr::select(region_chrom=gene_chrom, region_start=gene_start, region_end=gene_end, region_name=gene_id, region_subset),
      # tlx_offtargets_df %>% dplyr::distinct(region_chrom=offtarget_chrom, region_start=offtarget_start, region_end=offtarget_end, offtarget_name=offtarget_bait_name) %>% dplyr::mutate(region_subset="Off-target"),
      rdc_df %>% dplyr::select(region_chrom=rdc_chrom, region_start=rdc_start, region_end=rdc_end, region_name=rdc_name) %>% dplyr::mutate(region_subset="RDC")) %>%
    df2ranges(region_chrom, region_start, region_end) %>%
    innerJoinByOverlaps(groseq_sumdf %>% df2ranges(groseq_chrom, groseq_start, groseq_end)) %>%
    dplyr::group_by(region_chrom, region_start, region_end) %>%
    dplyr::mutate(groseq_score=sum(groseq_score, na.rm=T)) %>%
    df2ranges(region_chrom, region_start, region_end) %>%
    innerJoinByOverlaps(tlx_concentration_df %>% df2ranges(Rname, Junction, Junction)) %>%
    dplyr::arrange(match(concentration, c(0, 0.2, 0.3,0.4,0.6))) %>%
    dplyr::mutate(tlx_concentration=factor(tlx_concentration, unique(tlx_concentration))) %>%
    dplyr::inner_join(libfactors_centration_df$libfactors, by="tlx_sample") %>%
    dplyr::group_by(tlx_translocation, tlx_concentration, concentration, region_subset, region_name, region_chrom, region_start, region_end, groseq_score) %>%
    dplyr::summarise(junctions_count=dplyr::n(), junctions_norm_count=sum(library_factor)) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(tlx_translocation, region_subset, region_name, region_chrom, region_start, region_end, groseq_score) %>%
    dplyr::do((function(z){
      zz<<-z

      ret = data.frame(
        junctions_count=z$junctions_count,
        tlx_concentration=z$tlx_concentration,
        junctions_norm_count=z$junctions_norm_count,
        junctions_rel=z$junctions_count/mean(z$junctions_count),
        junctions_norm_mean=z$junctions_norm_count/mean(z$junctions_norm_count)
      )
      ref_concentration = "DMSO"
      ref_junctions_norm_count = ret$junctions_norm_count[ret$tlx_concentration==ref_concentration]

      if(any(ret$tlx_concentration==ref_concentration)) {
        ret$junctions_norm_rel = log2(ret$junctions_norm_count/ref_junctions_norm_count)
      } else {
        ret$junctions_norm_rel = NA_real_
      }
      ret
    })(.)) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(tlx_translocation, region_subset, region_name, region_chrom, region_start, region_end, groseq_score) %>%
    dplyr::filter(groseq_score>=500 & region_subset!="Gene < 100kb") %>%
    dplyr::ungroup()

  # tlx_count_df %>%
  #   dplyr::group_by(region_subset) %>%
  #   dplyr::summarise(n=length(unique(region_name)), min_length=min(region_end-region_start))

  pdf("reports/05-APH_concentration/APH_concentration.pdf", width=11.69, height=8.27, paper="a4r")
  #
  # 1a. Number of junctions increase with concentration (same as 1b but as line with standard error)
  #
  tlx_count_sumdf = tlx_count_df %>%
    dplyr::group_by(region_subset, tlx_translocation, tlx_concentration) %>%
    dplyr::summarise(junctions_norm_rel.se=sd(junctions_norm_rel, na.rm=T)/sqrt(sum(!is.na(junctions_norm_rel))), junctions_norm_rel=mean(junctions_norm_rel, na.rm=T))

  tlx_count_stat = tlx_count_df %>%
    dplyr::group_by(tlx_translocation, tlx_concentration) %>%
    rstatix::t_test(junctions_norm_rel ~ region_subset) %>%
    dplyr::filter((group1 == "Gene > 100kb" | group2 == "Gene > 100kb") & tlx_concentration != "DMSO") %>%
    dplyr::inner_join(tlx_count_sumdf %>% dplyr::rename(junctions_norm_rel.se.g1="junctions_norm_rel.se", junctions_norm_rel.g1="junctions_norm_rel"), by=c("tlx_translocation", "tlx_concentration", "group1"="region_subset")) %>%
    dplyr::inner_join(tlx_count_sumdf %>% dplyr::rename(junctions_norm_rel.se.g2="junctions_norm_rel.se", junctions_norm_rel.g2="junctions_norm_rel"), by=c("tlx_translocation", "tlx_concentration", "group2"="region_subset")) %>%
    dplyr::mutate(group1=as.factor(group1), group2=as.factor(group2)) %>%
    dplyr::group_by(tlx_translocation, tlx_concentration, group1, group2) %>%
    dplyr::mutate(
      ymin=pmin(junctions_norm_rel.g1, junctions_norm_rel.g2),
      ymax=pmax(junctions_norm_rel.g1, junctions_norm_rel.g2),
      x=(as.numeric(group1)+as.numeric(group2))[1]/2,
      p.signif=dplyr::case_when(p<=0.0001~"****", p<=0.001~"***", p<=0.01~"**", p<=0.05~"*", T~"ns")
    )

  ggplot(tlx_count_sumdf) +
      facet_wrap(~tlx_translocation, scales="free") +
      geom_line(aes(y=junctions_norm_rel, x=as.numeric(tlx_concentration), color=region_subset)) +
      geom_errorbar(aes(ymin=junctions_norm_rel-junctions_norm_rel.se, ymax=junctions_norm_rel+junctions_norm_rel.se, x=as.numeric(tlx_concentration), color=region_subset), width=0.1) +
      geom_hline(yintercept=0, linetype="dashed") +
      # ggprism::add_pvalue(tlx_count_stat, tip.length=0.005) +
      geom_segment(aes(x=as.numeric(tlx_concentration)+x/6,y=ymin,xend=as.numeric(tlx_concentration)+x/6,yend=ymax), data=tlx_count_stat, size=0.2) +
      geom_segment(aes(x=as.numeric(tlx_concentration)+x/6,y=ymin,xend=as.numeric(tlx_concentration)+x/6-0.05,yend=ymin), data=tlx_count_stat, size=0.2) +
      geom_segment(aes(x=as.numeric(tlx_concentration)+x/6,y=ymax,xend=as.numeric(tlx_concentration)+x/6-0.05,yend=ymax), data=tlx_count_stat, size=0.2) +
      geom_text(aes(x=as.numeric(tlx_concentration)+x/6, y=ymax/2+ymin/2, label=p.signif), data=tlx_count_stat, hjust=-0.1) +
      labs(y="Offtarget−normalized translocations count\ndevided by DMSO translocations count (per each RDC), log2", color="Subset") +
      scale_x_continuous(breaks=1:nlevels(tlx_count_sumdf$tlx_concentration), labels=levels(tlx_count_sumdf$tlx_concentration)) +
      scale_color_manual(values=region_palette) +
      theme_paper(base_size=12) +
      theme_x_factors(size=10) +
      theme(legend.key.size=unit(1.2, 'cm'))


  #
  # 1b. Number of junctions increase with concentration
  #
  tlx_count_rdc_df = tlx_count_df %>% dplyr::filter(region_subset=="RDC") %>% dplyr::mutate(junctions_norm_mean100=junctions_norm_mean*100)
  tlx_count_rdc_stat = tlx_count_rdc_df %>%
    dplyr::group_by(tlx_translocation) %>%
    rstatix::t_test(junctions_norm_mean100 ~ tlx_concentration) %>%
    rstatix::add_xy_position()
  ggplot(tlx_count_rdc_df, aes(y=junctions_norm_mean100, x=tlx_concentration)) +
      geom_hline(yintercept=100, linetype="dashed") +
      geom_boxplot(aes(fill=tlx_concentration), outlier.shape=NA, outlier.alpha=0, pattern_fill="#000000", pattern_color="#00000000", pattern_spacing=0.01) +
      geom_point(aes(fill=region_subset), position=position_jitterdodge(jitter.width=0.1, jitter.height=0, seed=2), alpha=0.2) +
      ggprism::add_pvalue(tlx_count_rdc_stat, tip.length=0.005) +
      labs(y="Offtarget-normalized translocations count\ndevided by mean translocations count (per each RDC), log2", fill="APH concentration", pattern="Junction location") +
      scale_fill_manual(values=group_palette, guide=guide_legend(override.aes=list(pattern="none"))) +
      facet_wrap(~tlx_translocation) +
      scale_y_continuous(breaks=function(y) seq(0, ceiling(max(y)), 50)) +
      theme_paper(base_size=12) +
      theme_x_factors(size=10) +
      theme(legend.key.size=unit(1.2, 'cm'))

  #
  # 2. Plot shift of junctions peak at each APH each concentration
  #
  tlx_relshift_df = tlx_shift_df %>% dplyr::filter(is.finite(tlx_strand_shift))
  tlx_relshift_stat = tlx_relshift_df %>%
    dplyr::group_by(tlx_translocation) %>%
    rstatix::t_test(tlx_strand_relshift ~ tlx_concentration) %>%
    rstatix::add_xy_position(scales="free")
  ggplot(tlx_relshift_df, aes(y=tlx_strand_relshift, x=tlx_concentration)) +
    geom_hline(yintercept=0, linetype="dashed") +
    geom_boxplot(aes(fill=tlx_concentration), outlier.shape = NA, outlier.alpha = 0) +
    geom_point(aes(size=junctions_sense_count+junctions_anti_count), position=position_jitter(width=0.2, height=0, seed=2), alpha=0.3) +
    # geom_text(aes(label=gsub("RDC_", "", rdc_name)), color="#FF0000", position=position_jitter(width=0.2, height=0, seed=2)) +
    ggprism::add_pvalue(tlx_relshift_stat, tip.length=0.005) +
    labs(y="Relative distance between centromeric(+) and telomeric(-)\n junctions average locations", size="No. of breaks", fill="APH concentration", pattern="Junction location") +
    scale_fill_manual(values=group_palette) +
    facet_wrap(~tlx_translocation, scales="free") +
    scale_y_continuous(labels=scales::percent, n.breaks=20) +
    theme_paper(base_size=12) +
    theme_x_factors(size=10) +
    theme(legend.key.size=unit(1.2, 'cm'))

  #
  # 3. Plot proportion of junctions in each concentration
  #
  tlx_shift_df %>%
    dplyr::filter(!is.na(rdc_gene_strand)) %>%
    reshape2::melt(measure.vars=c("junctions_sense_count", "junctions_anti_count")) %>%
    dplyr::group_by(tlx_translocation, rdc_name) %>%
    dplyr::mutate(junctions_count=(value)/sum(value)) %>%
    dplyr::mutate(variable=ifelse(variable=="junctions_sense_count", "+", "-")) %>%
    ggplot(aes(y=junctions_count, x=tlx_concentration)) +
      ggpattern::geom_boxplot_pattern(aes(fill=tlx_concentration, pattern=variable), outlier.shape=NA, outlier.alpha=0, pattern_fill="#000000", pattern_color="#00000000", pattern_spacing=0.01) +
      geom_point(aes(fill=paste0(tlx_concentration, variable)), position=position_jitterdodge(jitter.width=0.2, jitter.height=0, seed=2), alpha=0.3) +
      labs(y="Proportion of RDC breaks per APH concentration (each point is RDC)", size="No. of breaks", fill="APH concentration", pattern="Junction orientation") +
      scale_fill_manual(values=group_palette, guide=guide_legend(override.aes=list(pattern="none"))) +
      ggpattern::scale_pattern_manual(values=c("-"='none', "+"='stripe'), guide=guide_legend(override.aes=list(fill=group_palette[3], colour="#00000000"))) +
      facet_grid(~tlx_translocation) +
      scale_y_continuous(labels=scales::percent, n.breaks=20) +
      theme_paper(base_size=12) +
      theme_x_factors(size=10) +
      theme(legend.key.size=unit(1.2, 'cm'))

  #
  # 4. Plot proportion of junctions in each concentration
  #
  tlx_proportion_stat = tlx_proportion_df %>%
    dplyr::group_by(tlx_translocation, tlx_concentration) %>%
    rstatix::t_test(tlx_fc ~ collision) %>%
    dplyr::mutate(p.signif=dplyr::case_when(p<=0.0001~"****", p<=0.001~"***", p<=0.01~"**", p<=0.05~"*", T~"ns")) %>%
    rstatix::add_xy_position() %>%
    dplyr::mutate(group1=tlx_concentration, group2=tlx_concentration, xmin=as.numeric(tlx_concentration)-0.2, xmax=as.numeric(tlx_concentration)+0.2)
  ggplot(tlx_proportion_df, aes(y=tlx_fc, x=tlx_concentration)) +
    geom_hline(yintercept=0.5, linetype="dashed") +
    ggpattern::geom_boxplot_pattern(aes(fill=tlx_concentration, pattern=collision), outlier.shape=NA, outlier.alpha=0, pattern_fill="#000000", pattern_color="#00000000", pattern_spacing=0.01) +
    # geom_point(aes(x=tlx_concentration, color=rdc_gene_strand), position=position_jitterdodge(jitter.width=0.1, jitter.height=0), shape=1, alpha=0.8) +
    ggprism::add_pvalue(tlx_proportion_stat, xmin="xmin", xmax="xmax", tip.length=0.005) +
    labs(pattern="Region", y="Proportion of translocations in alignment with replication fork\n(aligned_translocations/all_translocations)", fill="APH concentration") +
    ggpattern::scale_pattern_manual(values=c("co-directional"='none', "head-on"='stripe'), guide=guide_legend(override.aes=list(fill=group_palette[3], colour="#00000000"))) +
    scale_fill_manual(values=group_palette, guide=guide_legend(override.aes=list(pattern="none"))) +
    scale_color_manual(values=c("+"="#000000", "-"="#000000"), guide=guide_legend(override.aes=list(pattern="none"))) +
    scale_y_continuous(labels=scales::percent, n.breaks=10) +
    guides(color="none") +
    facet_grid(~tlx_translocation) +
    theme_paper(base_size=12) +
    theme_x_factors(size=10) +
    theme(legend.key.size=unit(1.2, 'cm'))
  dev.off()


  #
  # Calculate coverage and export bedgraph
  #
  if(debug) {
    tlxcov_concentration_df = tlx_concentration_df %>%
      tlx_coverage(group="group", extsize=params_concentration$extsize, exttype=params_concentration$exttype, libfactors_df=libfactors_centration_df$libfactors, ignore.strand=F)
    tlxcov_concentration_df %>% tlxcov_write_bedgraph(path="reports/05-APH_concentration/bedgraph", group="group", ignore.strand=F)
  }
}