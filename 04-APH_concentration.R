library(readr)
library(dplyr)
library(ggplot2)
library(stringr)
library(ggpattern)
devtools::load_all('breaktools/')
source("00-utils.R")

APH_concentration = function()
{
  debug = F
  # group_palette = c("APH 0.2 uM 96h"="#C6DBEF", "APH 0.3 uM 96h"="#6DAACE", "APH 0.4 uM 96h"="#317BA5", "APH 0.6 uM 96h"="#335E9D", "DMSO"="#CCCCCC")
  group_palette = c("APH 0.2 uM 96h"="#C49A6C", "APH 0.3 uM 96h"="#C49A6C", "APH 0.4 uM 96h"="#C49A6C", "APH 0.6 uM 96h"="#C49A6C", "DMSO"="#CCCCCC")
  subset_palette = c("RDC"="#2DBFC4", "Gene"="#F28C69")
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
  replication_df = readr::read_tsv("data/replication_reduced_subsets.tsv")
  replication_ranges = replication_df %>%
    df2ranges(replication_chrom, pmin(replication_start, replication_end), pmax(replication_start, replication_end))

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
  samples_df = tlx_read_samples(annotation_path="data/htgts_samples.tsv", samples_path="data") %>%
    dplyr::filter(grepl("APH concentration", experiment)) %>%
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
    # dplyr::filter(!tlx_is_offtarget) %>%
    dplyr::mutate(tlx_group=paste0(tlx_group, " (", ifelse(tlx_is_bait_chrom, "Intra", "Inter"), ")")) %>%
    # dplyr::group_by(tlx_is_bait_chrom, Rname) %>%
    # dplyr::mutate(dbscan_cluster=dbscan::dbscan(matrix(Junction), minPts=20, eps=100)$cluster) %>%
    # dplyr::mutate(tlx_strand=ifelse(dbscan_cluster==0, tlx_strand, "cluster")) %>%
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
    readr::write_tsv("reports/04-concentration/APH_concentration_normalization.tsv")
  libfactors_centration_df$offtargets %>%
    readr::write_tsv("reports/04-concentration/APH_concentration_normalization_offtargets.tsv")

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
    innerJoinByOverlaps(replication_ranges) %>%
    dplyr::mutate(
      collision=dplyr::case_when(rdc_gene_strand==replication_strand~"co-directional", T~"head-on"),
      replication_overlap_start=pmax(rdc_start, pmin(replication_end, replication_start)),
      replication_overlap_end=pmin(rdc_end, pmax(replication_end, replication_start)),
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
      readr::write_tsv("APH_junctions_count.tsv")


  genes_nonrdc_df = genes_df %>%
    dplyr::filter(gene_length>=min(rdc_df$rdc_length)) %>%
    df2ranges(gene_chrom, gene_start, gene_end) %>%
    leftJoinByOverlaps(rdc_df %>% df2ranges(rdc_chrom, rdc_start, rdc_end)) %>%
    dplyr::filter(is.na(rdc_name)) %>%
    dplyr::select(dplyr::matches("gene_"))

  tlx_count_df = dplyr::bind_rows(
      genes_nonrdc_df %>% dplyr::select(region_chrom=gene_chrom, region_start=gene_start, region_end=gene_end, region_name=gene_id) %>% dplyr::mutate(region_subset="Gene"),
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
    dplyr::filter(groseq_score>=500) %>%
    dplyr::ungroup()

  # tlx_count_df %>%
  #   dplyr::group_by(region_subset) %>%
  #   dplyr::summarise(n=length(unique(region_name)), min_length=min(region_end-region_start))

  pdf("reports/04-concentration/APH_concentration.pdf", width=11.69, height=8.27, paper="a4r")
  #
  # 1a. Number of junctions increase with concentration (same as 1b but as line with standard error)
  #
  tlx_count_sumdf = tlx_count_df %>%
    # dplyr::mutate(tlx_concentration=factor(tlx_concentration, unique(tlx_concentration))) %>%
    dplyr::group_by(region_subset, tlx_translocation, tlx_concentration) %>%
    dplyr::summarise(junctions_norm_rel.se=sd(junctions_norm_rel, na.rm=T)/sqrt(sum(!is.na(junctions_norm_rel))), junctions_norm_rel=mean(junctions_norm_rel, na.rm=T))
  tlx_count_stat = tlx_count_df %>%
    dplyr::group_by(region_subset, tlx_translocation, tlx_concentration) %>%
    dplyr::mutate(s=sd(junctions_norm_rel, na.rm=T)/sqrt(sum(!is.na(junctions_norm_rel))), m=mean(junctions_norm_rel, na.rm=T), y.position=m+s) %>%
    dplyr::group_by(tlx_translocation, tlx_concentration) %>%
    dplyr::mutate(y.position=max(y.position)+0.1) %>%
    dplyr::group_by(tlx_translocation, tlx_concentration, y.position) %>%
    rstatix::t_test(junctions_norm_rel ~ region_subset) %>%
    # rstatix::remove_ns(col="p", signif.cutoff=0.05) %>%
    dplyr::group_by(tlx_translocation, tlx_concentration) %>%
    dplyr::mutate(group1=as.numeric(tlx_concentration)-0.1, group2=as.numeric(tlx_concentration)+0.1) %>%
    dplyr::ungroup()
  ggplot(tlx_count_sumdf) +
      geom_hline(yintercept=0, linetype="dashed") +
      geom_line(aes(y=junctions_norm_rel, x=as.numeric(tlx_concentration), color=region_subset)) +
      geom_errorbar(aes(ymin=junctions_norm_rel-junctions_norm_rel.se, ymax=junctions_norm_rel+junctions_norm_rel.se, x=as.numeric(tlx_concentration), color=region_subset), width=0.1) +
      ggprism::add_pvalue(tlx_count_stat, tip.length=0.005) +
      facet_wrap(~tlx_translocation, scales="free") +
      labs(y="Offtarget−normalized translocations count\ndevided by DMSO translocations count (per each RDC), log2", color="Subset") +
      scale_x_continuous(breaks=1:nlevels(tlx_count_sumdf$tlx_concentration), labels=levels(tlx_count_sumdf$tlx_concentration)) +
      scale_color_manual(values=subset_palette) +
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

  # x = tlx_proportion_df %>%
  #   # dplyr::mutate(n=ifelse(replication_strand=="+", n_anti, n_sense)) %>%
  #   # dplyr::mutate(n=ifelse(replication_strand=="+", n_sense, n_anti)) %>%
  #   reshape2::dcast(tlx_group+tlx_concentration+tlx_translocation+rdc_name+rdc_gene_strand~collision, value.var="n") %>%
  #   # dplyr::mutate(`co-directional`=tidyr::replace_na(`co-directional`, 1), `head-on`=tidyr::replace_na(`head-on`, 1)) %>%
  #   dplyr::mutate(prop=log2(`co-directional`/`head-on`))
  # ggplot(x, aes(y=prop, x=tlx_concentration)) +
  #   geom_hline(yintercept=0, linetype="dashed") +
  #   geom_boxplot(aes(fill=rdc_gene_strand), outlier.shape=NA, outlier.alpha=0, pattern_fill="#000000", pattern_color="#00000000", pattern_spacing=0.01) +
  #   # geom_point(aes(x=tlx_concentration), position=position_jitterdodge(jitter.width=0.1, jitter.height=0), shape=1, alpha=0.8) +
  #   labs(pattern="Gene strand", y="Proportion of right-moving fork breaks, cent/(cent+tel)\nMost genes will have a head-on area and co-directional area", size="No. of breaks", fill="APH concentration") +
  #   # scale_fill_manual(values=group_palette, guide=guide_legend(override.aes=list(pattern="none"))) +
  #   facet_grid(~tlx_translocation) +
  #   scale_y_continuous(labels=scales::percent, n.breaks=10) +
  #   theme_paper(base_size=12) +
  #   theme_x_factors(size=10) +
  #   theme(legend.key.size=unit(1.2, 'cm'))


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
    tlxcov_concentration_df %>% tlxcov_write_bedgraph(path="reports/04-concentration-bedgraph", group="group", ignore.strand=F)
  }
}


APH_concentration_crosscorrelation_unused = function()
{
  rdc2tlxcov_df = rdc_df %>%
    df2ranges(rdc_chrom, rdc_extended_start, rdc_extended_end) %>%
    innerJoinByOverlaps(tlxcov_concentration_df %>% df2ranges(tlxcov_chrom, tlxcov_start, tlxcov_end)) %>%
    dplyr::mutate(tlx_translocation=gsub(".*(Inter|Intra).*", "\\1", tlx_group), tlx_concentration=gsub("(.*96h|DMSO) .*", "\\1", tlx_group))

  #
  # Bootstrap this specific dataset to further fileter significant RDC
  #
  rdc_bootstrap_df = rdc_df %>%
    tidyr::crossing(rdc2tlxcov_df %>% dplyr::distinct(tlx_group, tlx_translocation, tlx_concentration)) %>%
    dplyr::group_by(tlx_group, tlx_translocation, tlx_concentration) %>%
    dplyr::do((function(z){
      zz<<-z
      bootstrap_data_overlaps(
        evaluated_ranges=z %>% df2ranges(rdc_chrom, rdc_extended_start, rdc_extended_end),
        data_ranges=tlx_concentration_df %>% dplyr::filter(tlx_group==z$tlx_group[1]) %>% df2ranges(Rname, Junction, Junction),
        genome_tiles_step=10000, genome_tiles_width=50000, n_samples=1000)
    })(.)) %>%
    dplyr::ungroup()
  rdc_bootstrap_sumdf = rdc_df %>%
    dplyr::inner_join(rdc_bootstrap_df, by=c("rdc_chrom"="bootstrap_chrom", "rdc_extended_start"="bootstrap_start", "rdc_extended_end"="bootstrap_end")) %>%
    # tidyr::crossing(tlxcov_concentration_df %>% dplyr::distinct(tlx_group)) %>%
    dplyr::group_by(tlx_group, tlx_translocation, tlx_concentration, rdc_name) %>%
    dplyr::filter(!all(bootstrap_data_count==0)) %>%
    dplyr::do((function(z){
      sg = z$bootstrap_data_count[z$bootstrap_type=="signal"]
      fit_gamma = fitdistrplus::fitdist(z$bootstrap_data_count[z$bootstrap_type=="background"], distr="nbinom", method="qme", probs=c(0.1, 0.9))
      fc = sg/fit_gamma$estimate["mu"]
      pvalue = pnbinom(sg, mu=fit_gamma$estimate["mu"], size=fit_gamma$estimate["size"], lower.tail=F)
      data.frame(concentration_bootstrap_pvalue=pvalue, concentration_bootstrap_fc=fc)
    })(.)) %>%
    dplyr::ungroup()

  #
  # Describe pileup to select pileups with enough data
  #
  rdc_peaks_sumdf = rdc2tlxcov_df %>%
    dplyr::group_by(tlx_group, tlx_translocation, tlx_concentration, rdc_name) %>%
    dplyr::do((function(z){
      z %>%
        dplyr::filter(tlxcov_pileup>=4 & tlxcov_pileup>=max(0.5*tlxcov_pileup)) %>%
        df2ranges(tlxcov_chrom, tlxcov_start, tlxcov_end, tlx_strand) %>%
        GenomicRanges::reduce(min.gapwidth=50e3) %>%
        as.data.frame()
    })(.)) %>%
    dplyr::group_by(tlx_group, tlx_translocation, tlx_concentration, rdc_name) %>%
    dplyr::summarise(peaks_telomeric_count=sum(strand=="+"), peaks_centromeric_count=sum(strand=="-"), peaks_telomeric_length=sum((end-start)[strand=="+"]), peaks_centromeric_length=sum((end-start)[strand=="-"]))

  rdc_concentration_df = rdc_df %>%
    dplyr::select(-dplyr::matches("^(rdc_bootstrap_|peaks_)")) %>%
    inner_join(rdc_bootstrap_sumdf %>% dplyr::select(tlx_group, rdc_name, concentration_bootstrap_pvalue), by="rdc_name") %>%
    inner_join(rdc_peaks_sumdf, by=c("tlx_group", "rdc_name")) %>%
    dplyr::mutate(concentration_is_significant=concentration_bootstrap_pvalue<=0.01) %>%
    dplyr::mutate(concentration_has_single_peak=peaks_telomeric_count==1 & peaks_centromeric_count==1 & pmin(peaks_telomeric_length, peaks_centromeric_length)>=50e3)
  rdc2tlxcov_df = rdc_df %>%
    df2ranges(rdc_chrom, rdc_extended_start-250e3, rdc_extended_end+250e3) %>%
    innerJoinByOverlaps(tlxcov_concentration_df %>% df2ranges(tlxcov_chrom, tlxcov_start, tlxcov_end)) %>%
    dplyr::mutate(tlx_translocation=gsub(".*(Inter|Intra).*", "\\1", tlx_group), tlx_concentration=gsub("(.*96h|DMSO) .*", "\\1", tlx_group))
  rdc2tlxcov_concentration_df = rdc2tlxcov_df %>%
    dplyr::inner_join(rdc_concentration_df %>% dplyr::select(tlx_group, rdc_name, dplyr::starts_with("concentration_")), by=c("tlx_group", "rdc_name")) %>%
    dplyr::filter(concentration_is_significant)

  #
  # Cross-correlation ()
  #
  ccs_df = rdc2tlxcov_concentration_df %>%
    dplyr::filter(tlxcov_pileup>1) %>%
    dplyr::group_by(tlx_group, tlx_concentration, tlx_translocation, rdc_name, rdc_gene_strand, rdc_chrom, rdc_extended_start, rdc_extended_end, concentration_is_significant) %>%
    dplyr::do((function(y){
      yy<<-y
      r = tlx_strand_crosscorrelation(tlxcov_df=y, step=1000, negative_lag=T, negative_correlation=F)
      r$pileup_sense = max(y$tlxcov_pileup[y$tlx_strand=="+"], na.rm=T)
      r$pileup_anti = max(y$tlxcov_pileup[y$tlx_strand=="-"], na.rm=T)
      r
    })(.)) %>%
    dplyr::ungroup()
  ccs_good_df = ccs_df %>%
    dplyr::filter(!is.na(crosscorrelation_lag) & crosscorrelation_cor>0.1 & crosscorrelation_cor_pvalue<=0.01 & pmin(pileup_anti, pileup_sense)>=5)

  #
  # 5. Cross-correlation
  #
  translocation_palette = c("Inter"="#00D1FF", "Intra"="#E40006")
  ggplot(ccs_good_df, aes(y=crosscorrelation_lag, x=tlx_concentration)) +
    geom_hline(yintercept=0, linetype="dashed") +
    geom_boxplot(aes(fill=tlx_concentration), outlier.shape = NA, outlier.alpha=0) +
    geom_point(aes(x=tlx_concentration, size=pileup_sense+pileup_anti), color="#FFFFFF", shape=16, position=position_jitter(width=0.2, height=0, seed=2), alpha=0.9) +
    geom_point(aes(x=tlx_concentration, size=pileup_sense+pileup_anti),  shape=21, position=position_jitter(width=0.2, height=0, seed=2), alpha=0.9) +
    labs(y="Relative distance between centromeric(+) and telomeric(-)\nstrand breaks peaks (cross-correlation)", size="Peak height", fill="APH concentration") +
    scale_fill_manual(values=group_palette) +
    scale_color_manual(values=translocation_palette) +
    facet_wrap(~tlx_translocation, scales="free") +
    theme_bw(base_size=12) +
    theme(axis.title.x=element_blank(), axis.text.x=element_text(angle=45, vjust=1, hjust=1, size=10), legend.key.size=unit(1.2, 'cm'))
}