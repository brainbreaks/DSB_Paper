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
  genes_cache = "tmp/genes.rda"
  if(file.exists(genes_cache)) {
    load(genes_cache)
  } else {
    genes_df = gtf_read('~/Workspace/genomes/mm10/annotation/mm10.ncbiRefSeq.gtf.gz')
    genes_ranges = genes_df %>%
      dplyr::filter(gene_length>=1e5) %>%
      df2ranges(gene_chrom, gene_start, gene_end)
    save(genes_df, genes_ranges, file=genes_cache)
  }

  #
  # Load samples data
  #
  samples_df = tlx_read_samples("~/Workspace/Datasets/HTGTS/samples/All_samples.tsv", "~/Workspace/Datasets/HTGTS") %>%
    dplyr::filter(tlx_exists & celltype=="NPC" & organism=="mouse" & sample!="VI035" & (
      grepl("(Csmd1|Ctnna2|Nrxn1) promoter/enhancer", experiment) |
      grepl("concentration", experiment) & concentration==0.4 |
      grepl("Wei", experiment))
    )

  #
  # Load TLX
  #
  tlx_all_df = tlx_read_many(samples_df, threads=10) %>%
    tlx_extract_bait(bait_size=19, bait_region=12e6) %>%
    tlx_calc_copynumber(bowtie2_index="~/Workspace/genomes/mm10/mm10", max_hits=100, threads=24) %>%
    tlx_mark_offtargets(offtargets_df, offtarget_region=1e5, bait_region=1e4)
  libfactors_df = tlx_all_df %>% tlx_libsizes()

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
    ))
  tlx_rdc_df = dplyr::bind_rows(
    tlx_rdc_df %>% dplyr::mutate(tlx_group=paste0(tlx_group, " (Wei+DKFZ)")),
    tlx_rdc_df %>% dplyr::mutate(tlx_group=paste0(tlx_group, " (", ifelse(grepl("Wei", experiment), "Wei", "DKFZ"), ")"))
  ) %>% dplyr::mutate(tlx_control=F)

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


  # save(tlx_rdc_df, tlxcov_rdc_df, tlxcov_rdc_strand_df, tlxcov_rdc_combined_df, params_rdc, genes_df, libfactors_df, libfactors_df, offtargets_df, bgmodel_df, macs_combined_rdc, params_rdc,  file="backup2.rda")
  # load("backup2.rda")

  #
  # Detect RDC
  #
  params_rdc = macs2_params(extsize=5e4, exttype="symmetrical", llocal=1e7, minpvalue=0.01, maxgap=100e3, minlen=100e3, seedlen=50e3, seedgap=10e3, baseline=2)
  tlxcov_rdc_df = tlx_rdc_df %>% tlx_coverage(group="group", extsize=params_rdc$extsize, exttype=params_rdc$exttype, libfactors_df=libfactors_df, ignore.strand=T, recalculate_duplicate_samples=F)
  tlxcov_rdc_strand_df = tlx_rdc_df %>% tlx_coverage(group="group", extsize=params_rdc$extsize, exttype=params_rdc$exttype, libfactors_df=libfactors_df, ignore.strand=F, recalculate_duplicate_samples=F)
  tlxcov_rdc_combined_df = dplyr::bind_rows(tlxcov_rdc_strand_df, tlxcov_rdc_df)
  bgmodel_combned_df = tlxcov_rdc_combined_df %>%
    dplyr::group_by(tlx_group) %>%
    dplyr::do((function(z) {
      zz<<-z
      z.ranges = df2ranges(z, tlxcov_chrom, tlxcov_start, tlxcov_end, tlx_strand)
      z.mask = coverage_find_empty_intervals(coverage_ranges=z.ranges, coverage_column="tlxcov_pileup", minlen=1e3, mincoverage=0)
      z.bgmodel_df = macs2_coverage_bgmodel(coverage_ranges=z.ranges, distr="nbinom", coverage_column="tlxcov_pileup", mask_ranges=z.mask, debug_plots=F)
      cbind(z[1, "tlx_group", drop=F], z.bgmodel_df)
    })(.))
  macs_combined_rdc = tlxcov_macs2(tlxcov_df=tlxcov_rdc_combined_df, group="group", params=params_rdc, bgmodel_df=bgmodel_combned_df, debug_plots=F, extended_islands=T, extended_islands_dist=1e6, extended_islands_significance=0.1)

  #
  # Find biggest qvalue among all strands for each interval
  #
  qvalues_corrected_df = macs_combined_rdc$qvalues %>%
    dplyr::group_by(tlx_group) %>%
    dplyr::do((function(tlx_g) {
      tlx_g %>%
        df2ranges(qvalue_chrom, qvalue_start, qvalue_end, qvalue_strand) %>%
        coverage_merge_strands(aggregate_fun=max, score_column="qvalue_score") %>%
        as.data.frame() %>%
        dplyr::select(qvalue_chrom=seqnames, qvalue_start=start, qvalue_end=end, qvalue_strand=strand, qvalue_score=score)
    })(.)) %>%
    dplyr::ungroup()

  #
  # Reduce strand-specific and non-specific extended peaks to a single set of detected RDC (stratified by tlx_group)
  #
  rdc_maxgap = 200e3
  rdc_minlen = 100e3
  islands_combined_reduced_df = macs_combined_rdc$islands %>%
    dplyr::group_by(tlx_group) %>%
    dplyr::do(GenomicRanges::reduce(df2ranges(., island_chrom, island_extended_start, island_extended_end), min.gapwidth=rdc_maxgap) %>% as.data.frame()) %>%
    dplyr::ungroup()  %>%
    dplyr::select(island_combined_group=tlx_group, island_combined_chrom=seqnames, island_combined_extended_start=start, island_combined_extended_end=end) %>%
    dplyr::group_by(island_combined_group) %>%
    dplyr::mutate(island_combined_name=paste0("ISLAND_", stringr::str_pad((0:(dplyr::n()))[-1], 3, pad="0"))) %>%
    dplyr::ungroup()

  #
  # Find overlap between stranded and non-stranded peak detection and reduce
  #
  rdc_genes_df = islands_combined_reduced_df %>%
    # dplyr::filter(island_combined_extended_end-island_combined_extended_start>=rdc_minlen) %>%
    df2ranges(island_combined_chrom, island_combined_extended_start, island_combined_extended_end) %>%
    leftJoinByOverlaps(macs_combined_rdc$islands %>% df2ranges(island_chrom, island_extended_start, island_extended_end)) %>%
    dplyr::filter(island_combined_group==tlx_group) %>%
    dplyr::group_by(tlx_group, rdc_chrom=island_combined_chrom, rdc_extended_start=island_combined_extended_start, rdc_extended_end=island_combined_extended_end) %>%
    dplyr::summarize(
      rdc_start=min(island_start),
      rdc_end=max(island_end),
      n=dplyr::n(),
      rdc_significant_nplus=sum(island_strand=="+"),
      rdc_significant_nminus=sum(island_strand=="-"),
      rdc_significant_ncombined=sum(island_strand=="*"),
      rdc_significant_strands=paste0(rdc_significant_nplus, "+", rdc_significant_nminus, "-", rdc_significant_ncombined, "*")) %>%
    df2ranges(rdc_chrom, rdc_extended_start, rdc_extended_end) %>%
    leftJoinByOverlaps(qvalues_corrected_df %>% dplyr::rename(qvalue_tlx_group="tlx_group") %>% df2ranges(qvalue_chrom, qvalue_start, qvalue_end)) %>%
    dplyr::filter(tlx_group==qvalue_tlx_group) %>%
    dplyr::mutate(qvalue_start=pmax(qvalue_start, rdc_extended_start), qvalue_end=pmin(qvalue_end, rdc_extended_end)) %>%
    dplyr::group_by(tlx_group, rdc_chrom, rdc_start, rdc_end) %>%
    dplyr::do((function(z){
      z %>%
        dplyr::filter(qvalue_score>=-log10(params_rdc$minsignif)) %>%
        dplyr::summarize(rdc_sumscore=sum(qvalue_score*(qvalue_end-qvalue_start), na.rm=T), rdc_maxscore=max(qvalue_score, na.rm=T), rdc_significant_sumarea=sum(qvalue_end-qvalue_start, na.rm=T)) %>%
        dplyr::bind_cols(z %>% dplyr::select(dplyr::starts_with("rdc_")) %>% dplyr::slice(1))
    })(.)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(rdc_subset=gsub(".*\\((.*)\\)", "\\1", tlx_group), tlx_group=gsub(" ?\\(.*", "", tlx_group)) %>%
    dplyr::group_by(tlx_group) %>%
    dplyr::mutate(rdc_maxscore.adjusted=-log10(p.adjust(10^-rdc_maxscore))) %>%
    dplyr::ungroup() %>%
    df2ranges(rdc_chrom, rdc_extended_start, rdc_extended_end) %>%
    leftJoinByOverlaps(genes_ranges) %>%
    dplyr::arrange(dplyr::desc(gene_length)) %>%
    dplyr::distinct(tlx_group, rdc_subset, rdc_chrom, rdc_extended_start, rdc_extended_end, .keep_all=T) %>%
    dplyr::mutate(rdc_gene=gene_id, gene_overlap=pmin(rdc_extended_end, gene_end)-pmax(rdc_extended_start, gene_start)) %>%
    dplyr::group_by(tlx_group, rdc_subset) %>%
    dplyr::mutate(rdc_name=paste0("RDC_", stringr::str_pad((0:(dplyr::n()))[-1], 3, pad="0")), rdc_length=rdc_end-rdc_start, rdc_extended_length=rdc_extended_end-rdc_extended_start) %>%
    dplyr::ungroup() %>%
    dplyr::select(tlx_group, dplyr::starts_with("rdc_"))

  #
  # Bootstrap background to find exact probability and signal-to-noise fold change
  #
  genome_tiles_step = 10000
  genome_tiles_width = 50000
  seqlengths = tlx2seqlengths(tlx_rdc_df)
  genome_tiles = seqlengths2tiles(seqlengths, genome_tiles_width, genome_tiles_step)
  tlx_rdc_ranges = tlx_rdc_df %>% df2ranges(Rname, Junction, Junction)
  rdc_bootstrap_df = rdc_genes_df %>%
    dplyr::group_by(tlx_group, rdc_subset) %>%
    dplyr::do((function(z){
      zz<<-z
      # asd()
      z = z %>% dplyr::mutate(rdc_tile_count=ceiling(rdc_extended_length / genome_tiles_width))
      z_ranges = z %>% df2ranges(rdc_chrom, rdc_extended_start, rdc_extended_end)
      z_tiles =  genome_tiles %>% leftJoinByOverlaps(z_ranges) %>% dplyr::filter(is.na(rdc_extended_start)) %>% dplyr::select(dplyr::starts_with("tile_"))
      n_samples = 1000
      bg_df = z_tiles %>%
        dplyr::inner_join(z, by=c("tile_chrom"="rdc_chrom")) %>%
        dplyr::group_by(rdc_name) %>%
        dplyr::sample_n(rdc_tile_count*n_samples, replace=T) %>%
        dplyr::group_by(rdc_name) %>%
        dplyr::mutate(bootstrap_sample_num=floor((1:dplyr::n()-1)/rdc_tile_count)+1) %>%
        dplyr::ungroup()
      tlx_rdc_ranges.f = tlx_rdc_ranges[tlx_rdc_ranges$tlx_group==paste0(z$tlx_group[1], " (", z$rdc_subset[1], ")")]
      bg_df$bootstrap_junctions_count = bg_df %>% df2ranges(tile_chrom, tile_start, tile_end) %>% GenomicRanges::countOverlaps(tlx_rdc_ranges.f)
      bg_df = bg_df %>%
        dplyr::group_by(rdc_name, bootstrap_sample_num) %>%
        dplyr::summarize(bootstrap_junctions_count=sum(bootstrap_junctions_count), bootstrap_type="background") %>%
        dplyr::group_by(rdc_name) %>%
        dplyr::mutate(bootstrap_junctions_mean=mean(bootstrap_junctions_count)) %>%
        dplyr::ungroup()
      sg_df = data.frame(rdc_name=z$rdc_name, bootstrap_sample_num=1, bootstrap_type="signal", bootstrap_junctions_count=z_ranges %>% GenomicRanges::countOverlaps(tlx_rdc_ranges.f)) %>%
        dplyr::mutate(bootstrap_junctions_mean=bootstrap_junctions_count)
      z_res = dplyr::bind_rows(bg_df, sg_df)
    })(.))

  rdc_bootstrap_sumdf = rdc_genes_df %>%
    dplyr::inner_join(rdc_bootstrap_df, by=c("tlx_group", "rdc_subset", "rdc_name")) %>%
    dplyr::group_by(tlx_group, rdc_subset, rdc_name) %>%
    dplyr::do((function(z){
      zz<<-z
      sg = z %>% dplyr::filter(bootstrap_type=="signal") %>% .$bootstrap_junctions_count
      bg = z %>% dplyr::filter(bootstrap_type=="background") %>% .$bootstrap_junctions_count
      # fit_gamma = fitdistrplus::fitdist(bg, distr="norm", method="qme", probs=c(0.1, 0.9))
      # pvalue = pnorm(sg, mean=fit_gamma$estimate["mean"], sd=fit_gamma$estimate["sd"], lower.tail=F)
      # fc = sg/fit_gamma$estimate["mean"]
      # zscore = (sg-fit_gamma$estimate["mean"])/fit_gamma$estimate["sd"]
      fit_gamma = fitdistrplus::fitdist(bg, distr="nbinom", method="qme", probs=c(0.1, 0.9))
      fc = sg/fit_gamma$estimate["mu"]
      pvalue = pnbinom(sg, mu=fit_gamma$estimate["mu"], size=fit_gamma$estimate["size"], lower.tail=F)
      zscore = (sg-fit_gamma$estimate["mu"])/fit_gamma$estimate["size"]
      data.frame(rdc_bootstrap_pvalue=pvalue, rdc_bootstrap_fc=fc, rdc_bootstrap_zscore=zscore)
    })(.))

  if(F) {
    pdf("reports/rdc_bootstrap_boxplots.pdf", width=8.27, height=11.69, paper="a4")
    rdc_genes_df %>%
      dplyr::inner_join(rdc_bootstrap_df, by=c("tlx_group", "rdc_subset", "rdc_name")) %>%
      dplyr::filter(rdc_is_significant) %>%
      dplyr::group_by(tlx_group, rdc_subset) %>%
      dplyr::group_split() %>%
      # dplyr::filter(tlx_group=="APH-Inter" & rdc_subset=="Wei+DKFZ")) %>% %>%
      lapply(function(z)
        ggplot(z) +
          geom_boxplot(aes(x=reorder(rdc_name, bootstrap_junctions_mean/rdc_extended_length), y=1e6*bootstrap_junctions_count/rdc_extended_length, color=bootstrap_type), outlier.shape=NA, position=position_identity()) +
          scale_y_log10(breaks=scales::trans_breaks("log10", function(x) 10^x), labels=scales::trans_format("log10", scales::math_format(10^.x))) +
          annotation_logticks() +
          facet_grid(rdc_chrom~., scales="free_y", space="free_y") +
          labs(x="", y="Junctions (per Mb)", title=paste0(z$tlx_group[1], " (", z$rdc_subset[1], ")"), fill="Bootstrap") +
          guides(fill="none", color="none") +
          theme_bw(base_size=8) +
          theme(axis.text.x=element_text(size=3), axis.title.x=element_text(size=16), strip.text.y.left=element_text(angle=0)) +
          coord_flip()
      )
    dev.off()
  }

  rdc_df = rdc_genes_df %>%
    dplyr::inner_join(rdc_bootstrap_sumdf, by=c("tlx_group", "rdc_subset", "rdc_name")) %>%
    # dplyr::mutate(rdc_is_significant=rdc_maxscore>=2 & rdc_extended_length>=300e3 & rdc_significant_sumarea>=100e3) %>%
    dplyr::mutate(rdc_is_significant=rdc_significant_sumarea>=100e3+params_rdc$extsize & rdc_bootstrap_pvalue<=0.01)
  readr::write_tsv(rdc_df, "data/rdc.tsv")


  #
  # Look at second junctions
  #
  if(F) {
    # baits_df = readr::read_tsv("~/Workspace/Datasets/HTGTS/dkfz_baits.tsv") %>%
    #   dplyr::mutate(cut_expected=ifelse(bait_strand=="+", bait_end, bait_start))
    # tlx_resected_df = tlx_rdc_df %>%
    #   dplyr::filter(bait_name %in% c("Chr17_41Mb", "Chr5_101Mb", "Chr6_70Mb") & tlx_group %in% c("APH-Inter (DKFZ)", "DMSO-Inter (DKFZ)")) %>%
    #   dplyr::filter(Rend-Rstart>=50) %>%
    #   dplyr::inner_join(baits_df, by=c("bait_name"="bait_name")) %>%
    #   dplyr::filter(bait_strand==bait_strand_sgRNA) %>%
    #   dplyr::mutate(cut_effective=ifelse(bait_strand=="+", B_Rend, B_Rstart)) %>%
    #   dplyr::mutate(cut_resection=ifelse(bait_strand_sgRNA=="+", cut_expected-cut_effective, cut_effective-cut_expected))
    #
    # ggplot(tlx_resected_df) +
    #   geom_histogram(aes(x=cut_resection, fill=bait_name), binwidth = 1) +
    #   geom_vline(xintercept=3+0:6*8, linetype="dotted") +
    #   facet_wrap(bait_name~tlx_group, scales="free", ncol=2)
    #
    # tlx_resected_ranges = tlx_resected_df %>%
    #   dplyr::filter(grepl("Chr6", bait_name)) %>%
    #   df2ranges(Rname, Junction, Junction)
    # x = rdc_df %>%
    #   dplyr::filter(rdc_is_significant & tlx_group=="APH-Inter" & rdc_subset=="DKFZ") %>%
    #   df2ranges(rdc_chrom, rdc_extended_start, rdc_extended_end) %>%
    #   innerJoinByOverlaps(tlx_resected_ranges) %>%
    #   dplyr::mutate(cut_resection_raw=cut_resection, cut_resection=cut_resection) %>%
    #   dplyr::group_by(cut_resection) %>%
    #   dplyr::mutate(cut_resection.min=min(cut_resection_raw)) %>%
    #   dplyr::group_by(rdc_chrom, rdc_extended_start, rdc_extended_end, cut_resection, cut_resection.min) %>%
    #   dplyr::summarize(count=dplyr::n()) %>%
    #   dplyr::group_by(rdc_chrom, rdc_extended_start, rdc_extended_end) %>%
    #   dplyr::mutate(count.sum=sum(count, na.rm=T), count.max=max(count, na.rm=T), rank=length(count) - rank(count) + 1) %>%
    #   dplyr::arrange(cut_resection) %>%
    #   dplyr::group_by(rdc_chrom, rdc_extended_start, rdc_extended_end, cut_resection) %>%
    #   dplyr::mutate(prop=count/count.sum) %>%
    #   dplyr::ungroup() %>%
    #   # dplyr::mutate(rdc_id=paste0(rdc_chrom, ":", rdc_start, "-", rdc_end)) %>%
    #   df2ranges(rdc_chrom, rdc_extended_start, rdc_extended_end) %>%
    #   innerJoinByOverlaps(genes_ranges) %>%
    #   dplyr::arrange(dplyr::desc(gene_length)) %>%
    #   dplyr::distinct(rdc_chrom, rdc_extended_start, rdc_extended_end, cut_resection, .keep_all=T) %>%
    #   dplyr::mutate(rdc_id=gene_id)
    # xa = x %>% dplyr::group_by(rdc_id) %>% dplyr::summarize(n=sum(count, na.rm=T)) %>% dplyr::mutate(x=as.numeric(grepl("Grid|Ctnna", rdc_id))) %>% tibble::column_to_rownames("rdc_id")
    # xx = x %>%
    #   reshape2::dcast(rdc_id~cut_resection, value.var="prop") %>%
    #   replace(is.na(.),0) %>%
    #   reshape2::melt(id.vars="rdc_id", variable.name="cut_resection", value.name="value") %>%
    #   dplyr::mutate(cut_resection=as.numeric(as.character(cut_resection))) %>%
    #   dplyr::group_by(rdc_id) %>%
    #   dplyr::mutate(value=zoo::rollmean(value, 5, align = "left", na.pad = T)) %>%
    #   dplyr::ungroup() %>%
    #   dplyr::filter(cut_resection>=0 & cut_resection<=60) %>%
    #   # dplyr::filter(count>=2) %>%
    #   reshape2::dcast(rdc_id~cut_resection, value.var="value") %>%
    #   tibble::column_to_rownames("rdc_id") %>%
    #   replace(is.na(.),0)
    #
    # pheatmap::pheatmap(xx, cluster_cols=F, cluster_rows=T, annotation_row=xa)
  }

  #
  # Write debugging information from RDC calling
  #
  if(debug)
  {
    tlx_rdc_df %>% tlx_write_bed(path="reports/detect_rdc/bed", group="group", ignore.strand=T)
    tlxcov_rdc_combined_df %>% tlxcov_write_bedgraph(path="reports/detect_rdc/bedgraph", group="group", ignore.strand=F)

    macs_combined_rdc$islands %>%
      dplyr::group_by(tlx_group, island_strand) %>%
      dplyr::do((function(df){
        dff<<-df
        df %>%
          dplyr::select(island_chrom, island_extended_start, island_extended_end, island_name, island_score, island_strand, island_start, island_end) %>%
          readr::write_tsv(paste0("reports/detect_rdc/islands-", df$tlx_group[1], "_", df$island_strand[1], ".bed"), col_names=F)
      })(.))

    islands_combined_reduced_df %>%
      dplyr::group_by(island_combined_group) %>%
      dplyr::do((function(df){
        dff<<-df
        # df = islands_combined_reduced_df %>% dplyr::filter(island_combined_group=="APH-Inter (DKFZ)")
        df %>%
          dplyr::mutate(l=island_combined_extended_end-island_combined_extended_start, thickStart=island_combined_extended_start, thickEnd=island_combined_extended_end) %>%
          dplyr::mutate(i=tidyr::replace_na(findInterval(l, seq(rdc_minlen, max(l), length.out=1000))+1, 1), score=1, strand="*", rgb=dplyr::case_when(l<=rdc_minlen~"#BBBBFF", T~colorRampPalette(c("#0000B2", "#FA00FF"))(1000)[i])) %>%
          # dplyr::filter(island_combined_name=="ISLAND_159") %>%
          dplyr::select(island_combined_chrom, island_combined_extended_start, island_combined_extended_end, island_combined_name, score, strand, thickStart, thickEnd, rgb) %>%
          readr::write_tsv(paste0("reports/detect_rdc/islands-reduced-", df$island_combined_group[1], ".bed"), col_names=F)
      })(.))

    rdc_df %>%
      dplyr::filter(rdc_is_significant) %>%
      dplyr::group_by(tlx_group, rdc_subset) %>%
      dplyr::do((function(df){
        dff<<-df
        df %>%
          dplyr::mutate(rdc_strand="*") %>%
          dplyr::select(rdc_chrom, rdc_extended_start, rdc_extended_end, rdc_name, rdc_significant_sumarea, rdc_strand, rdc_start, rdc_end) %>%
          readr::write_tsv(paste0("reports/detect_rdc/rdc-", df$tlx_group[1], " (", df$rdc_subset[1], ").bed"), col_names=F)
      })(.))

    macs_combined_rdc$qvalues %>%
      dplyr::group_by(tlx_group, qvalue_strand) %>%
      dplyr::do((function(df){
        writeLines('track color="255,102,102" altColor="255,0,0"', con=paste0("reports/detect_rdc/qvalues-", df$tlx_group[1], "_", df$qvalue_strand[1], ".bedgraph"))
        df %>%
          dplyr::select(qvalue_chrom, qvalue_start, qvalue_end, qvalue_score) %>%
          readr::write_tsv(paste0("reports/detect_rdc/qvalues-", df$tlx_group[1], "_", df$qvalue_strand[1], ".bedgraph"), col_names=F, append=T)
        writeLines('track color="255,153,51" altColor="255,0,0"', con=paste0("reports/detect_rdc/qvalue-", df$tlx_group[1], ".bedgraph"))
        df %>%
          dplyr::select(qvalue_chrom, qvalue_start, qvalue_end, bgmodel_signal) %>%
          readr::write_tsv(paste0("reports/detect_rdc/baseline-", df$tlx_group[1], "_", df$qvalue_strand[1], ".bedgraph"), col_names=F, append=T)
      })(.))
  }
}