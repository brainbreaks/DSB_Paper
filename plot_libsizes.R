setwd("~/Workspace/DSB_Paper")

library(readr)
library(dplyr)
library(ggplot2)
library(cowplot)
library(stringr)
# library(plyranges)
devtools::load_all('~/Workspace/breaktools/')


plot_libsizes_all = function()
{
  samples_df = readr::read_tsv("~/Workspace/Datasets/HTGTS/samples/All_samples.tsv", comment="#") %>%
    dplyr::mutate(library_set=dplyr::case_when(grepl("promoter/enhancer", experiment) ~ "Promoter/enhancer inhibition", T~experiment)) %>%
    dplyr::mutate(path=paste0("~/Workspace/Datasets/HTGTS/", path)) %>%
    dplyr::filter(file.exists(path) & !grepl("Tena|Wei",  experiment))

  tlx_df = tlx_read_many(samples_df, threads=30)
  tlx_df = tlx_remove_rand_chromosomes(tlx_df)

  pdf("reports/libsizes.pdf", width=11.69, height=8.27, paper="a4r")
  libsizes_df = tlx_df %>%
    dplyr::mutate(treatment=ifelse(tlx_control, "DMSO", "Aph")) %>%
    dplyr::group_by(library_set, celltype, treatment, organism, .drop=F) %>%
    dplyr::summarize(library_size=dplyr::n())
  ggplot(libsizes_df) +
    geom_bar(aes(x=library_set, y=library_size, fill=treatment), position="dodge", color="#EEEEEE", size=0.1, stat="identity") +
    # geom_bar(aes(x=tlx_sample, y=usefull_size, fill=tlx_control, group=tlx_sample), position="dodge", color="#EEEEEE", size=0.1, stat="identity") +
    # geom_text(aes(x=tlx_sample, y=0, group=tlx_sample, label=tlx_sample), hjust=1, position=position_dodge(width=0.9), size=4, angle=90) +
    # scale_fill_manual(values=color_scheme) +
    scale_y_continuous(labels=scales::label_number(accuracy=1, scale=1e-3, suffix="K")) +
    labs(y="Junctions", title="All TLX files") +
    theme_grey(base_size=14) +
    guides(fill=guide_legend(nrow=1, byrow=TRUE)) +
    ggpubr::theme_pubclean(base_size=12) +
    theme(legend.position="bottom", axis.title.x=element_blank(), axis.text.x=element_text(angle=90, hjust=1), axis.ticks.x=element_blank())
  dev.off()
}

plot_libsizes = function()
{
  # Load TLX
  samples_df = tlx_read_samples("~/Workspace/Datasets/HTGTS/samples/All_samples.tsv", "~/Workspace/Datasets/HTGTS") %>%
    dplyr::filter(!grepl("Sonic hedgehog|Hydroxyurea|Nocodazole", experiment)) %>%
    dplyr::mutate(library_set=dplyr::case_when(grepl("promoter/enhancer", experiment) ~ "Promoter/enhancer inhibition", T~experiment)) %>%
    dplyr::mutate(group=dplyr::case_when(
      control ~ "DMSO",
      grepl("Wei", experiment) & !control ~ "Treatment",
      grepl("Wei", experiment) & control ~ "DMSO",
      grepl("All", experiment) & control ~ "DMSO",
      T ~ gsub(" ?\\(.*", "", group)
    )) %>%
    dplyr::filter(experiment=="APH concentration" & bait_chrom=="chr6")

  samples_df.missing = samples_df %>%
    dplyr::filter(!tlx_exists) %>%
    dplyr::group_by(experiment) %>%
    dplyr::summarise(samples=paste0(sample, collapse=", "))

  writeLines(paste0("TLX files missing: ", paste0(samples_df.missing$experiment, ":\n    ", samples_df.missing$samples, collapse="\n")))

  tlx_df = tlx_read_many(samples_df %>% dplyr::filter(tlx_exists), threads=30)
  tlx_df = tlx_extract_bait(tlx_df, bait_size=19, bait_region=2e6)
  pdf("reports/inter-intra_proportion.pdf", width=8.27, height=11.69)
  interchrom_prop_df = tlx_df %>%
    dplyr::group_by(tlx_sample, experiment) %>%
    dplyr::summarise(prop=sum(tlx_is_bait_chrom)/dplyr::n())
  ggplot(interchrom_prop_df) +
    geom_density(aes(x=prop, fill=experiment), alpha=0.3) +
    labs(x="Intra-/inter- chromosome breaks", y="") +
    theme_bw(base_size=40)
  dev.off()

  # Display library sizes
  libsizes_df = tlx_df %>%
    dplyr::mutate(Junction=dplyr::case_when(tlx_is_bait_junction~"Bait", tlx_is_bait_chrom~"Intra-chromosome", T~"Inter-chromosome")) %>%
    dplyr::mutate(experiment_short=gsub(" \\(.*", "", experiment)) %>%
    dplyr::group_by(tlx_sample, run, bait_chrom, group, experiment, experiment_short, concentration, Junction, .drop=F) %>%
    dplyr::summarize(library_size=dplyr::n(), usefull_size=dplyr::n())

  plist = lapply(split(libsizes_df, f=libsizes_df$experiment_short), FUN=function(df) {
    # df = libsizes_df %>% dplyr::filter(grepl("concentr", experiment))
    dff <<- df
    palette = c(
      "Parental cell"="#FCBBA1", "Allelic deletion"="#EF3B2C", "Allelic+promoter deletion"="#67000D",
      "APH 0.2 uM 96h"="#C6DBEF", "APH 0.3 uM 96h"="#6DAACE", "APH 0.4 uM 96h"="#317BA5", "APH 0.6 uM 96h"="#08306B",
      "Treatment"="#4292C6", "DMSO"="#666666")[unique(df$group)]
    aplha = c("Bait"=0.4, "Intra-chromosome"=1, "Inter-chromosome"=0.8)
    if(any(grepl("concentration", df$experiment))) {
      df = df %>%
        dplyr::arrange(concentration) %>%
        dplyr::mutate(
          tlx_sample=paste0(tlx_sample, " (", bait_chrom, ")"),
          experiment=dplyr::case_when(grepl("experimental", experiment)~"Experimental", T~paste0("APH conc. (",  run, ")"))
        ) %>%
        dplyr::mutate(tlx_sample=factor(tlx_sample, unique(tlx_sample)))
    }

    ggplot(df) +
      geom_bar(aes(x=reorder(tlx_sample, concentration), y=library_size, fill=group, alpha=Junction), position="stack", color="#EEEEEE", size=0.1, stat="identity") +
      labs(y="# junctions", fill="Group") +
      facet_grid(~experiment, scales="free", space="free_x") +
      theme_grey(base_size=10) +
      theme(
        axis.text.x=ggplot2::element_text(angle=45, hjust=1, vjust=1, size=4),
        axis.title.x=ggplot2::element_blank(),
        legend.position=ifelse(!grepl("Csmd1|Ctnna2", df$experiment[1]), "bottom", "none"),
        legend.text=element_text(size=8),
        legend.key.size = unit(8, "pt")
      ) +
      scale_y_continuous(labels=scales::label_number(accuracy=1, scale=1e-3, suffix="K")) +
      scale_fill_manual(values=palette) +
      scale_alpha_manual(values=aplha)
  })

  cowplot::plot_grid(plotlist=plist, align="v", axis="t", ncol=1)


  pdf("reports/libsizes.pdf", width=8.27, height=11.69)
  libsizes_sumdf = tlx_df %>%
    dplyr::mutate(treatment=ifelse(tlx_control, "DMSO", "Aph")) %>%
    dplyr::group_by(library_set, celltype, treatment, organism, .drop=F) %>%
    dplyr::summarize(library_size=dplyr::n())
  ggplot(libsizes_sumdf) +
    geom_bar(aes(x=library_set, y=library_size, fill=treatment), position="dodge", color="#EEEEEE", size=0.1, stat="identity") +
    # geom_bar(aes(x=tlx_sample, y=usefull_size, fill=tlx_control, group=tlx_sample), position="dodge", color="#EEEEEE", size=0.1, stat="identity") +
    # geom_text(aes(x=tlx_sample, y=0, group=tlx_sample, label=tlx_sample), hjust=1, position=position_dodge(width=0.9), size=4, angle=90) +
    # scale_fill_manual(values=color_scheme) +
    scale_y_continuous(labels=scales::label_number(accuracy=1, scale=1e-3, suffix="K")) +
    labs(y="Junctions", title="All TLX files") +
    theme_grey(base_size=14) +
    guides(fill=guide_legend(nrow=1, byrow=TRUE)) +
    ggpubr::theme_pubclean(base_size=12) +
    theme(legend.position="bottom", axis.title.x=element_blank(), axis.text.x=element_text(angle=90, hjust=1), axis.ticks.x=element_blank())

  cowplot::plot_grid(plotlist=plist, align="v", axis="t", ncol=1)
  dev.off()
}