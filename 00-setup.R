install.packages("devtools")
devtools::install_cran(c("igraph", "baseline", "smoother", "dbscan"))
devtools::install_cran(c("ggvenn", "randomcoloR", "ggbeeswarm", "ggpmisc", "ggridges", "ggrepel", "units", "gridpattern", "ggpattern"))
devtools::install_bioc(c("GenomicFeatures", "GenomicRanges", "rtracklayer", "Biostrings", "ComplexHeatmap"))
devtools::install_deps("breaktools/")

unzip("data/TLX.zip", exdir="data")

dir.create("reports", recursive=T, showWarnings=F)
dir.create("genomes", recursive=T, showWarnings=F)
dir.create("tmp", recursive=T, showWarnings=F)
