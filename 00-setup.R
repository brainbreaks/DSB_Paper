install.packages("devtools")
devtools::install_cran(c("igraph", "baseline", "smoother", "dbscan"))
devtools::install_cran(c("ggvenn", "randomcoloR", "ggbeeswarm", "ggpmisc", "ggridges", "ggrepel"))
devtools::install_bioc(c("GenomicFeatures", "GenomicRanges", "rtracklayer", "Biostrings", "ComplexHeatmap"))
devtools::install_deps("breaktools/")

unzip("data/TLX.zip", exdir="data")

dir.create("reports", recursive=T, showWarnings=F)
dir.create("genomes", recursive=T, showWarnings=F)
dir.create("tmp", recursive=T, showWarnings=F)

# Download genomes
file.remove("download.py")
download.file("https://raw.githubusercontent.com/brainbreaks/genome_downloader/master/download.py", "download.py")
sys::exec_wait("python3", args=c("download.py", "mm10", "genomes"))

# Required binaries:
# cmake
# gcc
# blat
# bowtie
# samtools
