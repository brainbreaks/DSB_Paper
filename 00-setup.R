install.packages("devtools")
devtools::install_cran(c("readr", "dplyr", "data.table", "fitdistrplus", "rstatix", "bedr", "reshape2", "ggseqlogo", "randomColorR", "foreach"))
devtools::install_cran(c("igraph", "baseline", "smoother", "dbscan"))
devtools::install_cran(c("ggvenn", "randomcoloR", "ggbeeswarm", "ggpmisc", "ggridges", "ggrepel", "units", "gridpattern", "ggpattern"))
devtools::install_cran(c("ggprism", "ggpubr"))
devtools::install_bioc(c("GenomicFeatures", "GenomicRanges", "rtracklayer", "Biostrings", "ComplexHeatmap"))
devtools::install_deps("breaktools/")

# Install blast, bowtie2 and dustmasker if you don't have them. The easiest way is to use Anaconda installer. You might need to add Anaconda binaries folder to your PATH environment variable.
# If you are using OSX+Rstudion also add PATH to ~/.Renviron file:
#
# Sys.setenv(PATH=paste0(Sys.getenv("PATH"), ":/path/to/anaconda/bin"))
#
# https://anaconda.org/bioconda/blast
# https://anaconda.org/bioconda/bowtie2
# https://anaconda.org/bioconda/blat

unzip("data/TLX.zip", exdir="data")

dir.create("reports", recursive=T, showWarnings=F)
dir.create("genomes", recursive=T, showWarnings=F)
dir.create("tmp", recursive=T, showWarnings=F)
