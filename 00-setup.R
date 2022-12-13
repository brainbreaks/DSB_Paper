install.packages("devtools")
devtools::install_cran(c("readr", "dplyr", "data.table", "fitdistrplus", "rstatix", "bedr", "reshape2", "ggseqlogo", "randomColorR", "foreach"))
devtools::install_cran(c("igraph", "baseline", "smoother", "dbscan"))
devtools::install_cran(c("ggvenn", "randomcoloR", "ggbeeswarm", "ggpmisc", "ggridges", "ggrepel", "units", "gridpattern", "ggpattern"))
devtools::install_cran(c("ggprism", "ggpubr"))
devtools::install_bioc(c("GenomicFeatures", "GenomicRanges", "rtracklayer", "Biostrings", "ComplexHeatmap"))
devtools::install_deps("breaktools/")

#
# Create Python environment for Keras and Tensorflow
# To run Tensorboard (for validating model training) use command
#
# tensorboard --logdir logs/
#
reticulate::virtualenv_create("r-tensorflow", python="/usr/bin/python3.8")
reticulate::py_install("tensorflow-gpu")
tensorflow::install_tensorflow(envname="r-tensorflow")
keras::install_keras(envname="r-tensorflow", tensorflow="gpu")



# Install blast, bowtie2 and dustmasker if you don't have them. The easiest way is to use Anaconda installer. You might need to add Anaconda binaries folder to your PATH environment variable.
# If you are using OSX+Rstudion also add PATH to ~/.Renviron file:
#
# Sys.setenv(PATH=paste0(Sys.getenv("PATH"), ":/path/to/anaconda/bin"))
#
# https://anaconda.org/bioconda/blast
# https://anaconda.org/bioconda/bowtie2
# https://anaconda.org/bioconda/blat

dir.create("reports", recursive=T, showWarnings=F)
dir.create("genomes", recursive=T, showWarnings=F)
dir.create("tmp", recursive=T, showWarnings=F)

#
# Download TLX files from NCBI
# https://drive.google.com/u/0/uc?id=1gUVUePDl89nnYBTb4ZjL03l8NaLdytdB&export=download
#
unzip("data/data.zip", exdir="data/")
dir.create("data/TLX", recursive=T, showWarnings=F)
file.copy(Sys.glob("data/TLX_paper/*"), "data/TLX", overwrite=T, recursive=T)
file.copy(Sys.glob("data/TLX_public/*"), "data/TLX", overwrite=T, recursive=T)
# file.copy(Sys.glob("reports/00-upload_ncbi/TLX/*.tlx"), "data/TLX", overwrite=TRUE)


